/******************************************************************************/
/**
  @file vmc.cpp

  Solve the many-body problem for the given system nucleus.

  Syntax:
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vmc number_samples prng_seed number_protons number_neutrons particle_mass
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  @author Patrick J. Fasano
  University of Notre Dame

  + 12/5/16 (pjf): Created, based on radialutils from shell.

*******************************************************************************/

#include <mkl_vsl.h>
#include <omp.h>
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "eigen3/Eigen/Core"

#include "units.h"
#include "random_buffer.h"
#include "wavefunction.h"
#include "potential.h"


////////////////////////////////////////////////////////////////
// tunable parameters
/////////////////////////////////////////////////////////////////

const std::vector<double>::size_type kRandomCacheSize = 3000;
const double kStepSize = 1*units::kFM;
const int kDimension = 3;
const double kTargetAcceptance = 0.4;

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

// Stores simple parameters for run
struct RunParameters {
  // mode
  uint64_t number_samples;
  uint32_t prng_seed;
  int number_protons;
  int number_neutrons;
  double particle_mass;
  double oscillator_length;
};

void PrintUsage(char *argv[]) {
  std::cout << "Usage: " << argv[0]
            << " number_samples prng_seed number_protons number_neutrons particle_mass oscillator_length"
            << std::endl;
  std::cout << "  -- particle mass should be given in proton masses" << std::endl;
}

void ProcessArguments(int argc, char *argv[], RunParameters& run_parameters) {
  // usage message
  if (argc-1 != 6) {
    PrintUsage(argv);
    std::exit(EXIT_FAILURE);
  }

  // number of samples
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> run_parameters.number_samples;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid number of samples" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // random number generator seed
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> run_parameters.prng_seed;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid PRNG seed" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // number of protons
  {
    std::istringstream parameter_stream(argv[3]);
    parameter_stream >> run_parameters.number_protons;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid number of particles" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // number of neutrons
  {
    std::istringstream parameter_stream(argv[4]);
    parameter_stream >> run_parameters.number_neutrons;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid number of particles" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // particle mass
  {
    std::istringstream parameter_stream(argv[5]);
    parameter_stream >> run_parameters.particle_mass;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid particle mass." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // variational parameter
  {
    std::istringstream parameter_stream(argv[6]);
    parameter_stream >> run_parameters.oscillator_length;
    if (!parameter_stream) {
      PrintUsage(argv);
      std::cerr << "Invalid oscillator length." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

int main(int argc, char *argv[]) {
  // initialize MPI and OpenMP
  int num_mpi_ranks, mpi_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // process arguments
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  // initialize PRNG
  VSLStreamStatePtr stream;
  vslNewStream(&stream, VSL_BRNG_MT2203+mpi_rank, run_parameters.prng_seed);

  // initialize random number buffers
  buffer::RandomBuffer<double> random_step_buffer(kDimension*kRandomCacheSize, stream, -kStepSize, kStepSize);
  buffer::RandomBuffer<double> random_uniform_buffer(kRandomCacheSize, stream, 0, 1);

  // initialize single particle wave function and potential
  wf::SphericalOscillatorWF wavefunc(run_parameters.oscillator_length*units::kFM);
  potential::AnharmonicOscillator V(units::kFM);

  // initialize our variables for collecting observables
  double ke = 0, pe = 0;
  int64_t n = 0;

  // initialize Metropolis state
  Eigen::Vector3d state({0., 0., 0.});
  double prev_norm = 1e-6;
  int64_t  moves = 0;
  int64_t proposals = 0;
  double acceptance_rate;

  while (n < run_parameters.number_samples) {
    // main Metropolis algorithm loop (parallelized with OpenMP)
    // #pragma omp parallel for reduction(+:ke) reduction(+:pe) reduction(+:n)
    for (size_t i = 0; i < (kRandomCacheSize-1); ++i) {
      // propose a new position and calculate its norm
      Eigen::Vector3d proposed_state = state +
        Eigen::Vector3d({
          random_step_buffer[i],
          random_step_buffer[i+kRandomCacheSize],
          random_step_buffer[i+2*kRandomCacheSize]
        });
      double proposed_norm = wavefunc.norm(proposed_state);
      ++proposals;

      // if the new state doesn't pass the Metropolis test, continue to the next iteration
      if ((proposed_norm/prev_norm) < 1.) {
        if (random_uniform_buffer[i] > (proposed_norm/prev_norm)) { continue; }
      }
      // otherwise, move the walker
      state = proposed_state;
      prev_norm = proposed_norm;
      ++moves;

      if (moves > 1000) {
        double kinetic_energy = (-0.5 * units::kH_bar2 * units::kM_p) * wavefunc.local_laplacian(state);
        ke += kinetic_energy;
        pe += V(state);
        // pe += (units::kH_bar2 * std::pow(state.dot(state), 2)) / (2 * units::kM_p * std::pow(units::kFM,6));
        ++n;
      }
    }
    acceptance_rate = static_cast<double>(moves)/static_cast<double>(proposals);
    if (acceptance_rate > kTargetAcceptance) {
      random_step_buffer.min((acceptance_rate-kTargetAcceptance+1)*random_step_buffer.min());
      random_step_buffer.max((acceptance_rate-kTargetAcceptance+1)*random_step_buffer.max());
    } else if (acceptance_rate < kTargetAcceptance) {
      random_step_buffer.min((acceptance_rate-kTargetAcceptance+1)*random_step_buffer.min());
      random_step_buffer.max((acceptance_rate-kTargetAcceptance+1)*random_step_buffer.max());

    }
    random_uniform_buffer.Fill();
    random_step_buffer.Fill();
  }
  std::cout << "Kinetic energy:   " << ke / n / units::kMeV << " MeV" << std::endl;
  std::cout << "Potential energy: " << pe / n / units::kMeV << " MeV" << std::endl;
  std::cout << "Total energy:     " << (ke/n + pe/n) / units::kMeV << " MeV" << std::endl;
  std::cout << "Acceptance rate:  " << acceptance_rate * 100 << " %" << std::endl;

  /* code */
  return 0;
}
