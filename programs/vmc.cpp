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


////////////////////////////////////////////////////////////////
// tunable parameters
/////////////////////////////////////////////////////////////////

const std::vector<double>::size_type kRandomCacheSize = 3000;
const double kStepSize = 5*units::kFM;

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
};

void PrintUsage(char *argv[]) {
  std::cout << "Usage: " << argv[0]
            << " number_samples prng_seed number_protons number_neutrons particle_mass"
            << std::endl;
  std::cout << "  -- particle mass should be given in proton masses" << std::endl;
}

void ProcessArguments(int argc, char *argv[], RunParameters& run_parameters) {
  // usage message
  if (argc-1 != 5) {
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
  buffer::RandomBuffer<double> random_double_buffer(3*kRandomCacheSize, stream, -kStepSize, kStepSize);
  // buffer::RandomBuffer<int> random_int_buffer(kRandomCacheSize, stream, 0, run_parameters.number_neutrons);

  // initialize single particle wave function
  wf::SphericalOscillatorWF wavefunc(1*units::kFM);

  // initialize our variables for collecting observables
  double tot = 0, tot2 = 0;
  int n = 0;
  while (n < run_parameters.number_samples) {
    #pragma omp parallel for reduction(+:tot) reduction(+:tot2) reduction(+:n)
    for (size_t i = 0; i < kRandomCacheSize-3; ++n, i+=3) {
      Eigen::Vector3d vec({random_double_buffer[i], random_double_buffer[i+1], random_double_buffer[i+2]});
      double kinetic_energy = (-0.5 * units::kH_bar2 * units::kM_p) * wavefunc.local_laplacian(
        vec
      );
      tot += kinetic_energy;
      tot2 += (units::kH_bar2 * vec.dot(vec)) / (2 * units::kM_p * std::pow(units::kFM,4));
    }
    random_double_buffer.Fill();
  }
  std::cout << tot / n << std::endl;
  std::cout << tot2 / n << std::endl;
  std::cout << tot/n + tot2/n << std::endl;

  /* code */
  return 0;
}
