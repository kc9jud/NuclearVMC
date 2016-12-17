/************************************************************/
/**
  @file wavefunction.h

  VMC trial wavefunctions.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 12/14/16 (pjf): Created.

****************************************************************/
#ifndef WAVEFUNCTION_H_
#define WAVEFUNCTION_H_

#include <vector>
#include <cmath>
#include "eigen3/Eigen/Core"


namespace wf {

class WaveFunction {
 public:
  WaveFunction() {}

  virtual double value(const Eigen::Vector3d& r) { return 0; }
  virtual double norm(const Eigen::Vector3d& r) { return std::pow(value(r),2); }

  virtual Eigen::Vector3d local_grad(const Eigen::Vector3d& r) { return {0, 0, 0}; }
  virtual Eigen::Vector3d grad(const Eigen::Vector3d& r) {
    return local_grad(r) * value(r);
  }

  virtual double local_laplacian(const Eigen::Vector3d& r) { return 0; }
  virtual double laplacian(const Eigen::Vector3d& r) {
    return local_laplacian(r) * value(r);
  }
};

class SphericalOscillatorWF : public WaveFunction {
 public:
  explicit SphericalOscillatorWF(double b);

  double value(const Eigen::Vector3d& r);
  double norm(const Eigen::Vector3d& r);

  Eigen::Vector3d local_grad(const Eigen::Vector3d& r);
  double local_laplacian(const Eigen::Vector3d& r);

 private:
  double b2_;
  double b4_;
};

// class AlphaParticleWF : public WaveFunction {
//  public:
//   AlphaParticleWF(double b_protons, double b_neutrons);
//
//   double value(const std::vector<Eigen::Vector3d>& positions);
//   double norm(const std::vector<Eigen::Vector3d>& positions);
//
//   double value(const std::vector<Eigen::Vector3d>& positions);
//
// };

}  // end namespace wf

#endif  // WAVEFUNCTION_H_
