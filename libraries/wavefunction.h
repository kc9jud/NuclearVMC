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

class SingleParticleWF {
 public:
  SingleParticleWF() {}

  virtual double value(const Eigen::Vector3d& r) { return 0; }
  virtual double norm(const Eigen::Vector3d& r) { return 0; }

  virtual Eigen::Vector3d local_grad(const Eigen::Vector3d& r) { return {0, 0, 0}; }
  virtual Eigen::Vector3d grad(const Eigen::Vector3d& r) {
    return local_grad(r) * value(r);
  }

  virtual double local_laplacian(const Eigen::Vector3d& r) { return 0; }
  virtual double laplacian(const Eigen::Vector3d& r) {
    return local_laplacian(r) * value(r);
  }
};

class SphericalOscillatorWF : public SingleParticleWF {
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

}  // end namespace wf

#endif  // WAVEFUNCTION_H_
