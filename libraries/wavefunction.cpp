/****************************************************************
  wavefunction.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include "wavefunction.h"

namespace wf {

SphericalOscillatorWF::SphericalOscillatorWF(double b)
  : b2_(std::pow(b, 2)), b4_(std::pow(b,4)) {}

double SphericalOscillatorWF::value(const Eigen::Vector3d& r) {
  double r2 = r.dot(r);
  return std::exp(r2 / (-2 * b2_));
}

double SphericalOscillatorWF::norm(const Eigen::Vector3d& r) {
  double r2 = r.dot(r);
  return std::exp(r2 / (-1 * b2_));
}

Eigen::Vector3d SphericalOscillatorWF::local_grad(const Eigen::Vector3d& r) {
  return r * (-b2_);
}

double SphericalOscillatorWF::local_laplacian(const Eigen::Vector3d& r) {
  double r2 = r.dot(r);

  return r2/b4_ - 3/b2_;
}

}  // end namespace wf
