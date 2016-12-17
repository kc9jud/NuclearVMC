/****************************************************************
  potential.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <cmath>
#include "eigen3/Eigen/Core"

#include "units.h"
#include "potential.h"

namespace potential {

Potential::Potential() {}

HarmonicOscillator::HarmonicOscillator(double b) : b4_(std::pow(b, 4)) {}

double HarmonicOscillator::operator()(const Eigen::Vector3d& r) const {
  return (units::kH_bar2 * r.dot(r)) / (2 * units::kM_p * b4_);
}

AnharmonicOscillator::AnharmonicOscillator(double b) : b6_(std::pow(b, 6)) {}

double AnharmonicOscillator::operator()(const Eigen::Vector3d& r) const {
  return (units::kH_bar2 * std::pow(r.dot(r), 2)) / (2 * units::kM_p * b6_);
}


AV4::AV4(int T, int S) {}

double AV4::operator()(double r) {};
}  // end namespace potential
