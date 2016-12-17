/************************************************************/
/**
  @file wavefunction.h

  VMC trial wavefunctions.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 12/14/16 (pjf): Created.

****************************************************************/

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include "av18pot.h"

namespace potential {

class Potential {
 public:
  Potential();

  virtual double operator()(double r) { return 0; }
};

class HarmonicOscillator {
 public:
  explicit HarmonicOscillator(double b);

  double operator()(const Eigen::Vector3d& r) const;

 private:
  double b4_;
};

class AnharmonicOscillator {
 public:
  explicit AnharmonicOscillator(double b);

  double operator()(const Eigen::Vector3d& r) const;

 private:
  double b6_;
};

class AV4 : public Potential {
public:
  AV4(int T, int S);

  double operator()(double r);
};
}  // end namespace potential

#endif  // POTENTIAL_H_
