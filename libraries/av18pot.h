/******************************************************************************/
/**
  @file av18pot.h

  Provides a nicer interface to the Intel MKL random number generators.

  @author Patrick J. Fasano
  University of Notre Dame

  + 12/14/16 (pjf): Created.

*******************************************************************************/
#ifndef AV18POT_H_
#define AV18POT_H_

extern "C" {
  void av18pw_(int* lpot, int* l, int* s, int* j,
               int* t, int* t1z, int* t2z, double* r, double v[2][2]);
  void av18op_(int* lpot, double* r, double vnn[]);
  void empot_(int* lpot, double* r, double vem[]);
  void consts_(int* lpot, double* mpi0, double* mpic, double* mp,
               double* mn, double* alpha, double* mup, double* mun);
}

#endif  // AV18POT_H-
