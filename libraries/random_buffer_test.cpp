#include <iostream>

#include "random_buffer.h"

int main() {
  // const size_t kNumRandom = 2;
  // const int kNumBlocks = 512000000;
  const size_t kNumRandom = 134217728;
  const int kNumBlocks = 8;
  VSLStreamStatePtr stream;
  int i, j;

  /* Initializing */
  double s = 0.0;
  vslNewStream(&stream, VSL_BRNG_SFMT19937, 777);
  buffer::RandomBuffer<double> r(kNumRandom, stream, 0.0, 5.0);

  /* Generating */
  #pragma omp parallel for reduction(+:s)
  for ( i=0; i < kNumBlocks; i++ ) {
    r.Fill();

    double ps = 0.0;
    #pragma omp simd reduction(+:ps)
    for ( j=0; j < kNumRandom; j++ ) {
       ps += r[j];
    }
    s += ps / (kNumRandom);
  }
  s /= kNumBlocks;

  /* Deleting the stream */
  vslDeleteStream(&stream);

  /* Printing results */
  std::cout << "Sample mean of uniform distribution = " << s << std::endl;

  return 0;
}
