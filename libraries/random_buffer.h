/******************************************************************************/
/**
  @file random_buffer.h

  Provides a nicer interface to the Intel MKL random number generators.

  @author Patrick J. Fasano
  University of Notre Dame

  + 12/14/16 (pjf): Created.

*******************************************************************************/
#ifndef RANDOM_BUFFER_H_
#define RANDOM_BUFFER_H_

#include <mkl.h>
#include <vector>

namespace buffer {

template <typename T>
class RandomBuffer {
 public:
  RandomBuffer(size_t buffer_size, VSLStreamStatePtr& stream, T min, T max)
    : buffer_size_(buffer_size), stream_(stream), min_(min), max_(max)
  {
    buffer_.reserve(buffer_size);
    Fill();
  }

  // std::vector-like interfaces
  typedef typename std::vector<T>::size_type size_type;
        T& operator[](size_type idx)       { return buffer_[idx]; }
  const T& operator[](size_type idx) const { return buffer_[idx]; }
        T& at(size_type idx)       { return buffer_.at(idx); }
  const T& at(size_type idx) const { return buffer_.at(idx); }
  size_t size() const { return buffer_size; }

  void Fill();

 private:
  size_t buffer_size_;
  VSLStreamStatePtr stream_;
  std::vector<T> buffer_;
  T min_;
  T max_;
};

template <>
void RandomBuffer<double>::Fill() {
  vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, buffer_size_, buffer_.data(), min_, max_);
}

template <>
void RandomBuffer<float>::Fill() {
  vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, buffer_size_, buffer_.data(), min_, max_);
}

template <>
void RandomBuffer<int>::Fill() {
  viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream_, buffer_size_, buffer_.data(), min_, max_);
}

}  // namespace buffer
#endif  // RANDOM_BUFFER_H_
