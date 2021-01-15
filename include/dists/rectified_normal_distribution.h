//
// Created by Amir Masoud Abdol on 12/03/2020.
//

#ifndef BAARAAN_RECTIFIED_NORMAL_DISTRIBUTION_H
#define BAARAAN_RECTIFIED_NORMAL_DISTRIBUTION_H

#include <iostream>
#include <random>

using boost::math::normal;

namespace baaraan {

template <class RealType = double> class rectified_normal_distribution {
public:
  // types
  typedef RealType result_type;

  class param_type {
    result_type mean_;
    result_type stddev_;

  public:
    typedef rectified_normal_distribution distribution_type;

    explicit param_type(result_type mean = 0, result_type stddev = 1)
        : mean_(mean), stddev_(stddev) {}

    result_type mean() const { return mean_; }

    result_type stddev() const { return stddev_; }

    friend bool operator==(const param_type &x, const param_type &y) {
      return x.mean_ == y.mean_ && x.stddev_ == y.stddev_;
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  param_type p_;
  std::normal_distribution<> norm_;

public:
  // constructors and reset functions
  explicit rectified_normal_distribution(result_type mean = 0,
                                         result_type stddev = 1)
      : p_(param_type(mean, stddev)), norm_(std::normal_distribution<>::param_type{mean, stddev}) {}

  explicit rectified_normal_distribution(const param_type &p) : p_(p) {}

  void reset() { norm_.reset(); }

  // generating functions
  template <class URNG> result_type operator()(URNG &g) {
    return (*this)(g, p_);
  }

  template <class URNG> result_type operator()(URNG &g, const param_type &p);

  // property functions
  result_type mean() const { return p_.mean(); }

  result_type stddev() const { return p_.stddev(); }

  param_type param() const { return p_; }

  void param(const param_type &p) { p_ = p; }

  result_type min() const { return 0; }

  result_type max() const { return +std::numeric_limits<RealType>::infinity(); }

  friend bool operator==(const rectified_normal_distribution &x,
                         const rectified_normal_distribution &y) {
    return x.p_ == y.p_;
  }

  friend bool operator!=(const rectified_normal_distribution &x,
                         const rectified_normal_distribution &y) {
    return !(x == y);
  }

  template <class _CharT, class _Traits, class _RT>
  friend std::basic_ostream<_CharT, _Traits> &
  operator<<(std::basic_ostream<_CharT, _Traits> &os,
             const rectified_normal_distribution<_RT> &x);

  template <class _CharT, class _Traits, class _RT>
  friend std::basic_istream<_CharT, _Traits> &
  operator>>(std::basic_istream<_CharT, _Traits> &is,
             rectified_normal_distribution<_RT> &x);
};

template <class RealType>
template <class URNG>
RealType
rectified_normal_distribution<RealType>::operator()(URNG &g,
                                                    const param_type &parm) {
  return max(0, norm_(g, parm));
}

}; // namespace baaraan

#endif // BAARAAN_RECTIFIED_NORMAL_DISTRIBUTION_H
