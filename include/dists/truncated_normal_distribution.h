//
// Created by Amir Masoud Abdol on 12/03/2020.
//

#ifndef BAARAAN_TRUNCATED_NORMAL_DISTRIBUTION_H
#define BAARAAN_TRUNCATED_NORMAL_DISTRIBUTION_H

#include "boost/math/distributions/normal.hpp"
#include <iostream>
#include <random>

using boost::math::normal;

namespace baaraan {

template <class RealType = double> class truncated_normal_distribution {
public:
  // types
  typedef RealType result_type;

  class param_type {
    result_type mean_;
    result_type stddev_;
    result_type lower_;
    result_type upper_;

  public:
    typedef truncated_normal_distribution distribution_type;

    explicit param_type(result_type mean = 0, result_type stddev = 1,
                        result_type lower = -3, result_type upper = 3)
        : mean_(mean), stddev_(stddev), lower_(lower), upper_(upper) {}

    result_type mean() const { return mean_; }

    result_type stddev() const { return stddev_; }

    result_type lower() const { return lower_; }

    result_type upper() const { return upper_; }

    friend bool operator==(const param_type &x, const param_type &y) {
      return x.mean_ == y.mean_ && x.stddev_ == y.stddev_ &&
             x.lower_ == y.lower_ && x.upper_ == y.upper_;
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  param_type p_;
  normal unit_normal_;
  std::uniform_real_distribution<> uniform_;

public:
  // constructors and reset functions
  explicit truncated_normal_distribution(result_type mean = 0,
                                         result_type stddev = 1,
                                         result_type lower = -4,
                                         result_type upper = 4)
      : p_(param_type(mean, stddev, lower, upper)) {}

  explicit truncated_normal_distribution(const param_type &p) : p_(p) {}

  void reset() { uniform_.reset(); }

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

  result_type min() const { return p_.lower(); }

  result_type max() const { return p_.upper(); }

  friend bool operator==(const truncated_normal_distribution &x,
                         const truncated_normal_distribution &y) {
    return x.p_ == y.p_;
  }

  friend bool operator!=(const truncated_normal_distribution &x,
                         const truncated_normal_distribution &y) {
    return !(x == y);
  }

  template <class _CharT, class _Traits, class _RT>
  friend std::basic_ostream<_CharT, _Traits> &
  operator<<(std::basic_ostream<_CharT, _Traits> &os,
             const truncated_normal_distribution<_RT> &x);

  template <class _CharT, class _Traits, class _RT>
  friend std::basic_istream<_CharT, _Traits> &
  operator>>(std::basic_istream<_CharT, _Traits> &is,
             truncated_normal_distribution<_RT> &x);
};

template <class RealType>
template <class URNG>
RealType
truncated_normal_distribution<RealType>::operator()(URNG &g,
                                                    const param_type &parm) {
  double alpha = (parm.lower() - parm.mean()) / parm.stddev();
  double beta = (parm.upper() - parm.mean()) / parm.stddev();

  double alpha_cdf = cdf(unit_normal_, alpha);
  double beta_cdf = cdf(unit_normal_, beta);

  double u = uniform_(g);
  double xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf);
  double xi = quantile(unit_normal_, xi_cdf);

  double x = parm.mean() + parm.stddev() * xi;

  return x;
}

}; // namespace baaraan

#endif // BAARAAN_TRUNCATED_NORMAL_DISTRIBUTION_H
