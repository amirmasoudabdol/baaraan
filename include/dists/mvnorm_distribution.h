//
// Created by Amir Masoud Abdol on 2019-06-14.
//
/// \file
/// This file contains the implementation of the multi-variate normal
/// random distribution

#ifndef BAARAAN_MVNORM_DISTRIBUTION_H
#define BAARAAN_MVNORM_DISTRIBUTION_H

#include <armadillo>
#include <iostream>
#include <random>

namespace baaraan {

template <class RealType = double> class mvnorm_distribution {
public:
  // types
  typedef arma::Mat<RealType> matrix_type;
  typedef arma::Col<RealType> vector_type;

  class param_type {
    size_t dims_;
    vector_type means_;
    matrix_type sigma_;

  private:
    matrix_type covs_lower_;
    matrix_type inv_covs_lower_;
    matrix_type inv_covs_;

    void factorize_covariance() {
      covs_lower_ = arma::chol(sigma_, "lower");
      inv_covs_lower_ = arma::inv(arma::trimatl(covs_lower_));
      inv_covs_ = inv_covs_lower_.t() * inv_covs_lower_;
    }

  public:
    typedef mvnorm_distribution distribution_type;

    explicit param_type(vector_type means, matrix_type sigma)
        : dims_(means.n_elem), means_(means), sigma_(sigma) {

      if (!means.is_colvec())
        throw std::logic_error("Mean should be a column vector.");

      if (sigma.n_rows != dims_)
        throw std::length_error("Covariance matrix has the wrong dimension.");

      if (!sigma.is_symmetric() || !sigma.is_square())
        throw std::logic_error(
            "Covarinace matrix is not square or symmetrical.");

      factorize_covariance();
    }

    /// TODO: I need a copy assignment operator for sure

    size_t dims() const { return dims_; }

    vector_type means() const { return means_; }

    matrix_type sigma() const { return sigma_; }

    matrix_type covs_lower() const { return covs_lower_; }
    matrix_type inv_covs() const { return inv_covs_; }

    friend bool operator==(const param_type &x, const param_type &y) {
      return arma::approx_equal(x.means_, y.means_, "absdiff", 0.001) &&
             arma::approx_equal(x.covs_, y.covs_, "absdiff", 0.001);
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  std::normal_distribution<> norm_; // N~(0, 1)

  param_type p_;
  vector_type v_;

public:
  // constructor and reset functions
  explicit mvnorm_distribution(vector_type means, matrix_type sigma)
      : p_(param_type(means, sigma)) {}

  explicit mvnorm_distribution(const param_type &p) : p_(p) {}

  void reset() { norm_.reset(); };

  // generating functions
  template <class URNG> vector_type operator()(URNG &g) {
    return (*this)(g, p_);
  }

  template <class URNG> vector_type operator()(URNG &g, const param_type &p);

  // batch generation
  template <class URNG> matrix_type operator()(URNG &g, size_t n) {
    return (*this)(g, p_, n);
  }

  template <class URNG> matrix_type operator()(URNG &g, const param_type &p, size_t n);

  // property functions

  vector_type means() const { return p_.means(); }

  matrix_type sigma() const { return p_.sigma(); }

  param_type param() const { return p_; }

  void param(const param_type &p) {
    // TODO: This needs more checks.
    p_ = p;

    //    factorize_covariance();
  }

public:
  vector_type min() const {
    return vector_type(p_.dims()).fill(
        -std::numeric_limits<RealType>::infinity());
  }

  vector_type max() const {
    return vector_type(p_.dims()).fill(
        +std::numeric_limits<RealType>::infinity());
  }

  friend bool operator==(const mvnorm_distribution &x,
                         const mvnorm_distribution &y) {
    return x.p_ == y.p_;
  }

  friend bool operator!=(const mvnorm_distribution &x,
                         const mvnorm_distribution &y) {
    return !(x == y);
  }

  template <class charT, class traits>
  friend std::basic_ostream<charT, traits> &
  operator<<(std::basic_ostream<charT, traits> &os,
             const mvnorm_distribution &means);

  template <class charT, class traits>
  friend std::basic_istream<charT, traits> &
  operator>>(std::basic_istream<charT, traits> &is, mvnorm_distribution &means);
};

template <class RealType>
template <class URNG>
typename mvnorm_distribution<RealType>::vector_type
mvnorm_distribution<RealType>::operator()(
    URNG &g, const mvnorm_distribution<RealType>::param_type &p) {

  v_.resize(p.dims(), 1);
  v_.imbue([&]() { return norm_(g); });
  // if (p.is_covs_diagmat()) {
  //     return arma::sqrt(p.covs_diag()) % v_ + p.means();
  // } else {
  return p.covs_lower() * v_ + p.means();
  // }
}

template <class RealType>
template <class URNG>
typename mvnorm_distribution<RealType>::matrix_type
mvnorm_distribution<RealType>::operator()(
    URNG &g, const mvnorm_distribution<RealType>::param_type &p, size_t n) {

  arma::mat res(p.dims(), n);

  res.each_col([&](vector_type &col){col = (*this)(g, p); });

  return res;

}


} // namespace baaraan

#endif // BAARAAN_MVNORM_DISTRIBUTION_H
