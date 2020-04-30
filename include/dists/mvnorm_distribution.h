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
  typedef arma::mat result_type;
  typedef arma::Mat<RealType> matrix_type;
  typedef arma::Col<RealType> vector_type;

  class param_type {
    size_t dims_;
    result_type means_;
    result_type sigma_;

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
    }

    // TODO: I think I need a copy assignment operator for handling the
    // sizes and special cases

    size_t dims() const { return dims_; }

    vector_type means() const { return means_; }

    matrix_type sigma() const { return sigma_; }

    vector_type covs_diag() const { return sigma_.diag(); }

    // bool is_covs_diagmat() const { return covs_.is_diagmat(); }

    friend bool operator==(const param_type &x, const param_type &y) {
      return arma::approx_equal(x.means_, y.means_, "absdiff", 0.001) &&
             arma::approx_equal(x.covs_, y.covs_, "absdiff", 0.001);
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  matrix_type covs_lower;
  matrix_type inv_covs_lower;
  matrix_type inv_covs;
  std::normal_distribution<> norm; // N~(0, 1)

  param_type p_;
  result_type tmp_;

public:
  explicit mvnorm_distribution() : p_(param_type{}) { tmp_.resize(p_.dims()); };

  // constructor and reset functions
  explicit mvnorm_distribution(result_type means, result_type covs)
      : p_(param_type(means, covs)) {

    // TODO: check if it's diagonal, initiate the diag model

    tmp_.resize(p_.dims(), 1);

    // if (!p_.is_covs_diagmat())
    factorize_covariance();
  }

  explicit mvnorm_distribution(const param_type &p) : p_(p) {}

  void reset() { norm.reset(); };

  // generating functions
  template <class URNG> result_type operator()(URNG &g) {
    return (*this)(g, p_);
  }

  template <class URNG> result_type operator()(URNG &g, const param_type &parm);

  // property functions

  result_type means() const { return p_.means(); }

  result_type sigma() const { return p_.sigma(); }

  param_type param() const { return p_; }

  void param(const param_type &params) {
    // TODO: This needs more checks.
    p_ = params;

    tmp_.resize(p_.dims());

    // if (!p_.is_covs_diagmat())
    factorize_covariance();
  }

private:
  void factorize_covariance() {
    covs_lower = arma::chol(p_.sigma(), "lower");
    inv_covs_lower = arma::inv(arma::trimatl(covs_lower));
    inv_covs = inv_covs_lower.t() * inv_covs_lower;
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
mvnorm_distribution<double>::result_type
mvnorm_distribution<RealType>::operator()(
    URNG &g, const mvnorm_distribution<RealType>::param_type &parm) {

  tmp_.imbue([&]() { return norm(g); });
  // if (parm.is_covs_diagmat()) {
  //     return arma::sqrt(parm.covs_diag()) % tmp_ + parm.means();
  // } else {
  return covs_lower * tmp_ + parm.means();
  // }
}

} // namespace baaraan

#endif // BAARAAN_MVNORM_DISTRIBUTION_H
