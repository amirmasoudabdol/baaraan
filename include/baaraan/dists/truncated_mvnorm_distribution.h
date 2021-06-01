///
/// Created by Amir Masoud Abdol on 23/10/2019.
///
/// @file
/// This file contains the implementation of truncated normal random
/// distribution.
/// 
/// @defgroup   TruncatedDistributions Truncated Random Distributions
/// @brief      List of truncated random distributions
/// 

#ifndef BAARAAN_TRUNCATED_MVNORM_DISTRIBUTION_H
#define BAARAAN_TRUNCATED_MVNORM_DISTRIBUTION_H

#include "boost/math/distributions/normal.hpp"
#include <armadillo>
#include <iostream>
#include <random>

using boost::math::normal;

namespace baaraan {

///
/// @brief      Truncated Normal Distribution
///             
/// @tparam     RealType  Indicates the type of return values
/// 
/// @ingroup    TruncatedDistributions
/// @ingroup    MultivariateDistributions
///
template <class RealType = double> class truncated_mvnorm_distribution {
public:
  // types
  typedef arma::Col<RealType> vector_type;
  typedef arma::Mat<RealType> matrix_type;

  class param_type {
    size_t dims_;
    vector_type means_;
    matrix_type sigma_;
    vector_type lowers_;
    vector_type uppers_;

  public:
    typedef truncated_mvnorm_distribution distribution_type;

    explicit param_type(vector_type means, matrix_type sigma,
                        vector_type lowers, vector_type uppers)
        : means_(means), sigma_(sigma), lowers_(lowers), uppers_(uppers) {

      dims_ = means.n_elem;

      // @todo Check if lower is actually lower than the upper

      // Checking whether dimensions matches
      if (lowers_.n_elem != dims_ || uppers_.n_elem != dims_)
        throw std::length_error("Check your arrays size");

      if (!sigma.is_symmetric() || !sigma.is_square())
        throw std::logic_error("Covariance matrix is not symmetric.");
    }

    size_t dims() const { return dims_; }

    vector_type means() const { return means_; }

    matrix_type sigma() const { return sigma_; }

    vector_type lowers() const { return lowers_; }

    vector_type uppers() const { return uppers_; }

    friend bool operator==(const param_type &x, const param_type &y) {
      return arma::approx_equal(x.means_, y.means_, "absdiff", 0.001) &&
             arma::approx_equal(x.covs_, y.covs_, "absdiff", 0.001) &&
             arma::approx_equal(x.lowers_, y.lowers_, "absdiff", 0.001) &&
             arma::approx_equal(x.uppers_, y.uppers_, "absdiff", 0.001);
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  std::uniform_real_distribution<> uniform{};
  param_type p_;

  arma::mat sub1(arma::mat x, int i) {
    x.shed_col(i);
    x.shed_row(i);
    return x;
  }

  arma::mat sub2(arma::mat x, int a, int b) {
    x.shed_col(b);
    return (x.row(a));
  }

  arma::vec negSubCol(arma::vec x, int i) {
    x.shed_row(i);
    return (x);
  }

  arma::rowvec negSubRow(arma::rowvec x, int i) {
    x.shed_col(i);
    return (x);
  }

public:
  ///
  /// @brief      Constructs an instance of the truncated multivariate normal 
  /// random distribution by accepting its individual parameters
  ///
  /// @param[in]  means   The mean vector
  /// @param[in]  sigma   The covariance matrix
  /// @param[in]  lowers  The vector of values indicating lower truncation bound
  /// @param[in]  uppers  The vector of values indicating upper truncation bound
  ///
  explicit truncated_mvnorm_distribution(vector_type means, matrix_type sigma,
                                         vector_type lowers, vector_type uppers)
      : p_(param_type(means, sigma, lowers, uppers)) {}

  ///
  /// @brief      Constructs an instance of the truncated multivariate normal 
  /// random distribution by accepting an initialized 
  /// truncated_mvnorm_distribution::param_type.
  ///
  /// @param[in]  p     
  ///
  explicit truncated_mvnorm_distribution(const param_type &p) : p_(p) {}

  void reset() { uniform.reset(); };

  // generating functions
  template <class URNG> vector_type operator()(URNG &g) {
    return (*this)(g, p_);
  }

  template <class URNG> vector_type operator()(URNG &g, const param_type &p);

  // property functions
  vector_type means() const { return p_.means(); }

  matrix_type sigma() const { return p_.sigma(); }

  param_type param() const { return p_; };

  void param(const param_type &params) { p_ = params; }

  vector_type lowers() const { return p_.lowers(); }
  vector_type min() const { return p_.lowers(); }

  vector_type uppers() const { return p_.uppers(); }
  vector_type max() const { return p_.uppers(); }

  friend bool operator==(const truncated_mvnorm_distribution &x,
                         const truncated_mvnorm_distribution &y) {
    return x.p_ == y.p_;
  }

  friend bool operator!=(const truncated_mvnorm_distribution &x,
                         const truncated_mvnorm_distribution &y) {
    return !(x == y);
  }

  template <class charT, class traits>
  friend std::basic_ostream<charT, traits> &
  operator<<(std::basic_ostream<charT, traits> &os,
             const truncated_mvnorm_distribution &means);

  template <class charT, class traits>
  friend std::basic_istream<charT, traits> &
  operator>>(std::basic_istream<charT, traits> &is,
             truncated_mvnorm_distribution &means);
};

// Implementation of rejection algorithm
template <class RealType>
template <class _URNG>
typename truncated_mvnorm_distribution<RealType>::vector_type
truncated_mvnorm_distribution<RealType>::operator()(
    _URNG &g, const truncated_mvnorm_distribution<RealType>::param_type &p) {

  auto n{1};
  auto d = p.dims();
  arma::mat trace = arma::zeros(n, d); // trace of MCMC chain

  // draw from U(0,1)
  arma::vec U(n * d);
  U.imbue([&]() { return uniform(g); });

  auto l{0}; // iterator for U

  // calculate conditional standard deviations
  arma::vec sd(d);
  arma::cube P = arma::zeros(1, d - 1, d);

  for (int i = 0; i < d; i++) {
    // partitioning of sigma
    arma::mat Sigma = sub1(p.sigma(), i);
    double sigma_ii = p.sigma()(i, i);
    arma::rowvec Sigma_i = sub2(p.sigma(), i, i);

    P.slice(i) = Sigma_i * Sigma.i();
    double p_i = arma::as_scalar(P.slice(i) * Sigma_i.t());
    sd(i) = sqrt(sigma_ii - p_i);
  }

  arma::vec x = p.means();

  // run Gibbs sampler for specified chain length (MCMC chain of n samples)
  for (int j = 0; j < n; j++) {

    // sample all conditional distributions
    for (int i = 0; i < d; i++) {

      // calculation of conditional expectation and conditional variance
      arma::rowvec slice_i = P.slice(i);
      arma::vec slice_i_times = slice_i * (negSubCol(x, i) - negSubCol(x, i));
      double slice_i_times_double = arma::as_scalar(slice_i_times);
      double mu_i = p.means()(i) + slice_i_times_double;

      // transformation
      double Fa = cdf(normal{mu_i, sd(i)}, p.lowers()(i));
      double Fb = cdf(normal{mu_i, sd(i)}, p.uppers()(i));

      x(i) = mu_i + sd(i) * quantile(normal{0, 1}, U(l) * (Fb - Fa) + Fa);

      l = l + 1;
    }

    trace.row(j) = x.t();
  }

  return trace.t();
}

} // namespace baaraan

#endif // BAARAAN_TRUNCATED_MVNORM_DISTRIBUTION_H
