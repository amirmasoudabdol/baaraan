//
// Created by Amir Masoud Abdol on 2019-10-28
//

#ifndef BAARAAN_MV_T_DISTRIBUTION_H
#define BAARAAN_MV_T_DISTRIBUTION_H

#include <armadillo>
#include <iostream>
#include <random>

namespace baaraan {

template <class RealType = double> class mv_t_distribution {
public:
  // types
  typedef arma::Mat<RealType> matrix_type;
  typedef arma::Col<RealType> vector_type;

  class param_type {
    size_t dims_;
    int dof_;
    vector_type means_;
    matrix_type sigma_;

    arma::mat covs_lower_;
    arma::mat inv_covs_lower_;
    arma::mat inv_covs_;

    void factorize_covariance() {
      covs_lower_ = arma::chol(sigma_, "lower");
      inv_covs_lower_ = arma::inv(arma::trimatl(covs_lower_));
      inv_covs_ = inv_covs_lower_.t() * inv_covs_lower_;
    }

    matrix_type covs_lower() const { return covs_lower_; }
    matrix_type inv_covs() const { return inv_covs_; }

  public:
    typedef mv_t_distribution distribution_type;

    explicit param_type(double dof, vector_type means, matrix_type sigma)
        : dof_(dof), dims_(means.n_elem), means_(means), sigma_(sigma) {

      if (dof <= 0)
        throw std::logic_error("degress of freedom should be positive.");

      if (!means.is_colvec())
        throw std::logic_error("Mean should be a column vector.");

      if (sigma.n_rows != dims_)
        throw std::length_error("Covariance matrix has the wrong dimension.");

      if (!sigma.is_symmetric() || !sigma.is_square())
        throw std::logic_error(
            "Covarinace matrix is not square or symmetrical.");

      factorize_covariance();
    }

    size_t dims() const { return dims_; }

    double dof() const { return dof_; }

    vector_type means() const { return means_; }

    matrix_type sigma() const { return sigma_; }


    friend bool operator==(const param_type &x, const param_type &y) {
      return x.dof_ == y.dof_ &&
             arma::approx_equal(x.means_, y.means_, "absdiff", 0.001) &&
             arma::approx_equal(x.sigma_, y.sigma_, "absdiff", 0.001);
    }

    friend bool operator!=(const param_type &x, const param_type &y) {
      return !(x == y);
    }
  };

private:
  std::normal_distribution<> norm; // N~(0, 1)
  std::chi_squared_distribution<> chisq;

  param_type p_;
  vector_type v_;

public:
  // constructor and reset functions

  explicit mv_t_distribution(const param_type &p) : p_(p) {}

  explicit mv_t_distribution(double dof, vector_type means, matrix_type sigma)
      : p_(param_type(dof, means, sigma)) {}

  void reset() { norm.reset(); };

  // generating functions
  template <class URNG> vector_type operator()(URNG &g) {
    return (*this)(g, p_);
  }

  template <class URNG> vector_type operator()(URNG &g, const param_type &p);

  // property functions
  double dof() const { return p_.dof(); }

  vector_type means() const { return p_.means(); }

  vector_type sigma() const { return p_.sigma(); }

  param_type param() const { return p_; }

  void param(const param_type &params) {
    p_ = params;
  }

  vector_type min() const {
    return vector_type(p_.dims()).fill(
        -std::numeric_limits<RealType>::infinity());
  }

  vector_type max() const {
    return vector_type(p_.dims()).fill(
        +std::numeric_limits<RealType>::infinity());
  }

  friend bool operator==(const mv_t_distribution &x,
                         const mv_t_distribution &y) {
    return x.p_ == y.p_;
  }

  friend bool operator!=(const mv_t_distribution &x,
                         const mv_t_distribution &y) {
    return !(x == y);
  }

  template <class charT, class traits>
  friend std::basic_ostream<charT, traits> &
  operator<<(std::basic_ostream<charT, traits> &os,
             const mv_t_distribution &means);

  template <class charT, class traits>
  friend std::basic_istream<charT, traits> &
  operator>>(std::basic_istream<charT, traits> &is, mv_t_distribution &means);
};

template <class RealType>
template <class URNG>
typename mv_t_distribution<RealType>::vector_type
mv_t_distribution<RealType>::operator()(
    URNG &g, const mv_t_distribution<RealType>::param_type &p) {

  v_.resize(p_.dims(), 1);

  v_.imbue([&]() { return norm(g); });
  // if (p.is_covs_diagmat()) {
  return arma::sqrt(p.dof() / p.covs_diag()) * v_ + p.means();
  // } else {
  //     return covs_lower * v_ + p.means();
  // }
}

} // namespace baaraan

#endif // BAARAAN_MV_T_DISTRIBUTION_H
