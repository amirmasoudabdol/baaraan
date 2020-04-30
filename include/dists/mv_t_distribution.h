//
// Created by Amir Masoud Abdol on 2019-10-28
//

#ifndef BAARAAN_MV_T_DISTRIBUTION_H
#define BAARAAN_MV_T_DISTRIBUTION_H

#include <random>
#include <iostream>
#include <armadillo>

namespace baaraan {

	template<class RealType = double>
    class mv_t_distribution {
    public:
        // types
        typedef arma::mat result_type;

        class param_type {
            size_t dims_;
            int dof_;
            result_type means_;
            result_type sigma_;

        public:
            typedef mv_t_distribution distribution_type;

            explicit param_type() {
            	dof_ = 2;
                arma::mat means(1, 1, arma::fill::zeros);
                arma::mat covs(1, 1, arma::fill::eye);

                dims_ = means.n_elem;
                means_ = means;
                sigma_ = covs;
            }

            explicit param_type(double dof, result_type means, result_type sigma)
                    : dof_(dof), dims_(means.n_elem), means_(means), sigma_(sigma) {

                if (dof <= 0)
                	throw std::logic_error("degress of freedom should be positive.");

                if (!means.is_colvec())
                    throw std::logic_error("Mean should be a column vector.");

                if (sigma.n_rows != dims_)
                    throw std::length_error("Covariance matrix has the wrong dimension.");

                if (!sigma.is_symmetric() || !sigma.is_square())
                    throw std::logic_error("Covarinace matrix is not square or symmetrical.");
            }

            // TODO: I think I need a copy assignment operator for handling the sizes and special cases

            size_t dims() const { return dims_; }

            double dof() const {return dof_; }

            result_type means() const { return means_; }

            result_type sigma() const { return sigma_; }

            arma::vec covs_diag() const { return sigma_.diag(); }

            // bool is_covs_diagmat() const { return covs_.is_diagmat(); }

            friend
            bool operator==(const param_type &x, const param_type &y) {
                return x.dof_ == y.dof_ 
                	   && arma::approx_equal(x.means_, y.means_, "absdiff", 0.001)
                       && arma::approx_equal(x.sigma_, y.sigma_, "absdiff", 0.001);
            }

            friend
            bool operator!=(const param_type &x, const param_type &y) { return !(x == y); }
        };

    private:

        arma::mat covs_lower;
        arma::mat inv_covs_lower;
        arma::mat inv_covs;
        std::normal_distribution<> norm;  // N~(0, 1)
        std::chi_squared_distribution<> chisq;

        param_type p_;
        result_type tmp_;

    public:

        explicit mv_t_distribution() : p_(param_type{}) {
            tmp_.resize(p_.dims());
        };

        // constructor and reset functions
        explicit mv_t_distribution(double dof, result_type means, result_type covs)
                : p_(param_type(dof, means, covs)) {

            // TODO: check if it's diagonal, initiate the diag model

            tmp_.resize(p_.dims(), 1);

            // if (!p_.is_covs_diagmat())
            factorize_covariance();
        }

        explicit mv_t_distribution(const param_type &p)
                : p_(p) {}

        void reset() { norm.reset(); };

        // generating functions
        template<class URNG>
        result_type operator()(URNG &g) { return (*this)(g, p_); }

        template<class URNG>
        result_type operator()(URNG &g, const param_type &parm);

        // property functions
        double dof() const {return p_.dof(); }

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

        void factorize_covariance() {
            covs_lower = arma::chol(p_.covs(), "lower");
            inv_covs_lower = arma::inv(arma::trimatl(covs_lower));
            inv_covs = inv_covs_lower.t() * inv_covs_lower;
        }

        result_type min() const { return -std::numeric_limits<RealType>::infinity(); }

        result_type max() const { return std::numeric_limits<RealType>::infinity(); }

        friend bool operator==(const mv_t_distribution &x,
                               const mv_t_distribution &y) { return x.p_ == y.p_; }

        friend bool operator!=(const mv_t_distribution &x,
                               const mv_t_distribution &y) { return !(x == y); }

        template<class charT, class traits>
        friend
        std::basic_ostream<charT, traits> &
        operator<<(std::basic_ostream<charT, traits> &os,
                   const mv_t_distribution &means);

        template<class charT, class traits>
        friend
        std::basic_istream<charT, traits> &
        operator>>(std::basic_istream<charT, traits> &is,
                   mv_t_distribution &means);

    };

    template<class RealType>
    template<class URNG>
    mv_t_distribution<double>::result_type
    mv_t_distribution<RealType>::operator()(URNG &g, const mv_t_distribution<RealType>::param_type &parm) {

        tmp_.imbue([&]() { return norm(g); });
        // if (parm.is_covs_diagmat()) {
        return arma::sqrt(parm.dof() / parm.covs_diag()) * tmp_ + parm.means();
        // } else {
        //     return covs_lower * tmp_ + parm.means();
        // }

    }

}


#endif // BAARAAN_MV_T_DISTRIBUTION_H
