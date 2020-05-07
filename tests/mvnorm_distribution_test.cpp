//
// Created by Amir Masoud Abdol on 24/10/2019.
//

#define BOOST_TEST_MODULE MVNORM_DISTRIBUTION TEST
#define BOOST_TEST_DYN_LINK

#include <random>
#include <iostream>

#include "boost/test/unit_test.hpp"
#include "boost/histogram/histogram.hpp"

#include "dists/mvnorm_distribution.h"

using namespace baaraan;

BOOST_AUTO_TEST_CASE( mvnorm_constructor )
{

}

BOOST_AUTO_TEST_CASE( mvnorm_means_test )
{

  arma::vec tmeans {1, 1, 1};
  arma::mat tsigma{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  mvnorm_distribution<double> mvnorm{tmeans, tsigma};

  std::mt19937 gen(42);

  arma::mat sample(3, 10000);

  sample.each_col([&](arma::vec &v){v = mvnorm(gen);});

  arma::vec means = arma::mean(sample, 1);

  BOOST_CHECK( approx_equal(tmeans, means, "absdiff", 0.01) );

}

BOOST_AUTO_TEST_CASE( mvnorm_stddevs_test )
{
  arma::vec tmeans {1, 1, 1};
  arma::mat tsigma{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  mvnorm_distribution<double> mvnorm{tmeans, tsigma};

  std::mt19937 gen(42);

  arma::mat sample(3, 10000);

  sample.each_col([&](arma::vec &v){v = mvnorm(gen);});

  arma::vec stddevs = arma::stddev(sample, 1, 1);

  BOOST_CHECK( approx_equal(stddevs, tsigma.diag(), "absdiff", 0.01) );
}