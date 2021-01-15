[![macOS](https://github.com/amirmasoudabdol/baaraan/workflows/macOS/badge.svg)](https://github.com/amirmasoudabdol/baaraan/actions?query=workflow%3AmacOS)

# Baaraan

Baaraan is a collection of some of the missing random number distributions in C++ Standard Library, STL. For the start, I am planning to provide some of the multivariate distributions, e.g., Multi-variate Normal Distribution, and slowly hopefully expanding the list and even adding some noise functions to the library.

## Design

For the most part, Baaraan's random number distributions are sharing the same interface and implementation details as their STL counterparts. First advantages of this is that they will look and behave very familiar. Moreover, they can adapt to your setup, for example, they will work with different URNGs out of the box.

I use [Armadillo](http://arma.sourceforge.net) as the backend Linear Algebra library. Mainly, I'm using `arma::mat` and `arma::vec` as default `MatrixType`  and `VectorType` types. This is shared between all multi-variate distributions.

For the actual probability distribution implementation, i.e., density, probability, cumulative functions, I choose to follow [Boost](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref.html) interface. In Boost implementation, several [non-member functions](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref/nmp.html) are being used to query different properties of a distribution.

## Distributions

- [Multivariate Normal Distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
- [Multivariate t-Student Distibution](https://en.wikipedia.org/wiki/Multivariate_t-distribution?wprov=sfti1)

## Dependencies

- STL
- Boost
- Armadillo

## Build and Install

```bash
mkdir build && cd build
cmake ..
make -j4
make install
```