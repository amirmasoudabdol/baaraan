[![macOS](https://github.com/amirmasoudabdol/baaraan/workflows/macOS/badge.svg)](https://github.com/amirmasoudabdol/baaraan/actions?query=workflow%3AmacOS)
[![Ubuntu](https://github.com/amirmasoudabdol/baaraan/workflows/Ubuntu/badge.svg)](https://github.com/amirmasoudabdol/baaraan/actions?query=workflow%3AUbuntu)

<img src="docs/img/logo.png" width="250" align="right"/>

# Baaraan

Baaraan (Farsi: باران, baran | bârân |) is a collection of missing random number distributions for C++. For the start, some of the more important and useful distributions are implemented and tests, e.g., Multivariate Normal Distribution, and I am hopping that I can slowly expand the list of the distributions, and maybe even add some noise functions to the library as well.

> Baaraan is derived from [_baran_](https://en.wiktionary.org/wiki/باران), the Farsi word for _rain_. 

## Design

For the most part, Baaraan's random number distributions are sharing the same interface and implementation details as their STL counterparts. First advantages of this is that they will look and behave very familiar. Moreover, they can adapt to your setup, for example, they will work with different URNGs out of the box.

I use [Armadillo](http://arma.sourceforge.net) as the backend Linear Algebra library. Mainly, I'm using `arma::vec` and `arma::mat` as default `VectorType`  and `MatrixType` types. This is shared between all multi-variate distributions.

For the actual probability distribution implementation, i.e., density, probability, cumulative functions, I choose to follow [Boost](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref.html) interface. Boost uses several [non-member functions](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref/nmp.html) are being used to query different properties of a distribution.

## Available Distributions

### Multivariate 

- [Multivariate Normal](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
- [Multivariate t-Student Distibution](https://en.wikipedia.org/wiki/Multivariate_t-distribution?wprov=sfti1)


### Truncated

- [Truncated Normal](https://en.wikipedia.org/wiki/Truncated_normal_distribution)
- [Truncated Multivariate Normal](https://en.wikipedia.org/wiki/Truncated_normal_distribution)


### Rectified

- [Rectified Normal](https://en.wikipedia.org/wiki/Rectified_Gaussian_distribution)

## Installation

### Dependencies

- STL
- Boost
- Armadillo

### Build and Install

```bash
mkdir build && cd build
cmake ..
make -j4
make install
```


## Example