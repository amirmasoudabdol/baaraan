[![macOS](https://github.com/amirmasoudabdol/baaraan/workflows/macOS/badge.svg)](https://github.com/amirmasoudabdol/baaraan/actions?query=workflow%3AmacOS)
[![Ubuntu](https://github.com/amirmasoudabdol/baaraan/workflows/Ubuntu/badge.svg)](https://github.com/amirmasoudabdol/baaraan/actions?query=workflow%3AUbuntu)

<img src="docs/img/logo.png" width="250" align="right"/>

# Baaraan

Baaraan (Farsi: باران, baran | bârân |) is a collection of [STL's](https://en.cppreference.com/w/cpp/numeric/random) missing random number distributions for C++. For the start, some of the more important and useful distributions are implemented and tests, e.g., Multivariate Normal Distribution, and I am hopping that I can slowly expand the list of the distributions, and maybe even add some noise functions to the library as well.

> Baaraan is derived from [_baran_](https://en.wiktionary.org/wiki/باران), the Farsi word for _rain_. 

> 🚨 **Disclaimer:** I am in the process of testing and open-sourcing baaraan, so, please do not use it in production unless you have done the test yourself! 

For the most part, Baaraan's random number distributions are sharing the same interface and implementation details as their STL counterparts. First advantages of this is that they will look and behave very familiar. Moreover, they can adapt to your setup, for example, they will work with different URNGs out of the box.

I use [Armadillo](http://arma.sourceforge.net) as the backend Linear Algebra library. Mainly, I'm using `arma::vec` and `arma::mat` as default `VectorType`  and `MatrixType` types. This is shared between all multi-variate distributions.

For the actual probability distribution implementation, i.e., density, probability, cumulative functions, I choose to follow [Boost](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref.html) interface. Boost uses several [non-member functions](https://www.boost.org/doc/libs/1_62_0/libs/math/doc/html/math_toolkit/dist_ref/nmp.html) are being used to query different properties of a distribution.

## Available Distributions

**Multivariate:**
- [Multivariate Normal](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
- [Multivariate t-Student Distibution](https://en.wikipedia.org/wiki/Multivariate_t-distribution?wprov=sfti1)


**Truncated:**
- [Truncated Normal](https://en.wikipedia.org/wiki/Truncated_normal_distribution)
- [Truncated Multivariate Normal](https://en.wikipedia.org/wiki/Truncated_normal_distribution)


**Rectified:**
- [Rectified Normal](https://en.wikipedia.org/wiki/Rectified_Gaussian_distribution)

## Qucick Start

You can add baaraan to your project by copying it to your project folder, and setting it as a subfolder in your `CMakeLists.txt`. First you need to clone the repo, and copy it to your project folder:

```bash
git clone https://github.com/amirmasoudabdol/baaraan
cp -r baaraan your-project/
```

and now you can add it to your CMake project:

```cmake
add_subdirectory(baaraan)

# don't forget to link it to your executable or library
target_link_libraries(your-project baaraan) 
```

## Installation

Alternatively, you can build and install baaraan first, and use CMake's `find_package` command to add it to your project:

```bash
git clone https://github.com/amirmasoudabdol/baaraan
cd baaraan; mkdir build; cd build
cmake .. && make
make install
```

Now, you can simply use the `find_package` and if everything is setup correctly, you should be able to compile your projects and use baaraan in your project.

```cmake
find_package(baaraan)

# don't forget to link it to your executable or library
target_link_libraries(your-project baaraan) 
```

## Example

After installing and linking baaraan to your project, you should be able to simply `#include` any distributions and use it as follow:

```cpp
#include <iostream>
#include "baaraan/dists/mvnorm_distribution.h"

int main(int argc, char const *argv[])
{
	arma::mat sample(3, 10000);

	arma::vec means {1, 2, 3};
	arma::mat sigma{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	
	std::mt19937 gen(42);
	mvnorm_distribution<double> mvnorm{means, sigma};

	sample.each_col([&](arma::vec &v){v = mvnorm(gen);});

	std::cout << arma::mean(sample, 1);
	return 0;
}
```