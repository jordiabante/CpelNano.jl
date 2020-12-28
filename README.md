# CpelNano

[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelNano.jl/blob/master/LICENSE)

## Description

CpelNano is the first method designed to perform DNA methylation differential analysis
using Oxford nanopore data. The package is based on a recently published method [1].
CpelNano not only detects mean methylation level differences between groups, but it can
also detect significant differences in methylation entropy as well as in the probability
distribution of methylation (see [1] for technical details).

## Testing

CpelNano is tested against Julia `1.3.0` on the latest versions of Linux, macOS and Windows.

## Getting Started

### Prerequisites

* julia v1.3.0
* git.

### Installing

`CpelNano` and dependencies can be installed via the following command in julia's REPL:

```julia
(v1.3) pkg> add https://github.com/jordiabante/CpelNano.jl.git
```

## Running the tests

In order to test the package to ensure that it has been properly installed,
run the following command in a `julia` session:

```julia
(v1.3) pkg> test CpelNano
```

If the package has been properly installed, then all tests will be successful.

## Authors

* **Jordi Abante**
* **Sandeep Kambhampati**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md)
file for details.

## References

[1] Abante, J., Kambhampati, S., Feinberg, A.P., Goutsias, J., CpelNano: differential
DNA methylation analysis with Oxford nanopore data, *Journal* 2020 XYZ.
