## 2BSM-ScannerS-with-Lilith3
2BSM-ScannerS-with-Lilith3 is a modified public code from the famous ScannerS-2 & Lilith2 code which scans parameters and checks points of parameters for standard-scalar theories (BSM). We use this code to study primarily Theory and phenomenology of two-Higgs-doublet models.

## ScannerS

A modern c++ parameter scanner for models with extended scalar sectors.

**Take a look at the [ScannerS manual], the [online documentation] or, if you want to extend ScannerS, at the [contribution guide].**

## Authors
ScannerS is written by [Jonas Wittbrodt], partly based on the old ScannerS code
by Raul Coimbra, Marco Sampaio, and Rui Santos.

Some of the model and constraint implementations in ScannerS have been
contributed by the following kind people:
  - Philipp Basler (CxSMBroken, and large parts of the BSMPT interface and EWPT
    constraint)
  - Tizian Römer (CxSMDark)

### References

  - R. Coimbra, M. O. P. Sampaio and R. Santos  
    "ScannerS: Constraining the phase diagram of a complex scalar singlet at the LHC"  
    Eur. Phys. J. C (2013) 73:2428, arXiv:1301.2599 [hep-ph]

  - P.M. Ferreira, R. Guedes, M. O. P. Sampaio, R. Santos  
    "Wrong sign and symmetric limits and non-decoupling in 2HDMs"  
    JHEP 1412 (2014) 067, arXiv:1409.6723 [hep-ph]

  - R. Costa, M. Mühlleitner, M. O. P. Sampaio, R. Santos  
    "Singlet Extensions of the Standard Model at LHC Run 2: Benchmarks and Comparison with the NMSSM"  
    JHEP 1606 (2016) 034, arXiv:1512.05355 [hep-ph]

  - M. Mühlleitner, M. O. P. Sampaio, R. Santos, J. Wittbrodt  
   "The N2HDM under Theoretical and Experimental Scrutiny"  
    JHEP 1703 (2017) 094 arXiv:1612.01309 [hep-ph]

  - M. Mühlleitner, M. O. P. Sampaio, R. Santos, J. Wittbrodt  
   "ScannerS: Parameter Scans in Extended Scalar Sectors"  
    arXiv:2007.02985 [hep-ph]

Please make sure to also always cite the relevant references for the model and
constraints that you use.

## Installation
The code is compiled using CMake.

On most systems you should be able to compile a working ScannerS setup through
```
mkdir build && cd build
cmake ..
make
```

### Dependencies

#### Required - Manual

These dependencies are required to be installed by the user. However, these are
widely used and should be available on most systems.
  - Working compilers for C++, C and Fortran. The C++ compiler must support
    `c++-17` (e.g. `gcc-7` or newer). On Mac, you may need to set the `CC` and
    `CXX` environment variables to point to the compiler versions installed
    through homebrew (e.g. `CC=gcc-9 CXX=g++-9 cmake ..`). 
  - [CMake] >=3.11, download it through your package manager, through
    [`pip`][cmake_pip], or grab the latest [binary][cmake_bin].
  - [GSL], can be installed through the package manager on most unix systems.
    The package is called `libgsl-dev` on Ubuntu and `gsl` most everywhere else
    (e.g. on OpenSUSE/CentOS or homebrew).
  - [Eigen3] >=3.3.0, can be installed through the package manager on most unix
    systems. The package is called `libeigen3-dev` on Ubuntu, `eigen3` on
    OpenSUSE/CentOS and `eigen` in homebrew.


#### Required - Automatic

The recommended versions of the following dependencies are **downloaded
automatically** using the CMake [FetchContent] module.

  - [HiggsBounds] and [HiggsSignals] for checking against collider Higgs data
  - [AnyHdecay] for easy interfacing to the different HDECAY variants

Please see the [FetchContent] documentation on how to use a local version
instead.

#### Optional
  - [MicrOmegas] for dark matter constraints, set the variable
    `MicrOmegas_ROOT_DIR=/Path/to/micromegas/` when calling `cmake ..`. The
    interface is compatible with version `5.0.8/5.0.9` and `5.2.0`.
  - [EVADE] for vacuum stability constraints
  - [BSMPT] for requiring a first order EW phase transition

#### Tests
ScannerS contains a test suite that can be run with `make test` (or, for more
detailed results with `ctest --output-on-failure`). This test suite runs a
number of unit tests and also generates a few valid parameter points in each
model. These functional tests (called `... run`) can take a few minutes to
complete.

## Usage

Compilation creates several executables in the build directory. They provide
command line help explaining all of the arguments.

All executables support two run modes: `scan` and `check`. In `scan` mode the
parameter space is scanned randomly in the given parameter ranges to obtain
viable parameter points. In `check` mode a given set of parameter points is
rechecked against the implemented constraints.

E.g.:

```
./R2HDM output_file.tsv --config example_input/R2HDM_T1.ini scan -n 10
```
will generate 10 valid parameter points with the constraints and parameter
ranges defined in R2HDM_T1.ini. Any of the values in the config file can also be
passed as command line arguments.

```
./R2HDM rechecked_output_file.tsv check output_file.tsv
```

will recheck the resulting points and save any points that still pass the
constraints (this obviously only makes sense if the constraints have changed in
the meantime).

### Constraints

Constraints in ScannerS have three severity levels that can be configured
through config files or command line arguments:

 - **1 = enabled**, the constraint is applied and only parameter points that
   pass the constraint are saved, *default*
 - **0 = ignored**, the constraint is calculated and the result is saved in an
   entry called `valid_{constraintname}`, but points that fail the constraint
   are still kept in the output
 - **-1 = skip**, the constraint is completely ignored, no related data is
   present in the output

A list of the constraints for each model is available in the command line help.

<!-- Links -->
[ScannerS manual]: https://arxiv.org/abs/2007.02985
[online documentation]: https://jonaswittbrodt.gitlab.io/ScannerS/index.html
[GSL]: https://www.gnu.org/software/gsl/
[CMake]: https://cmake.org/
[cmake_pip]: https://pypi.org/project/cmake/
[FetchContent]: https://cmake.org/cmake/help/latest/module/FetchContent.html
[cmake_bin]: https://cmake.org/download/
[Eigen3]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[HiggsBounds]: https://gitlab.com/higgsbounds/higgsbounds
[HiggsSignals]: https://gitlab.com/higgsbounds/higgssignals
[AnyHdecay]: https://gitlab.com/jonaswittbrodt/anyhdecay
[MicrOmegas]: https://lapth.cnrs.fr/micromegas
[EVADE]: https://gitlab.com/jonaswittbrodt/EVADE
[BSMPT]: https://github.com/phbasler/BSMPT
[Jonas Wittbrodt]: mailto:jonas.wittbrodt@desy.de
[contribution guide]: https://gitlab.com/jonaswittbrodt/ScannerS/-/blob/master/CONTRIBUTING.md
