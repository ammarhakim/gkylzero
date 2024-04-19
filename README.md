![Linux Build](https://github.com/ammarhakim/gkylzero/actions/workflows/main.yml/badge.svg)
![Mac Build](https://github.com/ammarhakim/gkylzero/actions/workflows/osx-build.yml/badge.svg)

# About

This is the GkeyllZero layer of the Gkeyll code. The name is
pronounced as in the book "The Strange Case of Dr. Jekyll and
Mr. Hyde" which is required reading for all members of the Gkeyll
Team. GkeyllZero is written in C (and some C++ for GPUs).
GkeyllZero is developed at the Princeton Plasma Physics Laboratory (PPPL)
and is copyrighted 2020-2023 by Ammar Hakim and the Gkeyll Team.

**NOTE**: Though the code is called "Gkeyll" and "GkeyllZero" we use
the prefix "gkyl" in the source code itself. This may be confusing at
first but one gets used to it.

Documentation is available at http://gkeyll.rtfd.io.

# Installation

Installing GkeyllZero consists of three steps: building dependencies,
configuring, and compiling. How you do the first two of these steps varies
depending on whether you are installing on a computer we already have
installation scripts ("machine files") for, or if you are installing in
a new computer. Since the former is the most common we describe that first.

## On a known computer using machine files

We provide a set of "machine files" to ease the build process.
These are stored in the machines directory. For example, to
build on Traverse please run
```
./machines/mkdeps.traverse.sh
./machines/configure.traverse.sh
```
At this point you may need to manually load environment modules by hand
if you haven't done it yet; for example on Traverse and Stellar-amd you
need to ```module load cudatoolkit/11.6``` before proceeding. Next,
compile the code with:
```
make -j # install
```
where # is the number of cores you wish to use (if working on the login
node of a cluster it is preferable that you use fewer than the total number
of cores available, otherwise you will slow down the login node and irk
other users and admins).

The ```machine``` files default to installing in $HOME/gkylsoft/, and to
searching for dependencies there too. If you wish to install somewhere else
you'll need to indicate the install directory in ```--prefix``` in both
machine files, and make sure the paths to the dependencies are correct
in the ```machines/configure.<machine name>.sh``` file.

## On a new computer (no machine files available)

When installing on your own computer, or on a cluster we don't yet have
machine files for, you could choose to build machine files using existing
ones as guides, or you could do each step manually. We expand on each of
these steps below.

### Specifying/building dependencies

GkeyllZero has a few dependencies (e.g. OpenBLAS, SuperLU, MPI (optional)).
If you are on a cluster some or all of these libraries are likely already
installed. To install dependencies (either all or the ones your cluster
doesn't have) use the ```mkdeps.sh``` script in ```install-deps``` as
follows:
```
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes
```
In this example we opted to build OpenBLAS and SuperLU, assuming
neither is available in this new computer. This will install OpenBLAS and
SuperLU in the ```$HOME/gkylsoft``` directory; in order to install in a
different directory one must specify it with ``--prefix=``.

Note that OpenBLAS **requires you to have gfortran**
installed. Likewise, SuperLU **requires you to have cmake** installed.
**You** are responsible for installing this on your machine. Also keep
in mind that on Apple Macs you do not need to do anything special to use
LAPACK/BLAS as it comes with the developer tools.

At this step, you can download ADAS data for neutral interactions by adding
the flag ```--use-adas=yes``` after the other dependencies to be
installed.

### Configuring

The (optional) configure step is used when we need to specify the use of
specific compilers, dependency paths or other options. For example, if
a dependency is already available in your computer/cluster, you should
use the configure script to specify the location of the include and full
path to the static libraries for these dependencies. See configure help.
```
./configure --help
```
You can also see machine files in ```machines/``` for examples of how
this script is used.


### Compiling

Once you are done with installing/specifying the compiler and
dependencies simply type:
```
 make -j #
```
in the top-level directory (# is the number of cores, see previous
comment on this).

If you do not run configure (see previous section/step) you can also
specify the compiler when running make. For example:
```
 make CC=icc -j #
```
to build with the Intel C compiler (# is the number of cores, see previous
comment on this). If you use the NVIDIA nvcc compiler then the CUDA
specific parts of the code will be built:
```
make CC=nvcc -j #
```
(# is the number of cores, see previous comment on this).

If you want to use the code as a library (e.g. for use by
[gkyl](https://github.com/ammarhakim/gkyl/)) you should install it:
```
make install
```
or run compilation and library installation in one step as
```
make -j # install
```

Note that GkeyllZero is meant to be used as a *library*. You can use
it to create your own "app" for your particular problem. See that
various "app_*.c" files for examples. Full documentation is available
on the RTFD website linked above.

In order to test that your installation worked, you can compile one of
the unit or regression tests and run it. See instructions for that below.

# Developing for GkeyllZero

Built-in tests
--------------

GkeyllZero has built in unit (```/unit```) and regression (```/regression```) tests
that run automatically by our CI system or can be run manually by
developers. These tests are not compiled by the simple ```make```
command used to compile the code (see previous sections). Instead the
developer needs to build specific targets. For example we can build
the unit tests for gkyl_array on a CPU with
```
make build/unit/ctest_array
```
or build the same unit test on a NVIDIA GPU-accelerated node with
```
make cuda-build/unit/ctest_array
```
These produce an executable in the ```build``` or ```cuda-build``` directory,
respectively, that can be run. For example, to run the ```ctest_array``` executable
we just compiled use
```
./cuda-build/unit/ctest_array
```
Alternatively you could choose to build all unit tests with
```
make -j # unit
```
(# is the number of cores, see previous comment on this). A similar procedure
should be followed to compile and run regression tests in ```/regression```.
Regression test executables have a number of run-time options which you can list
with the `-h` flag, i.e. ```<executable> -h```.

You may also compile and run all the unit tests with one command:
```
make -j # check
```

If you wish to run a regression test from a different directory (e.g. scratch
in a cluster), copy the ```Makefile``` in the ```$HOME/gkylsoft/share``` directory
to your desired directory, change ```rt_twostream``` for the name of your test
there in (best to use search & replace), and run the makefile with ```make -j #```,
where ```#``` is a responsibly chosen number of cores.

Automatic regression testing
----------------------------

There is also an automatic regression testing system implemented within GkeyllZero,
included as part of our larger automatic CI framework, which can be built in serial
(for instance, for the gyrokinetics system) with:
```
make build/ci/gk_regression
```
or, in parallel, with:
```
make build/ci/gk_regression_parallel
```
and then run with:
```
./build/ci/gk_regression
```
or:
```
./build/ci/gk_regression_parallel
```
or (for instance, for the moment app) with:
```
make build/ci/moment_regression
```
or, in parallel, with (note that the parallel regression system requires you to have
built GkeyllZero using CUDA and MPI):
```
make cuda-build/ci/moment_regression_parallel
```
and then run with either:
```
./build/ci/moment_regression
```
or:
```
./cuda-build/ci/moment_regression_parallel
```
as necessary. By default, running the automatic regression testing system
presents the user with an interactive list of numerical options (1 to run the
full regression suite, 2 to view all regression results, 3 to run a specific
regression test, 4 to view a specific regression result, 5 to (re)generate
all accepted results, and 6 to (re)generate a specific accepted result).
However, these options can also be passed in as command line arguments, so that:
```
./build/ci/gk_regression 3 4
./build/ci/moment_regression 3 4
```
will specifically run regression test number 4 in serial, with no additional
interaction required from the user. The complete automatic CI framework can be
built in a single step using:
```
make ci
```

Development philosophy
---------------------

Out goal is to keep GkeyllZero as simple and dependency free as
possible. Some dependencies are unavoidable like MPI and linear
algebra libraries. However, we must avoid an exponentially increasing
dependency chain. Another goal is that GkeyllZero itself should be
pure modern (C99/C11) C. Some parts of the code need C++ (nvcc is a
C++ compiler for CUDA) but the core code itself should be in C.

Developing in C (and C++) requires very strong focus and
discipline. **Please consult** https://en.cppreference.com/w/ for
standards documentation for these languages and their
libraries. **Please use valgrind** to make sure all code is "valgrind
clean". Pay attention to all compiler warnings.

Most importantly, **please internalize and follow** the programming
philosophy outlined in the document ["A Minimalist Approach to
Software"](https://www.ammar-hakim.org/sj/pn/pn0/pn0-minimalism.html).