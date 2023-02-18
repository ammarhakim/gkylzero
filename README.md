![Linux Build](https://github.com/ammarhakim/gkylzero/actions/workflows/main.yml/badge.svg)
![Mac Build](https://github.com/ammarhakim/gkylzero/actions/workflows/osx-build.yml/badge.svg)

# About

This is the GkeyllZero layer of the Gkeyll code. The name is
pronounced as in the book "The Strange Case of Dr. Jekyll and
Mr. Hyde" which is required reading for all members of the Gkeyll
Team. GkeyllZero is written in C and some C++. GkeyllZero is developed
at Princeton Plasma Physics Laboratory (PPPL) and is copyrighted
2020-2022 by Ammar Hakim and the Gkeyll Team.

**NOTE**: Though the code is called "Gkeyll" and "GkeyllZero" we use
the prefix "gkyl" in the source code itself. This may be confusing at
first but one gets used to it.

Documentation is available at http://gkeyll.rtfd.io.

# Install

GkeyllZero has just a few dependencies: OpenBLAS and SuperLU. If you
are on a cluster these libraries are likely already installed. In that
case you should use the configure script to specify the location of
the include and full path to the static libraries for these. See
configure help.
```
   ./configure --help
```

# Building dependencies

Some parts of GkeyllZero rely on presence of LAPACK and BLAS. Note
that on Macs you do not need to do anything special to use LAPACK/BLAS
as it comes with the developer tools. GkeyllZero also relies on
SuperLU for sparse, direct linear solvers.

To install these yourself please use the mkdep.sh script to do so.
From the top-level directory do:
```
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes
```

This will install OpenBLAS and SuperLU in the $HOME/gkylsoft
directory. Note that OpenBLAS **requires you to have gfortran**
installed. Likewise, SuperLU **requires you to have cmake** installed.
**You** are responsible for installing this on your machine.

# Building Using Machine Files

We have provided a set of "machine files" to ease the build
process. These are stored in the machines directory. For example, to
build on Traverse please run
```
./machines/mkdeps.traverse.sh
./machines/configure.traverse.sh
```
After this is completed then just type:
```
make -j #
```
where # is the number of cores you wish to use (if working on the login
node of a cluster it is preferable that you use fewer than the total number
of cores available, otherwise you will slow down the login node and irk
other users).

Note: On Traverse and Stellar-amd you need to load cudatoolkit/11.6. 

# Building on your own machine

If you are on your own machine please use the mkdeps.sh shell-script
in the install-deps directory
```
   cd install-deps
   ./mkdeps.sh --build-openblas=yes --build-superlu=yes
```
to install dependencies. We have also included a number of pre-packaged
configure.[machine].sh and mkdeps.[machine].sh scripts in the machines directory.

Clone this repo and then optionally run the configure script to set the
compiler you wish to use and the various paths to the
dependencies. Once you are done with installing/specifying the compiler 
and dependencies simply type:
```
    make -j #
```
in the top-level directory (# is the number of cores, see previous comment
on this). To run all unit tests do:
```
    make check
```

If you do not run configure you can also specify the compiler when
running make. For example:
```
    make CC=icc -j #
```
to build with the Intel C compiler (# is the number of cores, see previous
comment on this). If you use the NVIDIA nvcc compiler then the CUDA
specific parts of the code will be built:
```
    make CC=nvcc -j #
```
(# is the number of cores, see previous comment on this). The
unit and regression test executables are written in the
`build/unit` and `build/regression` directories.

If you want to use the code as a library (e.g. for use by
[gkyl](https://github.com/ammarhakim/gkyl/)) you should install it:
```
  make install
```

Note that GkeyllZero is meant to be used as a *library*. You can use
it to create your own "app" for your particular problem. See that
various "app_*.c" files for examples. Full documentation is available
on the RTFD website linked above.


# Developing for GkeyllZero

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
