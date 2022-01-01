![Linux Build](https://github.com/ammarhakim/gkylzero/actions/workflows/main.yml/badge.svg)
![Mac Build](https://github.com/ammarhakim/gkylzero/actions/workflows/osx-build.yml/badge.svg)

# About

This is the GkeyllZero layer of the Gkeyll code. The name is
pronounced as in the book "The Strange Case of Dr. Jekyll and
Mr. Hyde" which is required reading for all members of the Gkeyll
Team. GkeyllZero is written in C and some C++. GkeyllZero is developed
at Princeton Plasma Physics Laboratory (PPPL) and is copyrighted
2020-2021 by Ammar Hakim and the Gkeyll Team.

**NOTE**: Though the code is called "Gkeyll" and "GkeyllZero" we use
the prefix "gkyl" in the source code itself. This may be confusing at
first but one gets used to it.

Documentation is available at http://gkeyll.rtfd.io.

# Install

GkeyllZero has just a few dependencies: OpenBLAS and SuperLU. If you
are on a cluster these libraries are likely already installed. In that
case you should use the configure to specify the location of the
include and full path to the static libraries for these. See
instructions below. If you are on your own machine please use the
mkdeps.sh shell-script to install them and then follow the
instructions below.

Clone this repo and then optinally run the configure script to set the
compiler you wish to use and the various paths to the
dependencies. Once you are done with installing/specifying the compiler 
and dependencies simply type:
```
    make
```
in the top-level directory. To run all unit tests do:
```
    make check
```

If you do not run configure you can also specify the compiler when
running make. For example:
```
    make CC=icc
```
to build with the Intel C compiler. If you use the NVIDIA nvcc
compiler then the CUDA specific parts of the code will be built:
```
    make CC=nvcc
```
Note that if your machine has more than one core (a highly likely
situation) you can run make in parallel, for example:
```
    make -j
```
will use several cores while compiling the code, and can be
potentially faster on most machines.

The unit and regression test executables are written in the
`build/unit` and `build/regression` directories.

If you want to use the code as a library you should install it:
```
  make install
```

Note that GkeyllZero is meant to be used as a *library*. You can use
it to create your own "app" for your particular problem. See that
various "app_*.c" files for examples. Full documentation is available
on the RTFD website linked above.

# LAPACK/BLAS and SuperLU

Some parts of GkeyllZero rely on presence of LAPACK and BLAS. Its
likely that highly optimized builds of these libraries are already
installed on your platform (Framework Accelerator on Mac, MKL on Intel
machines etc). In this case there is no need to install anything and
simply set the appropriate configure flags to the make command to
find/use your installed BLAS/LAPACK libraries. Note that on Macs you
do not need to do anything special to use LAPACK/BLAS as it comes with
the developer tools.

However, if you must install these please use the OpenBLAS
implementation. This is a highly optimized library that is
well-maintianed and easy to build and install. There are scripts to
build this already.  From the top-level directory do:
```
cd install-deps
./mkdeps.sh --build-openblas=yes --build-superlu=yes
```

This will install OpenBLAS and SuperLU in the $HOME/gkylsoft
directory. Note that OpenBLAS requires you to have gfortran
installed. You are responsible for installing this on your machine.

# Using precompiled OpenBLAS

If you can't (do not wish to) build the OpenBLAS library you can
instead install the pre-built ones that are available on
Ubuntu. (Similar libraries will be available on other Linux
distros). First, install the libraries as:
```
sudo apt-get install libopenblas-dev
sudo apt-get install liblapacke-dev
```

Then, configure the build system:
```
./configure --lapack-inc=/usr/include/x86_64-linux-gnu/openblas-pthread --lapack-inc="/usr/lib/x86_64-linux-gnu/liblapacke.a /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblas.a /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.a"
```

For other Linux distros you will need to adjust the paths shown above.

# Developing for GkeyllZero

Out goal is to keep GkeyllZero as simple and dependency free as
possible. Some dependencies are unavoidable like MPI and linear
algebra libraries. However, we must avoid an exponentially increasing
dependency chain. Another goal is that GkeyllZero itself should be
pure modern (C99) C. Some parts of the code need C++ (nvcc is a C++
compiler for CUDA) but the core code itself should be in C.

Developing in C (and C++) requires very strong focus and
discipline. **Please consult** https://en.cppreference.com/w/ for
standards documentation for these languages and their
libraries. **Please use valgrind** to make sure all code is "valgrind
clean". Pay attention to all compiler warnings.

Most importantly, **please internalize and follow** the programming
philosophy outlined in the document ["A Minimalist Approach to
Software"](http://ammar-hakim.org/minimalist-software.html).

# License

GkeyllZero and Gkeyll use the BSD-3 license. See below.

Copyright 2020-2021 Ammar Hakim and the Gkeyll Team

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the
   distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
