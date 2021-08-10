![Linux Build](https://github.com/ammarhakim/gkylzero/actions/workflows/main.yml/badge.svg)

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

Clone this repo, then type:
```
    make
```
in the top-level directory. To run all unit tests do:
```
    make check
```

You can specify the compiler when running make. For example:
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

# LAPACK and BLAS

Some parts of GkeyllZero rely on presence of LAPACK and BLAS. Very
likely highly optimized builds of these libraries are already
installed on your platform (Framework Accelerator on Mac, MKL on Intel
machines etc). In this case there is no need to install anything and
simply set the appropriate flags to the make command to find/use your
installed BLAS/LAPACK libraries.

However, if you must install these please use the OpenBLAS
implementation. This is a highly optimized library that is
well-maintianed and easy to build and install. There are scripts to
build this already.  From the top-level directory do:
```
cd install-deps
./mkdeps.sh --build-openblas=yes
```

This will install OpenBLAS in the $HOME/gkylsoft directory. Note that
OpenBLAS requires you to have gfortran installed. You are responsible
for installing this on your machine.

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
