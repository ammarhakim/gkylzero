![Build Dist](https://github.com/ammarhakim/gkylzero/actions/workflows/main.yml/badge.svg)

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

You can install the code in two ways, depending on what you want to do
with it. First, if you want to use the code as a library or run
simulations you should be using a released *deployment version*
shipped as a zip file. Alternately, clone this repo and run:
```
  ./mkdist.sh
  cd gkylzero-dist
  make install
```

GkeyllZero has minimal dependencies (none at all if you want to run in
serial; MPI for parallel and BLAS/LAPACK for linear algebra for some
solvers). You can set the compiler you want to use as:
```
    make CC=icc install 
```

You can run the unit tests from the deployment directory as:
```
    make check
```

Note that GkeyllZero is meant to be used as a *library*. You can use
it to create your own "app" for your particular problem. See that
various "app_*.c" files for examples. Full documentation is available
on the RTFD website linked above.

# Developing for GkeyllZero

Out goal is to keep GkeyllZero as simple and dependency free as
possible. Some dependencies are unavoidable like MPI, linear algebra
and FFT libraries. However, where possible we would like to avoid an
exponentially increasing dependency chain. Another goal is that
GkeyllZero itself should be pure modern (C99) C. Some tools used in
code generation need C++ (GiNaC CAS and Swig wrapper generator, for
example) but the *generated* code itself will be in C.

Developing in C (and C++) requires very strong focus and
discipline. Please consult https://en.cppreference.com/w/ for
standards documentation for these languages and their
libraries. Please use valgrind to make sure all code is "valgrind
clean". Pay attention to all compiler warnings.

Most importantly, *please internalize and follow* the programming
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
