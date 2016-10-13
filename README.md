# QuadProg++

A C++ library for Quadratic Programming which implements the [http://www.javaquant.net/papers/goldfarbidnani.pdf](Goldfarb-Idnani active-set dual method). 

At present it is limited to the solution of strictly convex quadratic programs.

Previous versions of the project were hosted on [https://sourceforge.net/projects/quadprog/?source=directory](sourceforge).

## Install

To build the library simply go through the `cmake .; make; make install` cycle.

In order to use it, you will be required to include in your code file the `Array.hh` header, which contains a handy C++ implementation of Vectors and Matrices.

## Contribution

Contributions and bug fixes are welcome.

Copyright (C) 2001-2016 Luca Di Gaspero, MIT License. 