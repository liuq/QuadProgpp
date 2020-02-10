# QuadProg++

A C++ library for Quadratic Programming which implements the [Goldfarb-Idnani active-set dual method](http://www.javaquant.net/papers/goldfarbidnani.pdf). 

At present it is limited to the solution of strictly convex quadratic programs.

Previous versions of the project were hosted on [sourceforge](https://sourceforge.net/projects/quadprog/?source=directory).

## Install

To build the library simply go through the `cmake .; make; make install` cycle.

In order to use it, you will be required to include in your code file the `Array.hh` header, which contains a handy C++ implementation of Vectors and Matrices.

## Contribution

Contributions and bug fixes are welcome.

Copyright (C) 2007-2016 Luca Di Gaspero, [MIT License](LICENSE). 

## References

D. Goldfarb, A. Idnani. A numerically stable dual method for solving strictly convex quadratic programs. *Mathematical Programming* 27 (1983) pp. 1-33.
