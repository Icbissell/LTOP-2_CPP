# LTOP2-CPP

These are the source files for LTOP-2, a program for modeling the steady-state orographic precipition field. The program has a number of dependencies, which are all additionally installed with the CMake installation. :

* [Gnuplot](http://www.gnuplot.info/) and [C++ Boost](https://www.boost.org/) are used for creating plots. Both are open source and very tractable. The necessary Boost files are also included with the LTOP2-CMAKE file.

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is used (albeit very briefly) to find the least-squares solution to a system of equations in the base_state() function. 

* [FFTW3](http://www.fftw.org/) is used in all the necessary Fourier transforms. 




