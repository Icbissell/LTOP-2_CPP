# LTOP-2_CPP

   By Mark Brandon and Immanuel Bissell of Yale University. 
   
 ![alt text](https://github.com/Icbissell/LTOP-2_CPP/blob/main/misc/Precipitation_Field.png)

## About
   LTOP-2 is a general program to simulate the steady state precipitation field with orographic lifting over a user-defined topography. This repository contains the C++ scripts under 'LTOP-2_CPP' as well as a [CMake](https://cmake.org/) file with two examples (a simple Gaussian hill and a section of the Himalayan shelf) under 'LTOP-2_CMAKE'. Instructions for setting up the CMake script are included in its respective README file. 
   In order to run the example usage program with CMake, you must have cmake on your machine, which can be installed on Linux with:
	
	 sudo apt-get install cmake 
  or macOS:
  
    brew install cmake
   
Plots are generated internally using the open-source [Gnuplot](http://www.gnuplot.info/), which can be installed on Linux with
   
    sudo apt-get install gnuplot 
    
 or macOS:
  
    brew install gnuplot
   
   The program has a number of additional dependencies, which are all installed with the CMake installation:

* [C++ Boost](https://www.boost.org/) is used for creating plots. The libraries are open source and very tractable. The necessary Boost files are included with the LTOP2-CMAKE file.

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is used (albeit very briefly) to find the least-squares solution to a system of equations in the base_state() function. The necessary header files are also included with the LTOP2-CMAKE file.

* [FFTW3](http://www.fftw.org/) is used in all the necessary Fourier transforms. 


## Abstract
   
## Getting the code
   You can download both the CMake file and LTOP-2 source files by navigating to your desired folder and running:
   
         git clone https://github.com/Icbissell/LTOP-2_CPP.git 
