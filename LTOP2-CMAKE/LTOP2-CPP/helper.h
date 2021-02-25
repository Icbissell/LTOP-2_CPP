#ifndef HELPER_H
#define HELPER_H
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <fftw3.h>


#define REAL 0
#define IMAG 1

void fft(int xsize, int ysize, fftw_complex *in, fftw_complex *out); //calculates the fft for a 2d array of size xsize * ysize

void ifft(int xsize, int ysize, fftw_complex *in, fftw_complex *out); /* calculates the inverse fft for a 2d array of size 'xsize' * 'ysize',
                                                                         rescales the array by size 'xsize' * 'ysize' */

double find_max_min(double *inarray, std::string choice, int size); //finds the max or min of a 1d array of size 'size'

double find_max_min_2d(double **inarray, std::string choice, int xsize, int ysize); //finds the max or min of a 2d array of size 'xsize' * 'ysize'

void cum_trapz_integrate(double **arr, double **int_sum, int xsize, int ysize); //calculates the trapezoidal integration down the rows of a matrix (compare to MatLab's cumtrapz)

double mean1(std::vector <double> array, int size);

#endif // HELPER_H
