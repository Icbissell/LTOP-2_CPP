#ifndef INTERP_1D_H
#define INTERP_1D_H
#include <cmath>
#include <iostream>

class Interp
{
public:
    int n, mm, jsav, corSearch, dj;
    double *xArray, *yArray; //x index vector and y vector to be interpolated on

    Interp();
    Interp(double *xArray, double *yArray, int m, int arrSize);

    double interp(double x); //decide which bracketing algorithm to use and interpolate on found bracket.
    int search(double x); //find j such that x is bracketed in indices [j] -[j+mm-1]
    int track(double x); //find j such that x is bracketed in indices [j] -[j+mm-1], using correlated search

    double virtual primInterp(int offSet, double x) = 0; //primitive linear interpolation method
};

class LinearInterp : public Interp
{

public:
    LinearInterp(double *xArray, double *yArray, int arrSize); //Linear interpolation
    double primInterp(int j, double x);
};

/* BilinInterp takes as inputs size of x dimension, size of y dimension, equally spaced x axis array,
   equally spaced y axis, and matrix to be queried */

class BilinInterp
{
public:
    int m,n; //dimensions of matrix
    double **yMatrix; //matrix to be interpolated
    LinearInterp interpx1, interpx2; //axes for indexing
    BilinInterp(double x1Size, double x2Size, double *x1Array, double *x2Array, double **yMatrix); //bilinear interpolation object

    double interp(double x1p, double x2p); //interpolation function, takes as input point in x and y dimension
};


#endif // INTERP_1D_H

