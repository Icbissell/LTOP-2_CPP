#ifndef INTERP_1D_H
#define INTERP_1D_H
#include <cmath>
#include <iostream>

class BaseInterp
{
public:
    int n, mm, jsav, cor, dj;
    double *xArray, *yArray;

    BaseInterp();
    BaseInterp(double *xArray, double *yArray, int m, int arrSize);

    double interp(double x);

    int locate(double x);

    int hunt(double x);

    double virtual rawInterp(int jlo, double x) = 0;
};

class LinearInterp : public BaseInterp
{

public:
    LinearInterp(double *xArray, double *yArray, int arrSize);

    double rawInterp(int j, double x);

};

/* BilinInterp takes as inputs size of x dimension, size of y dimension, equally spaced x axis array,
   equally spaced y axis, and matrix to be queried */

class BilinInterp
{
public:
    int m,n;
    double **yMatrix;
    LinearInterp x1terp, x2terp;
    BilinInterp(double x1Size, double x2Size, double *x1Array, double *x2Array, double **yMatrix);

    double interp(double x1p, double x2p);
};


#endif // INTERP_1D_H

