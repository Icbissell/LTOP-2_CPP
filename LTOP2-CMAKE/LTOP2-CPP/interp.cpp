#include "interp.h"

BaseInterp::BaseInterp(double *xArray, double *yArray, int m, int arrSize){
    this->n = arrSize;
    this->mm = m;
    this->jsav = 0;
    this->cor = 0;
    this->xArray = xArray;
    this->yArray = yArray;
    this->dj = std::min(1,(int)pow((double)n,0.25));
}

double BaseInterp::interp(double x){
    int jlo = cor ? hunt(x) : locate(x);
    return rawInterp(jlo,x);
}


int BaseInterp::locate(double x){
    int ju, jm, jl;
    if (n < 2 || mm < 2 || mm > n) {
        std::cerr << "locate size error" << std::endl;
        while (true) {}
    }

    bool ascnd = (xArray[n-1] >= xArray[0]); //true if ascending order of table, false otherwise.

    jl=0; //initialize lower
    ju = n-1; //and upper limits.

    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (x >= xArray[jm] == ascnd){
            jl=jm;
        }
        else{
            ju=jm;
        }
    }

    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
}


/* Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
increasing or decreasing. The returned value is not less than 0, nor greater than n-1. */

int BaseInterp::hunt(double x){
    int jl = jsav, jm, ju, inc = 1000;

    if (n < 2 || mm < 2 || mm > n) {
        std::cerr << "locate size error" << std::endl;
        while (!false) {}
    }

    bool ascnd = (xArray[n-1] >= xArray[0]); //True if ascending order of table, false otherwise.

    if (jl < 0 || jl > n-1) {
        jl=0;
        ju=n-1;
    }

    else{
        if (x >= xArray[jl] == ascnd) { //Hunt up:
            while(true) {
                ju = jl + inc;
                if (ju >= n-1) { // Off end of table.
                    ju = n-1; // Found bracket.
                    break;
                }
                else if(x < xArray[ju] == ascnd){
                    break;
                }
                else{
                    jl = ju;
                    inc += inc;
                }
            }
        }
        else{ // Hunt down
            ju = jl;
            while(true){
                jl = jl - inc;
                if (jl <= 0) {
                    jl = 0;
                    break;
                }
                else if (x >= xArray[jl] == ascnd) {
                    break;
                }
                else{
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }

    while(ju-jl > 1){
        jm = (ju+jl) >> 1;
        if (x >= xArray[jm] == ascnd){
            jl=jm;
        }
        else{
            ju=jm;
        }
    }

    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return std::max(0, std::min(n-mm, jl-((mm-2)>>1)));

}


LinearInterp::LinearInterp(double *xArray, double *yArray, int arrSize) : BaseInterp(xArray, yArray, 2, arrSize){}

double LinearInterp::rawInterp(int j, double x){

    if(xArray[j] == xArray[j+1]){ //Table is defective, but we can recover.
        return yArray[j];
    }

    else{
        return yArray[j] + ((x-xArray[j])/(xArray[j+1]-xArray[j]))*(yArray[j+1]-yArray[j]);
    }

}

BilinInterp::BilinInterp(double x1Size, double x2Size, double *x1Array, double *x2Array, double **yMatrix) :
    m(x1Size), n(x2Size), yMatrix(yMatrix), x1terp(x1Array,x1Array, x1Size), x2terp(x2Array,x2Array, x2Size) {}

double BilinInterp::interp(double x1p, double x2p){

    int i,j;
    double yy, t, u;
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);

    t = (x1p-x1terp.xArray[i])/(x1terp.xArray[i+1]-x1terp.xArray[i]);
    u = (x2p-x2terp.xArray[j])/(x2terp.xArray[j+1]-x2terp.xArray[j]);

    yy = (1.-t)*(1.-u)*yMatrix[i][j] + t*(1.-u)*yMatrix[i+1][j] + (1.-t)*u*yMatrix[i][j+1] + t*u*yMatrix[i+1][j+1];
    return yy;

}

