#include "interp.h"

Interp::Interp(double *xArray, double *yArray, int m, int arrSize){

    this->n = arrSize;
    this->mm = m;
    this->jsav = 0;
    this->corSearch = 0;
    this->xArray = xArray;
    this->yArray = yArray;
    this->dj = std::min(1,(int)pow((double)n,0.25));
}

double Interp::interp(double x){
    int offSet = corSearch ? track(x) : search(x);
    return primInterp(offSet,x);
}


int Interp::search(double x){
    int up, mid, low; //bracketing indices
    if (n < 2 || mm < 2 || mm > n) {
        std::cerr << "locate size error" << std::endl;
        while (true) {}
    }

    bool ascnd = (xArray[n-1] >= xArray[0]); //true if ascending order of table, false otherwise.

    low=0; //initialize lower
    up = n-1; //and upper limits.

    /* adjust limits to find bracket */
    while (up-low > 1) {
        mid = (up+low) >> 1;
        if (x >= xArray[mid] == ascnd){
            low=mid;
        }
        else{
            up=mid;
        }
    }

    corSearch = abs(low-jsav) > dj ? 0 : 1; //update correlated search decisor
    jsav = low;
    return std::max(0,std::min(n-mm,low-((mm-2)>>1)));
}

int Interp::track(double x){
    int low = jsav, mid, up, inc = 1000;

    if (n < 2 || mm < 2 || mm > n) {
        std::cerr << "locate size error" << std::endl;
        while (!false) {}
    }

    bool ascnd = (xArray[n-1] >= xArray[0]); //True if ascending order of table, false otherwise.

    if (low < 0 || low > n-1) {
        low=0;
        up=n-1;
    }

    else{
        if (x >= xArray[low] == ascnd) { //track up:
            while(true) {
                up = low + inc;
                if (up >= n-1) { // Off end of table.
                    up = n-1; // Found bracket.
                    break;
                }
                else if(x < xArray[up] == ascnd){
                    break;
                }
                else{
                    low = up;
                    inc += inc;
                }
            }
        }
        else{ // track down
            up = low;
            while(true){
                low = low - inc;
                if (low <= 0) {
                    low = 0;
                    break;
                }
                else if (x >= xArray[low] == ascnd) {
                    break;
                }
                else{
                    up = low;
                    inc += inc;
                }
            }
        }
    }

    while(up-low > 1){
        mid = (up+low) >> 1;
        if (x >= xArray[mid] == ascnd){
            low=mid;
        }
        else{
            up=mid;
        }
    }

    corSearch = abs(low-jsav) > dj ? 0 : 1;
    jsav = low;
    return std::max(0, std::min(n-mm, low-((mm-2)>>1)));

}


LinearInterp::LinearInterp(double *xArray, double *yArray, int arrSize) : Interp(xArray, yArray, 2, arrSize){} //constructor

double LinearInterp::primInterp(int j, double x){

    if(xArray[j] == xArray[j+1]){ //If same value, return value
        return yArray[j];
    }

    else{
        return yArray[j] + ((x-xArray[j])/(xArray[j+1]-xArray[j]))*(yArray[j+1]-yArray[j]); //weigh and return value
    }

}

BilinInterp::BilinInterp(double x1Size, double x2Size, double *x1Array, double *x2Array, double **yMatrix) :
    m(x1Size), n(x2Size), yMatrix(yMatrix), interpx1(x1Array,x1Array, x1Size), interpx2(x2Array,x2Array, x2Size) {} //constructor

double BilinInterp::interp(double x1p, double x2p){

    int i,j;
    double yy, t, u;
    i = interpx1.corSearch ? interpx1.track(x1p) : interpx1.search(x1p);
    j = interpx2.corSearch ? interpx2.track(x2p) : interpx2.search(x2p);

    t = (x1p-interpx1.xArray[i])/(interpx1.xArray[i+1]-interpx1.xArray[i]);
    u = (x2p-interpx2.xArray[j])/(interpx2.xArray[j+1]-interpx2.xArray[j]);

    yy = (1.-t)*(1.-u)*yMatrix[i][j] + t*(1.-u)*yMatrix[i+1][j] + (1.-t)*u*yMatrix[i][j+1] + t*u*yMatrix[i+1][j+1];
    return yy;

}

