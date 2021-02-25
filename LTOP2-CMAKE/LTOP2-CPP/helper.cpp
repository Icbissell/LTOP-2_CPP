#include "helper.h"

void fft(int xsize, int ysize, fftw_complex *in, fftw_complex *out){

    /* create a plan for the Fourier transform */
    fftw_plan plan = fftw_plan_dft_2d(xsize, ysize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    /* execute and cleanup plan */
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
}

void ifft(int xsize, int ysize, fftw_complex *in, fftw_complex *out){

    /* create a plan for the inverse Fourier transform */
    fftw_plan plan = fftw_plan_dft_2d(xsize, ysize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* execute and cleanup plan */
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();

    /* values are unnormalized (scaled by xsize*ysize), so we now divide by the size of array to normalize */
    int scale = xsize*ysize;

    for(int i = 0 ; i < scale; i ++){
        out[i][REAL] /= scale;
        out[i][IMAG] /= scale;
    }
}

double find_max_min(double *inarray, std::string choice, int size){

    /* set reference value to the first in array */
    double bound = inarray[0];

    /* search through array. If a larger value is found, update bound */
    if(choice == "max"){
        for(int i =0; i < size; i++){
            if(bound < inarray[i]){
                bound = inarray[i];
            }
        }
    }

    /* search through array. If a smaller value is found, update bound */
    if(choice == "min"){
        for(int i =0; i < size; i++){
            if(bound > inarray[i]){
                bound = inarray[i];
            }
        }
    }

    return bound;
}

double find_max_min_2d(double **inarray, std::string choice, int xsize, int ysize){

    /* set reference value to the first in array */
    double bound = inarray[0][0];

    /* search through array. If a larger value is found, update bound */
    if(choice == "max"){
        for(int i =0; i < xsize; i++){
            for(int j =0; j < ysize; j++){
                if(bound < inarray[i][j]){
                    bound = inarray[i][j];
                }
            }
        }
    }

    /* search through array. If a smaller value is found, update bound */
    if(choice == "min"){
        for(int i =0; i < xsize; i++){
            for(int j =0; j < ysize; j++){
                if(bound > inarray[i][j]){
                    bound = inarray[i][j];
                }
            }
        }
    }

    return bound;

}

void cum_trapz_integrate(double **arr, double **int_sum, int xsize, int ysize){

    /* Set elements in first row of integrated matrix to 0 (a single element has no area) */
    for(int i = 0; i < ysize; i++){
        int_sum[0][i] = 0;
    }

    /* Move down columns of array calculating a trapazoidal integration, i.e. (hi+h(i+1))/(x(i+1)-xi)*/
    for(int i = 1; i < xsize; i++){
        for(int j = 0; j < ysize; j++){
            double stepTrap = (arr[i][j] + arr[i-1][j])/2.;
            int_sum[i][j] = int_sum[i-1][j] + stepTrap;
        }
    }

}

double mean1(std::vector <double> array, int size){

    double arrSum = 0;
    for(int i = 0; i < size; i++){
        arrSum += array[i];
    }

    return arrSum/size;

}

