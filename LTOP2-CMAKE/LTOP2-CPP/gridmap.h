#ifndef GRIDMAP_H
#define GRIDMAP_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <fftw3.h>
#include <ctime>
#include "helper.h"
#include "gnuplot-iostream.h"
#include "interp.h"
#include "farfield.h"

#define REAL 0
#define IMAG 1
#define sind(x) (sin(fmod((x),360) * M_PI / 180.))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180.))

class GridMap
{
public:

    const double g = 9.81; //constant acceleration due to gravity
    const double TC2K = 273.15; //conversion from Celsius to Kelvin, [K]
    const double U = 10.0; //wind speed
    const double EPS = 2.2204e-16; //machine epsilon
    const double kappa = 0.0; //eddy diffusion (m/s^2)
    const double radiusEarth = 6371e3;
    const double omega = 7.2921e-5; //rotation rate of earth
    double azimuth = 90; //wind azimuth degrees
    double tauC = 1000.0;  //characteristic time for cloud water, [s]
    double tauF = 1000.0;  //characteristic time for precipitation fallout, [s]
    double Ns = 0.005; //saturated bouyancy frequency, [rad/s]
    double fC = 0;  //coriolis frequency, [rad/s]
    double T0 = TC2K + 15.0; //sea level Temperature [k]
    double gammaEnv;//envirnomental lapse rate, [K/m]
    double gammaSat; //saturated lapse rate, [K/m]
    double hRho; //density scale height
    double rhoS0; //saturated density at base [kg/m^3]
    double gammaRatio; //ratio of gammaSat/gammaEnv
    double hV; //scale height for saturated vapor (scalar, m)
    double rho0; //total air density at z=0 (scalar, kg/m^3)
    double fP = 0.97; //fraction of precipitation leaving the model

    /*Convert needed variables to complex for later calculations */

    std::complex <double> img = std::complex<double>(0.0,1.0); //complex number 0 + 1i, for complex calculations
    std::complex <double> comp_one = std::complex<double>(1.0,0.0); // 1, for complex calculations

    /* Set up dimensions for topography */

    int n0; //dimension of topography
    int n1; //dimension of topography
    int n; //points in atmospheric profile
    int nXPad; //padded dimensions of topography
    int nYPad; //padded dimensions of topography
    int rowLength = n0*n1;
    int bigRowLength = nXPad * nYPad;
    int sLength; //dimension of azimuth transfomed topography
    int tLength; //dimension of azimuth transfomed topography
    int nSPad; //padded dimension of azimuth transfomed topography
    int nTPad; //padded dimension of azimuth transfomed topography
    int stRowLength; //total elements in sLength * tLength array
    int padstRowLength; //total elements in nSPad * nTPad array
    double dX; //node spacing on grid in x-direction
    double dY; //node spacing on grid in y-direction
    double dkS; //spacing of kS wavenumber
    double dkT; //spacing of kT wavenumber
    double dS; //node spacing on grid in s-direction
    double dT; //node spacing on grid in t-direction
    double mostPrec;
    double mostPrec_ma; //maximum precipitation in meters/year
    double mostPrec_mmhr; //maximum precipitation in mm/hr
    double max_el;

    /* Set arrays for calculation */
    double *kS; //wavenumbers in kS direction
    double *kT; //wavenumbers in kT direction
    double *Sxy; //s coordinates in xy frame of reference
    double *Txy; //t coordinates in xy frame of reference
    double *Xst; //x coordinates in st frame of reference
    double *Yst; //y coordinates in st frame of reference
    double *T; //temperature profile
    double *zBar; //height of isosurface relative to topography
    double *gammaSatArr; //saturated lapse rate profile
    double *gammaEnvArr; //environmental lapse rate profile
    double **pWind; //precipitation grid in wind frame of reference
    double **zInterpolate; //topography in wind frame of reference
    double **pGrid; //precipitation grid in geographical frame of reference
    double **grid; //topographical grid

    std::vector<double> s; //s coordinates
    std::vector<double> t; //t coordinates
    std::vector<double> denominator; //denominators of vertical wavenumbers
    std::vector<std::vector<double>> topoVec; //vector for transfering topography
    std::vector<std::complex<double>> kZ; //wavenumbers in Z(vertical) direction
    std::vector<std::complex<double>> hHatPair;
    std::vector<std::complex<double>> GVHat; //conversion from water vapor to cloud water
    std::vector<std::complex<double>> GCHat; //transport and release of cloudwater
    std::vector<std::complex<double>> GFHat; //the fall of precipitation to the ground
    std::vector<double> pStarWindVec; //reduced precitation (no moisture balance)
    std::vector<double> QFStarWindVec; //column density for falling precipitation
    std::vector<double> QCStarWindVec; //column density for cloud water

    fftw_complex *hHat; //globally used Fourier transform of topography

    fftw_plan plan_r2c, plan_c2r;  //plans for fourier transforms

    GridMap(double azimuth); //constructor

    ~GridMap(); //destructor

    void get_input(); //read topography from the .csv file

    void calc_fc(); //calculate Coriolis frequency from latitude and longitude arrays

    void fill_base_state(double Ns, double T0); //add constants and profiles calculated in farfield.h to gridmap class

    void calc_tauF(); /* calculation of time constant for fallout. Invokes a calculation of the 258 K isosurface, which
                        marks the midpoint of the 268 - 248 K range for freezing in the atmosphere
                        (WBF zone, Cias and Jouzel, 1994).*/

    void fill_grid_topo(); //fill grid with topographical values from DEM file

    void interp_axes(int xsize, int ysize, double angle); /* calculates axes (s,t,x,y) and coordinates (Xst, Yst, Sxy, Txy) for original and rotated topography
                                                            as well as grid spacing dS and dT */
    void rotate(); //rotates topography by azimuth

    void reverse_rotation(double** inGrid, double** outGrid); //rotates topography back to original orientation

    void fill_waveNumber(); //calculates kT and kS wavenumbers

    void topo_to_wave_domain(); //converts Topography into fourier space/ wave domain

    void treat_singularities(); //treats singularities caused by abs(U*kS)==abs(fC) by setting the vertical velocities to 0

    void calc_kZ(); //calculates the vertical wave number

    void calc_response_functions(); /* calculates three response functions, GVHat( conversion from water vapor to cloud water),
                                       GCHat(transport and release of cloud water), and GFHat(fall of hydrometeors to the ground) */

    void convolution_functions(); /* calculates the pStarHat, the reference (reduced) precipition rate,
                                    and column densities QPStarHat and QCStarHat in the wave domain, and transforms to space domain
                                    QPStarWindVec, pStarWindVec, and QCStarWindVec carry these values to next function*/

    void final_prec_calc(); //calculates precipitation rate in rotated frame of reference

    void calc_most(); //calculate the maximum precipitation

    void LTOP_calc(); //container function for LTOP Calculations

    void make_plots(); /* creates heat plots of topography and precipitation rate
                        using gnuplot and a gnuplot pipeline (gnuplot-iostream.h) */
};

#endif // GRIDMAP_H

