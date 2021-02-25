#ifndef FARFIELD_H
#define FARFIELD_H
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include "Eigen/Dense"

class FarField {

public:
    double hRho; //scale height for total air density (scalar, m)
    double rho0; //total air density at z=0 (scalar, kg/m^3)
    double Ns; //saturated bouyancy frequency
    double T0; //sea-level temperature (scalar, K)
    double gammaRatio; //saturated buoyancy frequency (scalar, rad/s)
    double rhoS0; //water-vapor density at z=0 (scalar, kg/m^3)
    double hV; //scale height for saturated vapor (scalar, m)
    int nSize;

    double *T; //environmental temperature (vector, K)
    double *gammaEnv; //environmental lapse rate (vector, K/m)
    double *zBar; //saturated lapse rate (vector, K/m)
    double *gammaSat; //weight mean of gammaSat/gammaEnv

    void base_state(); //container function for base state calculations
    FarField(double nS, double T0); //constructor
    ~FarField();


private:

    const double cPD = 1005.7; //heat capacity dry air, constant pressure & 0 C (J/kg-K)
    const double cPV = 1952; //heat capacity water vapor, constant pressure & 0 C (J/kg-K)
    const double cC = 4218; //heat capacity of condensed phases (J/kg-K) (assuming liquid water at 0 C, as is common practice)
    const double RD = 287.0; //specific gas constant dry air (J/kg-K)
    const double p0 = 101325; //standard sea-level pressure (Pa)
    const double epsilon = 0.622; //molecular mass of water relative to air (kg/kg)
    const double dZ = 100;  //height increment (m) (recommended: 100 m)
    const double zMax = 10e3; //maximum height (m)
    const double g = 9.81;
    double L = 2.501e6; //latent heat of vaporization (water -> vapor) (J/kg)
    double eS;
    int n = round(zMax/dZ) + 1;

    std::vector<double> p; //environmental density profile (vector)
    std::vector<double> rhoS; //saturated density vapor profile
    std::vector<double> rS; //saturated mixing ratio (vapor mass/dry air mass)


    double fDiff(double gammaEnvEst, int indx); //define nested function for numerically solving for gammaEnv
    void saturated_vapor_pressure(double Temp);
    void find_gamma_env(); //calculate the atmospheric profile for the environmental lapse rate
    void calc_gamma_ratio(); //calculate the ratio of environmental to saturate lapse rate
    void calc_rhoS0_hV();
    void calc_rho0_hRho();

    /* Necessary functions for root finding */

    double x1, x2, root, guess;
    void bracket(double l, double u, int indx); //bracket the root of a function given an initial guess
    void decker_brent(double l, double u, double tol, int indx); //find the root of a function using the Decker-Brent method
    double sign(double a, double b); //the magnitude of a times the sign of b
    double brack_and_find(double guess, int indx); //container function for root-finding

};


#endif // FARFIELD_H

