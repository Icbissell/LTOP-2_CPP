#include "farfield.h"

FarField::FarField(double Ns, double T0)
{
    /* Make user input variables globally usable */
    this->Ns = Ns;
    this->T0 = T0;

    /* Allocate memory for arrays*/
    gammaEnv = new double[n];
    gammaSat = new double[n];
    zBar = new double[n];
    T = new double[n];

}

FarField::~FarField(){

    /* Deallocate memory*/
    delete[] gammaEnv;
    delete[] gammaSat;
    delete[] zBar;
    delete[] T;
}

double FarField::fDiff(double gammaEnvEst, int indx){

    // Vertical derivative of saturation mixing ratio (DK82, eq. 12)
    double dRS_dZ = rS[indx]*(1+rS[indx]/epsilon)/(RD*T[indx])*(-epsilon*L*gammaEnvEst/T[indx] + g*(1+rS[indx])/(1+rS[indx]/epsilon));

    // Saturated buoyancy frequency using DK82, eq 5, with rT = rS
    double NS2Est = (g/T[indx])*(gammaSat[indx] - gammaEnvEst)
            *(1 + L*rS[indx]/(RD*T[indx])) - g/(1 + rS[indx])*dRS_dZ;

    double diff = pow(Ns,2) - NS2Est;

    return diff;
}

void FarField::saturated_vapor_pressure(double Temp){

    /* Saturated vapor pressure over water and ice, from Goff and Gratch (1946), as modified by Goff (1965).
       See review by Murphy and Koop (2005) for details. The Wegener-Bergeron-Findeisen (WBF) zone, defined
       by the temperature range 268 to 248 K where water and ice coexist. For this zone, eS is set using
       weighted averages for water and ice components. */

    double LV = 2.501e6; // latent heat of vaporization (water -> vapor) (J/kg)
    double LS = 2.834e6; // latent heat of sublimation (ice -> vapor) (J/kg)

    /* Saturated vapor pressure (Pa) in equilibrium with water */
    double eS_Water = pow(10,-7.90298*(373.16/Temp - 1.) + 5.02808 * log10(373.16/Temp)
                          - 1.3816e-7*(pow(10,(11.344*(1. - Temp/373.16))) - 1)
                          + 8.1328e-3*(pow(10,(3.49149*(1. - 373.16/Temp))) - 1) + log10(101324.6));

    /* Saturation vapor pressure (Pa) in equilibrium with ice */
    double eS_Ice = pow(10, -9.09718*((273.16/Temp)-1)- 3.56654*log10(273.16/Temp)
                        + 0.876793*(1 - (Temp/273.16)) + log10(610.71));


    /* Mix values using a factor that goes from zero for T < 248 K
       to one for T > 268 K, corresponding to the WBF zone. */
    double factor = (Temp-248)/(268-248);
    if(factor < 0){factor = 0;} //all ice
    if(factor > 1){factor = 1;} //all water
    eS = eS_Ice*(1-factor) + eS_Water*factor;
    L = LV*factor + LS*(1.-factor);

}


void FarField::find_gamma_env(){

    /* Allocate memry for arrays */
    p.reserve(n);
    rS.reserve(n);
    std::vector<double> rho(n);

    for(int i = 0; i < n; i++){
        zBar[i] = i*dZ;
    }

    /* set base values for temperature and pressure profiles  */
    T[0] = T0;
    p[0] = p0;

    for(int i = 0; i < n; i++){

        /* partial pressure and mixing ratio for saturated water vapor */
        saturated_vapor_pressure(T[i]);
        rS[i] = epsilon*eS/(p[i] - eS);

        /* saturated adiabatic lapse rate, DK82, eq. 19, with 2`Q */
        gammaSat[i] = g*(1 + rS[i])*(1 + L*rS[i]/(RD*T[i]))
                /(cPD + cPV*rS[i] + pow(L,2)*rS[i]*(epsilon + rS[i])/(RD*pow(T[i],2)));

        double gammaEnvGuess = gammaSat[i] - pow(Ns,2)*T[i]/g;

        /* implicitly solve for gammaEnv by finding root of fDiff function*/
        gammaEnv[i] = brack_and_find(gammaEnvGuess, i);

        /* total density, eq. 3.15 in Wallace and Hobbs, 2006 */
        rho[i] = p[i]*(1 + rS[i])/(RD*T[i]*(1. + rS[i]/epsilon));

        /* forward difference for temperature, pressure, and saturated adiabat at next elevation */
        if(i < n-1){
            T[i+1] = T[i] - gammaEnv[i]*dZ;
            double dLnP_dZ = -g*rho[i]/p[i];
            p[i+1] = p[i]*exp(dLnP_dZ*dZ);
        }
    }

    /* check gammaEnv */
    for(int i = 0; i < n; i++){
        if(gammaEnv[i] < 1e-3){
            std::cerr << "Not allowed: Selected NS and T0 for the base state, local regions where gammaEnv < 1 C/km." << std::endl;
            break;
        }
    }

}

void FarField::calc_gamma_ratio(){

    std::vector<double> weights(n);
    rhoS.reserve(n);

    /* calculate water-vapor density */
    for(int i = 0; i < n; i++){
        rhoS[i] = p[i]*rS[i]/(RD*T[i]*(1+rS[i]/epsilon));
    }


    /* Estimate gammaRatio, which is the mean of the
       ratio of the saturated lapse rate over the environmental lapse rate.
       The mean is weighted proportional to the saturated-vapor density. */

    gammaRatio = 0;

    double rhoSSum = 0;
    for(int i = 0; i < n; i++){
        rhoSSum += rhoS[i];
    }

    for(int i = 0; i < n; i++){
        weights[i] = rhoS[i]/rhoSSum;
        gammaRatio += (weights[i]*gammaSat[i]/gammaEnv[i]);
    }

}

void FarField::calc_rhoS0_hV(){

    using namespace Eigen; //linear operations in this section are done using the open-source Eigen library
    /* set Eigen objects for calculations */
    MatrixXf A(n,2);
    VectorXf C(n);
    std::vector<double> weights(n);
    double rhoSSquareSum = 0;

    /* Estimate rhoS0 and hV for vertical water-vapor density.
       These parameters are estimated using an exponential fit, with
       weighting set proportional to the saturated-vapor density, and
       to account for the log transformation used for the fit. */

    for(int i = 0; i < n; i++){
        rhoSSquareSum += pow(rhoS[i],2);
    }

    for(int i = 0; i < n; i++){
        weights[i] = pow(rhoS[i],2)/rhoSSquareSum;
    }

    for(int i = 0; i < n; i++){
        A(i,0) = weights[i];
        A(i,1) = weights[i]*zBar[i];
        C(i) = weights[i]*log(rhoS[i]);
    }

    /* Find the least squares solution of AB = C. The operation done here
       (for efficency) is equivalent to solving the normal equation A'AB = A'C. */
    MatrixXf B = (A.transpose() * A).ldlt().solve(A.transpose() * C);
    rhoS0 = exp(B(0));
    hV = -1/B(1);

}

void FarField::calc_rho0_hRho(){

    using namespace Eigen;

    /* set Eigen objects for calculations */
    MatrixXf A(n,2);
    VectorXf C(n);
    std::vector<double> rhoTotal(n);
    std::vector<double> weights(n);
    double rhoTotalSum = 0;

    /* Exponential fit to estimate rho0 and hRho for total density
       Total density is based on eq. 3.15 in Wallace and Hobbs, 2006.
       Weights account for log transformation used for the fit. */
    for(int i = 0; i < n; i++){
        rhoTotal[i] = p[i]*(1+rS[i])/(RD*T[i]*(1+rS[i]/epsilon));
    }

    for(int i = 0; i < n; i++){
        rhoTotalSum += rhoTotal[i];
    }

    for(int i = 0; i < n; i++){
        weights[i] = rhoTotal[i]/rhoTotalSum;
    }

    for(int i = 0; i < n; i++){
        A(i,0) = weights[i];
        A(i,1) = weights[i]*zBar[i];
        C(i) = weights[i]*log(rhoTotal[i]);
    }

    MatrixXf B = (A.transpose() * A).ldlt().solve(A.transpose() * C);

    rho0 = exp(B(0));
    hRho = -1/B(1);

}

void FarField::base_state(){
    find_gamma_env();
    calc_gamma_ratio();
    calc_rhoS0_hV();
    calc_rho0_hRho();
    nSize = n;
}


/* Root finding functions */
void FarField::bracket(double l, double u, int indx){
    const int pers=1000; //define number of tries to bracket root before giving up
    const double FACTOR=1.6; //size of step when extending bracket

    double f1 = fDiff(l, indx);
    double f2 = fDiff(u, indx);
    bool found = false;

    if (l == u){
        std::cerr << "Bad initial range in zbrac" << std::endl;
    };

    for (int j=0; j<pers;j++) {
        /* check if f1*f2 < 0.0, implying there is a sign change within bracket and
           that the root is bracketed (hello intermediate value theorem) */
        if (f1*f2 < 0.0){
            found = true;
            break;
        }

        /* extend bracket if sign change is not yet within bracket */
        if (std::abs(f1) < std::abs(f2)){
            f1=fDiff(l += FACTOR*(l-u), indx);
        }
        else{
            f2=fDiff(u += FACTOR*(u-l), indx);
        }
    }

    /* update values to match found bracket */
    if(found){
        x1 = l;
        x2 = u;
    }

    else{
        std::cerr << "Range became too big, try again with closer guess: " << std::endl;
        exit(0);
    }
}

double FarField::sign(double a, double b){

    if(a < 0){
        a *= -1;
    }

    if(b < 0){
        a *= -1;
    }

    return a;
}

void FarField::decker_brent(double l, double u, double tol, int indx){

    const int ITMAX=200;
    const double EPS=1.0e-16;
    double a=l,b=u,c=u,d,e,fa=fDiff(a, indx),fb=fDiff(b, indx),fc,p,q,r,s,tol1,xm;
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
        std::cerr << ("Root must be bracketed") << std::endl;
    }
    fc=fb;
    for (int iter=0;iter<ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a;
            fc=fa;
            e=d=b-a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }

        tol1=2.0*EPS*std::abs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (std::abs(xm) <= tol1 || fb == 0.0) {
            root = b;
        }
        if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) {
                q = -q;
            }
            p=std::abs(p);
            double min1=3.0*xm*q-std::abs(tol1*q);
            double min2=std::abs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d;
                d=p/q;
            } else {
                d=xm;
                e=d;
            }
        } else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if (std::abs(d) > tol1){
            b += d;}
        else{
            b += sign(tol1,xm);}
        fb=fDiff(b, indx);
    }
}

double FarField::brack_and_find(double guess, int indx){
    double l = guess;
    double u = guess + 1;
    bracket(l, u, indx);
    decker_brent(x1, x2, 1e-7, indx);
    return root;
}
