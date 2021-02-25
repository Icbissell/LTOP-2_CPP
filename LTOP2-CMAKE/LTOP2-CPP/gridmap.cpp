#include "gridmap.h"

GridMap::GridMap(double azimuth)
{

    get_input();

    n0 = topoVec.size(); // Alter n0 & n1 sizes to match array
    n1 = topoVec[0].size();

    /* Grid-spacing dX and dY is set here to 1000m by default, but can
       be adjusted as desired */

    dX = 1000;
    dY = 1000;

    nXPad = 2*n0;
    nYPad = 2*n1;
    rowLength = n0*n1;
    bigRowLength = nXPad * nYPad;

    hHat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    zInterpolate = new double* [nXPad];
    for(int i = 0; i<nXPad; i++){
        zInterpolate[i] = new double[nYPad];
    }

    pWind = new double* [nXPad];
    for(int i = 0; i < nXPad; i++){
        pWind[i] = new double[nYPad];
    }

    pGrid = new double*[n0];
    for(int i = 0; i<n0; i++){
        pGrid[i] = new double[n1];
    }

    grid = new double* [n0];
    for(int i = 0; i<n0; i++){
        grid[i] = new double[n1];
    }

    this-> azimuth = azimuth;

}


GridMap::~GridMap() {

    delete[] kS;
    delete[] kT;

    delete[] Sxy;
    delete[] Txy;
    delete[] Xst;
    delete[] Yst;

    delete[] T;
    delete[] zBar;
    delete[] gammaEnvArr;
    delete[] gammaSatArr;

    for(int i =0; i < n0; i++){
        delete[] grid[i];
    }

    delete[] grid;

    for(int i =0; i < n0; i++){
        delete[] pGrid[i];
    }

    delete[] pGrid;

    for(int i =0; i < nXPad; i++){
        delete[] zInterpolate[i];
    }

    delete[] zInterpolate;

    for(int i =0; i < nXPad; i++){
        delete[] pWind[i];
    }

    delete[] pWind;
}


void GridMap::get_input(){

    /* Fill in topography grid */
    std::string topoPath = "LTOP2-CPP/himilaya.csv"; // link topography DEM here (square)
    std::ifstream inGrid(topoPath);
    std::string line, field;

    while (getline(inGrid,line)){ // get next line in file
        std::vector<std::string> v;
        std::stringstream ss(line);

        while (getline(ss,field,','))  // break line into comma delimitted fields
        {
            v.push_back(field);  // add each field to the 1D array
        }

        std::vector<double> doubleV;
        int vSize = v.size();
        for(int i = 0; i < vSize; i ++){
            doubleV.push_back(std::stod(v[i]));
        }
        topoVec.push_back(doubleV);  // add vector to 2D array
    }

    double iN = topoVec.size();
    double jN = topoVec[0].size();

    for(int i = 0; i < iN; i++){
        for(int j=0; j < jN; j++){
            if(std::isnan(topoVec[i][j])){
                std::cerr << "Topography file contains nans, please input real values \n";
                throw;
            }
        }
    }

    /* Remove topography below sea-level*/
    for(int i = 0; i < iN; i++){
        for(int j=0; j < jN; j++){
            if(topoVec[i][j] < 0){
                topoVec[i][j] = 0;
            }
        }
    }

}


void GridMap::fill_base_state(double Ns, double T0){

    FarField FarField(Ns, T0);
    FarField.base_state();

    gammaRatio = FarField.gammaRatio;
    rhoS0 = FarField.rhoS0;
    hV = FarField.hV;
    rho0 = FarField.rho0;
    hRho = FarField.hRho;

    n = FarField.nSize;
    zBar = new double[n];
    gammaEnvArr = new double[n];
    gammaSatArr = new double[n];
    T = new double[n];

    for(int i = 0; i < n; i++){
        zBar[i] = FarField.zBar[i];
        gammaSatArr[i] = FarField.gammaSat[i];
        gammaEnvArr[i] = FarField.gammaEnv[i];
        T[i] = FarField.T[i];
    }
}

void GridMap::calc_tauF(){

    LinearInterp zBarInt(T, zBar, n);
    double zBar258 = zBarInt.interp(258);

    /* Base-state environmental lapse rate at zBar258 */

    LinearInterp gammaEnvInt(zBar, gammaEnvArr, n);
    double gammaEnvBar258 = gammaEnvInt.interp(zBar258);

    /* Base-state Saturated lapse rate at zBar258 */

    LinearInterp gammaSatInt(zBar, gammaSatArr, n);
    double gammaSatBar258 = gammaSatInt.interp(zBar258);

    std::complex<double> *zHat = new std::complex<double>[padstRowLength];

    for (int i = 0; i < padstRowLength; i++){
        zHat[i] = std::complex<double>(hHat[i][REAL], hHat[i][IMAG])*exp((img*kZ[i]+comp_one/(2*hRho))*zBar258);
    }

    zHat[0] = 0;

    fftw_complex *zPrimeHat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);
    fftw_complex *zPrime = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    for (int i = 0; i < padstRowLength; i++){
        zPrimeHat[i][REAL] = real(zHat[i]);
        zPrimeHat[i][IMAG] = imag(zHat[i]);
    }

    delete[] zHat;

    ifft(nSPad, nTPad, zPrimeHat, zPrime);

    fftw_free(zPrimeHat);

    std::vector<double> zPrimeVec;

    /* Remove 0 padding from arrays */
    for (int i = 0; i < nSPad; i++){
        for(int j = 0; j < nTPad; j++){
            int indx = nTPad*i + j;
            if(i < sLength && j < tLength){
                zPrimeVec.push_back(zPrime[indx][REAL]);
            }

        }
    }

    fftw_free(zPrime);
    double *z258 = new double[stRowLength];

    for (int i = 0; i < stRowLength; i++){
        z258[i] = zBar258 + zPrimeVec[i]*(2.-gammaSatBar258/gammaEnvBar258);
    }

    /* Calculate mean fall time. Velocities wR and wS (m/s) for rain and snow, respectively.
       References: Langleben, 1954; White et al., 2002; Yuter and Houze, 2003;
       Barstad and Schuller, 2011. */

    double wFSnow = -1.;
    double wFRain = -6.;

    std::vector<double> tauFVec(stRowLength);

    /* Find tauF over grid and calculate mean fall time*/

    for(int i =0; i < sLength; i++){
        for(int j =0; j < tLength; j++){
            int indx = i*tLength + j;
            int lessEq;
            int more;

            (z258[indx] <= zInterpolate[i][j]) ? (lessEq = 1) : (lessEq = 0);
            (z258[indx] > zInterpolate[i][j]) ? (more = 1) : (more = 0);

            tauFVec[indx] = (lessEq*(-hV/wFSnow)+(more)*-((z258[indx]-zInterpolate[i][j])/wFRain
                                                          + hV*exp(-(z258[indx] - zInterpolate[i][j])/hV)/wFSnow));
        }
    }

    tauF = mean1(tauFVec, stRowLength);

    delete[] z258;
}

void GridMap::fill_grid_topo(){

    /* Populate topography with values from the DEM file*/

    for(int i =0; i < n0; i++){
        for(int j =0; j < n1; j++){
            grid[i][j] = topoVec[i][j];

        }
    }

}

void GridMap::interp_axes(int xsize, int ysize, double angle){

    /* Allocate arrays */
    double *xCoord = new double[n1];
    double *yCoord = new double[n0];
    Sxy = new double[rowLength];
    Txy = new double[rowLength];

    /* Calculate sin and cos of azimuth*/
    double sinAng = sind(angle), cosAng = cosd(angle);
    double orthSinAng = sind(angle + 90), orthCosAng = cosd(angle + 90);

    /* Set up coordinates in the x-direction in meters*/
    for(int i = 0; i < n1; i++){
        xCoord[i] = (double)i*dX;
    }

    /* Set up coordinates in the y-direction in meters*/
    for(int i = 0; i < n0; i++){
        yCoord[i] = (double)i*dY;
    }

    /* Calculate s,t coordinates for the x,y nodes of geographic grid */
    for(int i = 0; i < xsize; i++){
        for(int j = 0; j < ysize; j++){
            int indx = i*ysize+j;
            Sxy[indx] = xCoord[j]*sinAng + yCoord[i]*cosAng;
            Txy[indx] = -xCoord[j]*cosAng + yCoord[i]*sinAng;
        }
    }

    /* Deallocate arrays*/
    delete[] xCoord;
    delete[] yCoord;

    double sMin = find_max_min(Sxy, "min", rowLength);
    double sMax = find_max_min(Sxy, "max", rowLength);
    double tMin = find_max_min(Txy, "min", rowLength);
    double tMax = find_max_min(Txy, "max", rowLength);

    /* Calculate dS and dT, the node spacing for the new grid
       that matches the nodal spacing of the original grid given the
       new azimuth direction using the ellipse equation */

    dS = sqrt(1./(pow(sinAng/dX,2) + pow(cosAng/dY,2)));
    dT = sqrt(1./(pow(orthSinAng/dX,2) + pow(orthCosAng/dY,2)));

    dS = (sMax - sMin)/ceil((sMax - sMin)/dS);
    dT = (tMax - tMin)/ceil((tMax - tMin)/dT);

    /* Calculate the lengths of the s axis and t axis*/
    sLength = (sMax - sMin)/dS + 1;
    tLength = (tMax - tMin)/dT + 1;

    /* Set size of s and t vectors*/
    s.reserve(sLength);
    t.reserve(tLength);

    /* Populate s and t axes with coordinates*/
    for(int i = 0; i < sLength; i++){
        s[i] = sMin + i*dS;
    }

    for(int i = 0; i < tLength; i++){
        t[i] = tMin + i*dT;
    }

    /* Calculate length of padded arrays*/
    nSPad = 2*sLength;
    nTPad = 2*tLength;

    stRowLength = sLength*tLength;
    padstRowLength = nSPad*nTPad;

    Xst = new double[stRowLength];
    Yst = new double[stRowLength];

    dS = s[1] - s[0];
    dT = t[1] - t[0];

    /* Calculate x,y coordinates for the s,t nodes of geographic grid */
    for(int i = 0; i < sLength; i++){
        for(int j= 0; j < tLength; j++){
            int indx = i*tLength+j;
            Xst[indx] = s[i]*sinAng - t[j]*cosAng;
            Yst[indx] = s[i]*cosAng + t[j]*sinAng;
        }
    }

}

void GridMap::rotate(){

    /* Allocate memory for arrays */
    double *xCoord = new double[n0];
    double *yCoord = new double[n1];

    /* Calculate size, in meters, in x and y directions*/
    int xSize = n0*dX;
    int ySize = n1*dY;

    /* Populate coordinates on x and y axes*/
    for(int i = 0; i < n0; i++){
        xCoord[i] = (double)i*dX;
    }

    for(int i = 0; i < n1; i++){
        yCoord[i] = (double)i*dY;
    }

    /* Initialize the interpolator */
    BilinInterp bilinterp(n0, n1, xCoord, yCoord, grid);

    /* Fill in array 'zInterpolate' with values from 'grid' array. Interpolation is used
       when the queried point is not in the original grid e.g. when the grid has been
       rotated. If the query point is outside the bounds of the original grid, the value
       is filled with 0 */

    for(int i = 0; i < sLength; i++){
        for(int j = 0; j < tLength; j++){

            int indx = i*tLength+j;

            if((Xst[indx] >= 0) && (Xst[indx] <= ySize) && (Yst[indx] >= 0) && (Yst[indx] <= xSize)){
                zInterpolate[i][j] = bilinterp.interp(Yst[indx],Xst[indx]);
            }

            else{
                zInterpolate[i][j] = 0;
            }
        }
    }

    /* Deallocate memory for arrays */
    delete[] xCoord;
    delete[] yCoord;
}

void GridMap::reverse_rotation(double **inGrid, double **outGrid){

    /* Allocate memory for arrays */
    double *sCoord = new double[sLength];
    double *tCoord = new double[tLength];

    /* Populate coordinates on x and y axes*/
    for(int i = 0; i < sLength; i++){
        sCoord[i] = s[i];
    }

    for(int i = 0; i < tLength; i++){
        tCoord[i] = t[i];
    }

    /* Initialize the interpolator */
    BilinInterp bilinterp(sLength, tLength, sCoord, tCoord, inGrid);

    /* Fill in array 'precGrid' with values from 'pWind' array. Interpolation is used
       when the queried point is not in the original grid e.g. when the grid has been
       rotated. */

    for(int i = 0; i < n0; i++){
        for(int j =0; j < n1; j++){
            int indx = i*n1+j;
            outGrid[i][j] = bilinterp.interp(Sxy[indx],Txy[indx]);
        }
    }

    /* Deallocate memory for arrays */
    delete[] sCoord;
    delete[] tCoord;

}

void GridMap::topo_to_wave_domain() {

    /* Allocate memory for the complex array rGridIn*/
    fftw_complex *rGridIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    /* Simultaneously populate rGridIn with values from topography and add 0 padding.
       Padding is twice the length of the rotated array and accounts for the response
       of lifting over the topography by lessening wrap-around in the wave domain */

    for (int i = 0; i < nSPad; i++){
        for(int j = 0; j < nTPad; j++){
            int indx = nTPad*i+j;
            if(i < sLength && j < tLength){
                rGridIn[indx][REAL] = zInterpolate[i][j];
                rGridIn[indx][IMAG] = 0;
            }

            else{
                rGridIn[indx][REAL] = rGridIn[indx][IMAG] = 0;
            }
        }
    }

    /* Transform topography to frequency domain*/
    fft(nSPad, nTPad, rGridIn, hHat);

    /* Deallocate array*/
    fftw_free(rGridIn);

}

void GridMap::fill_waveNumber(){

    /* Allocate memory for arrays*/
    kS = new double[rowLength];
    kT = new double[rowLength];

    /* t wavenumber calculation */

    int i_kTMostNeg = ceil(tLength/2)+1;

    for(int i=0; i < nTPad; i++){
        kT[i] = (double)i/nTPad;
    }

    for(int i=i_kTMostNeg; i < nTPad; i++){
        kT[i] = kT[i] - 1;
    }

    for(int i=0; i < nTPad; i++){
        kT[i] = 2*M_PI*kT[i]/dT;
    }

    dkT = kT[1] - kT[0];

    /* s wavenumber calculation */

    int i_kSMostNeg = ceil(nSPad/2)+1;

    for(int i=0; i < nSPad; i++){
        kS[i] = (double)i/nSPad;
    }

    for(int i=i_kSMostNeg; i < nSPad; i++){
        kS[i] = kS[i] - 1;
    }

    for(int i = 0; i < nSPad; i++){
        kS[i] = 2*M_PI*kS[i]/dS;
    }

    dkS = kS[2] - kS[1];

}

void GridMap::treat_singularities(){

    /* Set size of denominator vector */
    denominator.reserve(padstRowLength);

    /* Fill in values of denominator according to (U*kS)^2-fC^2 for the calculation of
       m or kZ, the vertical wave number */
    for(int i = 0; i < nSPad; i++){
        for(int j = 0; j < nTPad; j++){
            int indx = i*nTPad + j;
            denominator[indx] = pow((U*kS[i]),2.) - pow(fC,2.);
        }
    }

    /* Treat singularities caused when abs(U*kS)==abs(fC) or denominator = 0.
       The justification is that at this singularity, the vertical velocities go to zero.
       The modification to hHat removes the excitation at the singularity, and
       the modification to the demoninator avoids errors due to irrelevant nans,
       from Queney, 1947, p. 46-48 */

    for(int i = 0; i < padstRowLength; i++){
        if(abs(denominator[i]) < pow(U*dkS/2.,2)){
            hHat[i][REAL] = hHat[i][IMAG] = 0;
            denominator[i] = EPS;
        }
    }

}

void GridMap::calc_kZ(){

    /* Declare vector with filled squares of vertical wave numbers */
    std::vector<std::complex <double>> kZ2 (padstRowLength);

    /* Set size of kZ vector */
    kZ.reserve(padstRowLength);

    /* Calculate kZ, the vertical wavenumber. If kZ2>0, then sqrt(kZ2) is real (propagating wave).
       Simultaneously set values of kZ to match the sign of its associated kS value.
       This ensures that progagating waves have an upwind phase tilt. */

    for (int i = 0; i < nSPad; i++){
        for (int j = 0; j < nTPad; j++){
            int indx = i*nTPad+j;
            kZ2[indx] = std::complex<double>(pow(kS[i],2.) + pow(kT[j],2.),0)*((pow(Ns,2.) - pow(U*kS[i],2.))/denominator[indx]) - 1./(4.*pow(hRho,2.)); //check in on complex double part
            kZ[indx] = sqrt(kZ2[indx]);

            if((real(kZ2[indx]) > 0) && (kS[i] < 0)){
                kZ[indx] = -kZ[indx];
            }
        }
    }
}

void GridMap::calc_response_functions(){

    /* Set size of vectors */
    hHatPair.reserve(padstRowLength);
    GVHat.reserve(padstRowLength);
    GCHat.reserve(padstRowLength);
    GFHat.reserve(padstRowLength);

    /* Convert hHat from fftw_complex to vector for calculations */
    for (int i = 0; i < padstRowLength; i++){
        hHatPair[i] = std::complex<double>(hHat[i][REAL], hHat[i][IMAG]);
    }

    fftw_free(hHat);

    /* Populate values of three response functions: GVHat, GCHat, GFHat */
    for (int i = 0; i < nSPad; i++){
        for (int j = 0; j < nTPad; j++){
            int indx = i*nTPad + j;
            GVHat[indx] = gammaRatio*rhoS0*img*U*kS[i]/(comp_one-img*kZ[indx]*hV);
            GCHat[indx] = (pow(tauC*(kappa*std::complex<double>(pow(kS[i],2) + pow(kT[j],2),0) - img*U*kS[i]) + comp_one, -1));
            GFHat[indx] = (pow(tauF*(kappa*std::complex<double>(pow(kS[i],2) + pow(kT[j],2),0) - img*U*kS[i])+ comp_one, -1));
        }
    }

}

void GridMap::convolution_functions(){

    /* Allocate memory for arrays */
    std::complex<double> *pStarHatPair = new std::complex<double>[padstRowLength];
    std::complex<double> *QCStarPair = new std::complex<double>[padstRowLength];
    std::complex<double> *QPStarPair = new std::complex<double>[padstRowLength];

    fftw_complex *QCStarWindHat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);
    fftw_complex *QFStarWindHat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);
    fftw_complex *pStarHat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);
    fftw_complex *pStarWind = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);
    fftw_complex *QCStarWind = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);
    fftw_complex *QFStarWind = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * padstRowLength);

    /* Calculate PStarHat, reference precipitation rate for wave domain.
       Calculate vapor ratio and moisture-corrected precipitation rate */
    for (int i = 0; i < padstRowLength; i++){
        pStarHatPair[i] = GVHat[i]*GCHat[i]*GFHat[i]*hHatPair[i];
        pStarHat[i][REAL] = real(pStarHatPair[i]);
        pStarHat[i][IMAG] = imag(pStarHatPair[i]);

        QCStarPair[i] = tauC*GVHat[i]*GCHat[i]*hHatPair[i];
        QCStarWindHat[i][REAL] = real(QCStarPair[i]);
        QCStarWindHat[i][IMAG] = imag(QCStarPair[i]);

        QPStarPair[i] = tauF*GVHat[i]*GCHat[i]*GFHat[i]*hHatPair[i];
        QFStarWindHat[i][REAL] = real(QPStarPair[i]);
        QFStarWindHat[i][IMAG] = imag(QPStarPair[i]);
    }

    /* Deallocate memory */
    delete[] pStarHatPair;
    delete[] QPStarPair;
    delete[] QCStarPair;

    /* Perform inverse Fourier transform to return to space domain */
    ifft(nSPad, nTPad, pStarHat, pStarWind);
    ifft(nSPad, nTPad, QCStarWindHat, QCStarWind);
    ifft(nSPad, nTPad, QFStarWindHat, QFStarWind);

    /* Deallocate memory */
    fftw_free(pStarHat);
    fftw_free(QFStarWindHat);
    fftw_free(QCStarWindHat);

    /* Remove 0 padding from arrays */
    for (int i = 0; i < nSPad; i++){
        for(int j = 0; j < nTPad; j++){
            int indx = nTPad*i + j;
            if(i < sLength && j < tLength){
                pStarWindVec.push_back(pStarWind[indx][REAL]);
                QFStarWindVec.push_back(QFStarWind[indx][REAL]);
                QCStarWindVec.push_back(QCStarWind[indx][REAL]);
            }

        }
    }

    for (int i = 0; i < stRowLength; i++){

        if(pStarWindVec[i] < 0){
            pStarWindVec[i] = 0;
        }

        if(QCStarWindVec[i] < 0){
            QCStarWindVec[i] = 0;
        }

        if(QFStarWindVec[i] < 0){
            QFStarWindVec[i] = 0;
        }
    }

    /* Deallocate memory */
    fftw_free(QFStarWind);
    fftw_free(pStarWind);
    fftw_free(QCStarWind);

}


void GridMap::final_prec_calc(){

    /* Allocate memory for arrays */
    double **cumTrapzArrIn = new double* [sLength];
    for(int i = 0; i < sLength; i++){
        cumTrapzArrIn[i] = new double[tLength];
    }

    double **cumTrapzArrOut = new double* [sLength];
    for(int i = 0; i < sLength; i++){
        cumTrapzArrOut[i] = new double[tLength];
    }

    std::vector<double> QTStarWind(stRowLength);

    for (int i = 0; i < stRowLength; i++){
        QTStarWind[i] = rhoS0*hV+QCStarWindVec[i]+QFStarWindVec[i];
    }

    for (int i = 0; i < sLength; i++){
        for (int j = 0; j < tLength; j++){
            int indx = i*tLength+j;
            cumTrapzArrIn[i][j] = (pStarWindVec[indx]*dS)/(QTStarWind[indx]);
        }
    }

    cum_trapz_integrate(cumTrapzArrIn, cumTrapzArrOut, sLength, tLength);

    std::vector<double> fVWind(stRowLength);

    for (int i = 0; i < sLength; i++){
        for (int j = 0; j < tLength; j++){
            int indx = i*tLength + j;
            fVWind[indx] = (rhoS0*hV/QTStarWind[indx])*exp(-(fP/U)*cumTrapzArrOut[i][j]);
        }
    }

    /* Calculate pWind, moisture-corrected precipitation-rate field */
    for (int i = 0; i < sLength; i++){
        for(int j = 0; j < tLength; j++){
            int indx = i*tLength + j;
            pWind[i][j] = fVWind[indx]*pStarWindVec[indx];
        }
    }

    /* Deallocate memory */
    for(int i =0; i < sLength; i++){
        delete[] cumTrapzArrIn[i];
    }

    delete[] cumTrapzArrIn;

    for(int i =0; i < sLength; i++){
        delete[] cumTrapzArrOut[i];
    }

    delete[] cumTrapzArrOut;

}


void GridMap::calc_most(){

    /* Calculate the maximum rate of precipitation */
    mostPrec = find_max_min_2d(pGrid, "max", n0, n1);

    /* Adjust to meters/year and millimeters/hour units */
    mostPrec_ma = mostPrec*3600*24*365.25/1e3;
    mostPrec_mmhr = mostPrec*3.6e3;
    max_el = find_max_min_2d(grid, "max", n0, n1);

}

void GridMap::LTOP_calc(){

    fill_base_state(Ns, T0);
    fill_grid_topo();
    interp_axes(n0,n1,azimuth);
    fill_waveNumber();
    rotate();
    topo_to_wave_domain();
    treat_singularities();
    calc_kZ();
    calc_tauF();
    calc_response_functions();
    convolution_functions();
    final_prec_calc();
    reverse_rotation(pWind, pGrid);
    calc_most();

}


void GridMap::make_plots(){

    /* Make topography Plot using gnuplot pipeline */

    /* Convert to vector for easy plotting */

    std::vector<std::vector<double>> gridVec;

    for(int i = 0; i < n0; i++){
        std::vector<double> gridVecHelper;
        for(int j = 0; j < n1; j++){
            gridVecHelper.push_back(grid[i][j]);
        }
        gridVec.push_back(gridVecHelper);
    }

    /* Plot using Gnuplot pipeline */

    Gnuplot gp1;
    gp1 << "set view map \n";
    gp1 << "set title 'Topography' \n";
    gp1 << "load 'LTOP2-CPP/parula.pal' \n";
    gp1 << "splot" << gp1.file1d(gridVec) << "matrix with image \n";
    gp1 << "set term postscript \n";
    gp1 << "set output 'Topography.ps' \n";
    gp1 << "replot" << std::endl;


    /* Make precipitation plot */

    /* Convert to vector for easy plotting */

    std::vector<std::vector<double>> precVec;

    for(int i = 0; i < n0; i++){
        std::vector<double> precVecHelper;
        for(int j = 0; j < n1; j++){
            precVecHelper.push_back(pGrid[i][j]*(double)3600*24*365.25/1e3);
        }
        precVec.push_back(precVecHelper);
    }

    Gnuplot gp2;
    gp2 << "set view map \n";
    gp2 << "set title 'Precipitation with Moisture Balance(m/a)' \n";
    gp2 << "load 'LTOP2-CPP/parula.pal' \n";
    gp2 << "splot" << gp2.file1d(precVec) << "matrix with image \n";
    gp2 << "set term postscript \n";
    gp2 << "set output 'Precipitation_Grid.ps' \n";
    gp2 << "replot" << std::endl;

#ifdef _WIN32
    // For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
    // the gnuplot window doesn't get closed.
    std::cout << "Press enter to exit." << std::endl;
    std::cin.get();
#endif


}
