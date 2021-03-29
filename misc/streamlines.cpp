void streamline()
{

    double zL0Map = 500;
    int zS0Prev = 1;
    //if (isFirst==true || zL0Map!=zS0Prev){
    zS0Prev = zL0Map;

    std::complex<double> *zPrimeHat = new std::complex<double>[padstRowLength];
    for (int i = 0; i < padstRowLength; i++)
    {
        zPrimeHat[i] = std::complex<double>(hHat[i][REAL], hHat[i][IMAG]) * std::exp((img * kZ[i] + comp_one / (2 * hRho)) * zL0Map);
    }
    std::complex<double> *sPrimeHat = new std::complex<double>[padstRowLength];
    for (int i = 0; i < nSPad; i++)
    {
        for (int j = 0; j < nTPad; j++)
        {
            sPrimeHat[i * nTPad + j] = ((img * kS[i] + fC * kT[j]) / (U * kS[i])) / (pow(kS[i], 2.) + pow(kT[j], 2.)) * (img * kZ[i * nTPad + j] + 1. / (2. * hRho)) * zPrimeHat[i * nTPad + j];
        }
    }

    sPrimeHat[0] = 0;

    std::complex<double> *tPrimeHat = new std::complex<double>[padstRowLength];
    for (int i = 0; i < nSPad; i++)
    {
        for (int j = 0; j < nTPad; j++)
        {
            tPrimeHat[i * nTPad + j] = ((img * kS[i] - fC * kT[j]) / (U * kS[i])) / (pow(kS[i], 2.) + pow(kT[j], 2.)) * (img * kZ[i * nTPad + j] + 1. / (2. * hRho)) * zPrimeHat[i * nTPad + j];
        }
    }

    tPrimeHat[0] = 0;

    fftw_complex *zSWindPrime = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);
    fftw_complex *zSWind = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    for (int i = 0; i < padstRowLength; i++)
    {
        zSWindPrime[i][REAL] = real(zPrimeHat[i]);
        zSWindPrime[i][IMAG] = imag(zPrimeHat[i]);
    }
    delete[] zPrimeHat;

    ifft(nSPad, nTPad, zSWindPrime, zSWind);
    fftw_free(zSWindPrime);

    for (int i = 0; i < padstRowLength; i++)
    {
        zSWind[i][REAL] += zL0Map;
    }

    double **zSwindArr = new double *[n0];
    for (int i = 0; i < n0; i++)
    {
        zSwindArr[i] = new double[n1];
    }

    for (int i = 0; i < n0; i++)
    {
        for (int j = 0; j < n1; j++)
        {
            zSwindArr[i][j] = zSWind[i * n1 + j][REAL];
        }
    }

    fftw_free(zSWind);

    /* s direction calculation */

    fftw_complex *sSWindPrime = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);
    fftw_complex *sSWind = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    for (int i = 0; i < padstRowLength; i++)
    {
        sSWindPrime[i][REAL] = real(sPrimeHat[i]);
        sSWindPrime[i][IMAG] = imag(sPrimeHat[i]);
    }
    delete[] sPrimeHat;

    ifft(nSPad, nTPad, sSWindPrime, sSWind);
    fftw_free(sSWindPrime);

    for (int i = 0; i < nSPad; i++)
    {
        for (int j = 0; j < nTPad; j++)
        {
            sSWind[i * nTPad + j][REAL] += s[i];
        }
    }

    double **sSwindArr = new double *[n0];
    for (int i = 0; i < n0; i++)
    {
        sSwindArr[i] = new double[n1];
    }

    for (int i = 0; i < n0; i++)
    {
        for (int j = 0; j < n1; j++)
        {
            sSwindArr[i][j] = sSWind[i * n1 + j][REAL];
        }
    }

    fftw_free(sSWind);

    /* t-direction calculation */

    fftw_complex *tSWindPrime = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);
    fftw_complex *tSWind = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 4 * bigRowLength);

    for (int i = 0; i < padstRowLength; i++)
    {
        tSWindPrime[i][REAL] = real(tPrimeHat[i]);
        tSWindPrime[i][IMAG] = imag(tPrimeHat[i]);
    }
    delete[] tPrimeHat;

    ifft(nSPad, nTPad, tSWindPrime, tSWind);
    fftw_free(tSWindPrime);

    for (int i = 0; i < nSPad; i++)
    {
        for (int j = 0; j < nTPad; j++)
        {
            tSWind[i * nTPad + j][REAL] += t[j];
        }
    }

    double **tSwindArr = new double *[n0];
    for (int i = 0; i < n0; i++)
    {
        tSwindArr[i] = new double[n1];
    }

    for (int i = 0; i < n0; i++)
    {
        for (int j = 0; j < n1; j++)
        {
            tSwindArr[i][j] = tSWind[i * n1 + j][REAL];
        }
    }

    fftw_free(tSWind);

    //}

    std::vector<double> sBarPath(rowLength);
    std::vector<double> tBarPath(rowLength);
    std::vector<double> xBarPath(n0);
    std::vector<double> yBarPath(n1);

    for (int i = 0; i < n0; i++)
    {
        for (int j = 0; j < n1; j++)
        {
            sBarPath[i * n1 + j] = xBarPath[i] * sind(azimuth) + yBarPath[j] * cosd(azimuth);
            tBarPath[i * n1 + j] = -xBarPath[i] * cosd(azimuth) + yBarPath[j] * sind(azimuth);
        }
    }

    std::vector<double> sL(rowLength);
    std::vector<double> tL(rowLength);
    std::vector<double> xL(rowLength);
    std::vector<double> yL(rowLength);

    double *xCoord = new double[n1];
    double *yCoord = new double[n0];

    /* Set up coordinates in the x-direction in meters*/
    for (int i = 0; i < n1; i++)
    {
        xCoord[i] = (double)i * dX;
    }

    /* Set up coordinates in the y-direction in meters*/
    for (int i = 0; i < n0; i++)
    {
        yCoord[i] = (double)i * dY;
    }

    BilinInterp sTerp(n0, n1, xCoord, yCoord, sSwindArr);
    BilinInterp tTerp(n0, n1, xCoord, yCoord, sSwindArr);

    for (int i = 0; i < rowLength; i++)
    {
        sL[i] = sTerp.interp(tBarPath[i], sBarPath[i]);
        tL[i] = tTerp.interp(tBarPath[i], sBarPath[i]);
        xL[i] = sL[i] * sind(azimuth) - tL[i] * cosd(azimuth);
        yL[i] = sL[i] * cosd(azimuth) + tL[i] * sind(azimuth);
    }

    BilinInterp zL(n0, n1, xCoord, yCoord, zSwindArr);

    for (int i = 0; i < n0; i++)
    {
        delete[] tSwindArr[i];
    }

    delete[] tSwindArr;

    for (int i = 0; i < n0; i++)
    {
        delete[] sSwindArr[i];
    }

    delete[] sSwindArr;

    for (int i = 0; i < n0; i++)
    {
        delete[] zSwindArr[i];
    }

    delete[] zSwindArr;

    delete[] xCoord;
    delete[] yCoord;
}

void windPath(double xPath0, double yPath0, double azimuth, int xSize, int ySize, double *x, double *y)
{
    int rowLength = xSize * ySize;
    std::vector<double> xyLimits(4);
    double find_max_min(double *x, std::string choice, int xSize);

    double dX = x[1] - x[0];
    double dY = y[1] - y[0];
    double dSPath = sqrt(1. / (pow((sind(azimuth) / dX), 2.) + pow((cosd(azimuth) / dY), 2.)));

    double xBox[] = {find_max_min(x, 'min', xSize), find_max_min(x, 'max', xSize)};
    double yBox[] = {find_max_min(y, 'min', ySize), find_max_min(y, 'max', ySize)};

    double dMax;
    double *dMaxComp = new double[rowLength];
    for (int i = 0; i < xSize; i++)
        for (int j = 0; j < ySize; j++)
        {
            {
                dMaxComp[i * ySize + j] = (sqrt(pow((x[i] - xPath0), 2) + pow((y[j] - yPath0), 2.)));
            }
        }

    dMax = find_max_min(dMaxComp, 'max', rowLength);

    delete[] dMaxComp;

    double xPath[] = {xPath0 + sind(azimuth) * (-dMax), xPath0 + sind(azimuth) * (dMax)};
    double yPath[] = {yPath0 + cosd(azimuth) * (-dMax), yPath0 + cosd(azimuth) * (dMax)};

    std::vector<std::pair<double, double>> crossVec;
    crossVec = polySect(xBox, yBox, xPath, yPath, 2);

    if (crossVec.size() == 0)
    {
        return;
    }

    if (((xPath[1] - xPath[0]) * sind(azimuth) + (yPath[1] - yPath[0]) * cosd(azimuth)) < 0)
    {
        std::swap(xPath[1], xPath[0]);
        std::swap(yPath[1], yPath[0]);
    }

    std::vector<double> sPath;
    sPath.push_back(sind(azimuth)*xPath[0] + cosd(azimuth)*yPath[0]);
    sPath.push_back(sind(azimuth)*xPath[1] + cosd(azimuth)*yPath[1]);

    int nS = floor((sPath[1] - sPath[0])/dSPath) + 1;

    double xInc = (xPath[1] - xPath[0])/(nS - 1);
    double yInc = (yPath[1] - yPath[0])/(nS - 1);
    double sInc = (sPath[1] - sPath[0])/(nS - 1);

    double* xPathVec = new double[nS];
    double* yPathVec = new double[nS];
    double* sPathVec = new double[nS];

    for (int i = 0; i < nS; i++)
    {
        xPathVec[i] = xPath[0] + xInc*i;
        yPathVec[i] = yPath[0] + yInc*i;
        sPathVec[i] = sPath[0] + sInc*i;
    }

    std::vector<std::pair<double, double> > upCrossVec;
    upCrossVec = polySect(xBox, yBox, xPathVec, yPathVec, 2);
    std::vector<std::pair<double, double> > sLimits;

    if (upCrossVec.size() < 2)
    {
        sLimits.push_back(std::make_pair(sPathVec[0], sPathVec[nS]));
        return;
    }

    std::vector<double> xLimDiff;
    std::vector<double> yLimDiff;

    for (int i = 1; i < upCrossVec.size(); i++)
    {
        xLimDiff.push_back(upCrossVec[i].first - upCrossVec[i-1].first);
        yLimDiff.push_back(upCrossVec[i].second - upCrossVec[i-1].second);
    }

    std::vector<bool> check;


    if(abs(xLimDiff[1])>abs(yLimDiff[1])){
        LinearInterp sLims(xPathVec, sPathVec, nS);
        sLimits.push_back(sLims.interp(upCrossVec[0].first));
        sLimits.push_back(sLims.interp(upCrossVec[1].first));
    }

    else{
        LinearInterp sLims(yPathVec, sPathVec, nS);
        sLimits.push_back(sLims.interp(upCrossVec[0].second));
        sLimits.push_back(sLims.interp(upCrossVec[1].second));
    }

    delete[] xPathVec;
    delete[] yPathVec;
    delete[] sPathVec;
}

bool sign(double a)
{
    if (a < 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

double func(double x, double m, double b)
{
    double y = m * x + b;
    return y;
}

std::vector<std::pair<double, double>> polySect(double *xBox, double *yBox, double *xVec, double *yVec, int vecSize)
{
    std::vector<std::pair<double, double>> crossVec;

    for (int i = 1; i < vecSize; i++)
    {
        double m = (yVec[i] - yVec[i - 1]) / (xVec[i] - xVec[i - 1]);
        double b = yVec[i] - m * xVec[i];

        //check if function crosses lower limit of box

        if ((sign(func(xVec[i], m, b) - yBox[0]) != sign(func(xVec[i - 1], m, b) - yBox[0])) && ((xVec[i] > xBox[0] && xBox[i] < xBox[1]) || (xVec[i - 1] > xBox[0] && xBox[i - 1] < xBox[1])))
        {
            crossVec.push_back(std::make_pair((yBox[0] - b) / m, yBox[0]));
        }

        //check if function crosses upper limit of box

        if ((sign(func(xVec[i], m, b) - yBox[1]) != sign(func(xVec[i - 1], m, b) - yBox[1])) && ((xVec[i] > xBox[0] && xBox[i] < xBox[1]) || (xVec[i - 1] > xBox[0] && xBox[i - 1] < xBox[1])))
        {
            crossVec.push_back(std::make_pair((yBox[1] - b) / m, yBox[1]));
        }

        //check if function crosses leftmost limit of box

        if ((sign(xVec[i] - xBox[0]) != sign(xVec[i - 1] - xBox[0])) && (((func(xVec[i], m, b) > yBox[0]) && (func(xVec[i], m, b) < yBox[1])) || ((func(xVec[i - 1], m, b) > yBox[0]) && (func(xVec[i - 1], m, b) < yBox[1]))))
        {
            crossVec.push_back(std::make_pair(xBox[0], (m * xBox[0] + b)));
        }

        //check if function crosses rightmost limit of box

        if ((sign(xVec[i] - xBox[1]) != sign(xVec[i - 1] - xBox[1])) && (((func(xVec[i], m, b) > yBox[0]) && (func(xVec[i], m, b) < yBox[1])) || ((func(xVec[i - 1], m, b) > yBox[0]) && (func(xVec[i - 1], m, b) < yBox[1]))))
        {
            crossVec.push_back(std::make_pair(xBox[1], (m * xBox[1] + b)));
        }

        if ((xVec[i] == xVec[i - 1]) && (sign(yVec[i] - yBox[1]) != sign(yVec[i - 1] - yBox[1])))
        {
            crossVec.push_back(std::make_pair(xVec[i], yBox[1]));
        }

        if ((xVec[i] == xVec[i - 1]) && (sign(yVec[i] - yBox[0]) != sign(yVec[i - 1] - yBox[0])))
        {
            crossVec.push_back(std::make_pair(xVec[i], yBox[0]));
        }

        if ((yVec[i] == yVec[i - 1]) && (sign(xVec[i] - xBox[1]) != sign(xVec[i - 1] - xBox[1])))
        {
            crossVec.push_back(std::make_pair(xBox[1], yVec[i]));
        }

        if ((yVec[i] == yVec[i - 1]) && (sign(xVec[i] - xBox[0]) != sign(xVec[i - 1] - xBox[0])))
        {
            crossVec.push_back(std::make_pair(xBox[0], yVec[i]));
        }
    }
    return crossVec;
}