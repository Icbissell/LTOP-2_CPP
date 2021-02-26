#include "gridmap.h"

int main()
{
    double azimuth;
    std::cout<< "Select an Azimuth for the Wind: " << std::endl;
    std::cin >> azimuth;

    std::clock_t    start;
    start = std::clock();

    GridMap gridMap(azimuth);
    //gridMap.calc_fc(); //uncomment this line if latitude and longitude arrays are included and calculation of Coriolis effects are desired.

    gridMap.LTOP_calc();
    double runTime =  (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);

    std::cout << "The max precipitation (m/a) is : " << gridMap.mostPrec_ma << std::endl;
    std::cout << "The max precipitation (mm/hr) is : " << gridMap.mostPrec_mmhr << std::endl;
    std::cout << "\n------------------------------------------------------ \n";
    std::cout << "          Solutions saved to LTOP2.txt file" << std::endl;

    /* Make plots */

    gridMap.make_plots();

    std::ofstream LTOP2Stats;
    LTOP2Stats.open("LTOP2.txt");
    LTOP2Stats << "Program runtime: " << runTime << " ms \n";
    LTOP2Stats << "--------------------------- Topography File --------------------------- \n";
    LTOP2Stats << "Grid dimensions nX and nY: " << gridMap.n0 << ", " << gridMap.n1 << std::endl;
    LTOP2Stats << "Grid spacing dX and dY: " << gridMap.dX << ", " << gridMap.dY << std::endl;
    LTOP2Stats << "Coriolis frequency at centroid latitude (mrad/s): " << gridMap.fC*1e3 << std::endl;
    LTOP2Stats << "\n---------------------- Solution --------------------\n";
    LTOP2Stats << "coriolis frequency, fC (rad/s): " << gridMap.fC << std::endl;
    LTOP2Stats << "Time constant for hydrometeor fallout, tauF (s): " << gridMap.tauF << std::endl;
    LTOP2Stats << "Average residence time for cloud water, tauC (s): " << gridMap.tauC << std::endl;
    LTOP2Stats << "Scale height for water vapor, hV (m): " << gridMap.hV << std::endl;
    LTOP2Stats << "Wind speed, U (m/s): " << gridMap.U << std::endl;
    LTOP2Stats << "Wind azimuth (degrees): " << gridMap.azimuth << std::endl;
    LTOP2Stats << "Saturated bouyancy frequency, Ns (rad/s): " << gridMap.Ns << std::endl;
    LTOP2Stats << "Sea-level temperature, T0 (K): " << gridMap.T0 << std::endl;
    LTOP2Stats << "Scale height for total density, hRho (m): " << gridMap.hRho << std::endl;
    LTOP2Stats << "Sea-level saturated density, rhoS0 (g/m^3): " << gridMap.rhoS0 << std::endl;
    LTOP2Stats << "Eddy diffusion, kappa (m/s^2): " << gridMap.kappa << std::endl;
    LTOP2Stats << "Average lapse-rate ratio, gammaSat/gammaEnv, gammaRatio: " << gridMap.gammaRatio << std::endl;
    LTOP2Stats << "Total density at sea level, rho0 (g/m^3) " << gridMap.rho0 << std::endl;
    LTOP2Stats << "Mountain-height number (hMax*Ns/U):" << gridMap.max_el*gridMap.Ns/gridMap.U << std::endl;
    LTOP2Stats << "Maximum Precipitation (m/a): " << gridMap.mostPrec_ma << std::endl;
    LTOP2Stats << "Maximum Precipitation (mm/hr): " << gridMap.mostPrec_mmhr << std::endl;
    LTOP2Stats << "Upper limit for Ns (rad/s): " << gridMap.U/gridMap.max_el << std::endl;
    LTOP2Stats.close();

    return 0;
}
