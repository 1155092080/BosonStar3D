#ifndef __OPENFILE_H__
#define __OPENFILE_H__
#include <fstream>

void openFileHydro();
void closeFileHydro();

extern std::ofstream foutDensityDM;
extern std::ofstream foutEpsilonDM;
extern std::ofstream foutPressureDM;
extern std::ofstream foutdpdrhoDM;
extern std::ofstream foutdpdepsilonDM;
extern std::ofstream foutCentralDensityDM;
extern std::ofstream foutMassDM;
extern std::ofstream foutEnergyDM;
extern std::ofstream foutMechEnergyDM;
extern std::ofstream foutIntEnergyDM;
extern std::ofstream foutVelocityDM;

extern std::ofstream foutDensityNM;
extern std::ofstream foutEpsilonNM;
extern std::ofstream foutPressureNM;
extern std::ofstream foutVelocityNM;
extern std::ofstream foutdpdrhoNM;
extern std::ofstream foutdpdepsilonNM;
extern std::ofstream foutCentralDensityNM;
extern std::ofstream foutMassNM;
extern std::ofstream foutEnergyNM;
extern std::ofstream foutMechEnergyNM;
extern std::ofstream foutIntEnergyNM;

#endif

