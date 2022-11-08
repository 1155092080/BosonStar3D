#include <fstream>
#include "fileIO.h"
#include "../parameter.h" 
using namespace std;

ofstream foutDensityDM;
ofstream foutEpsilonDM;
ofstream foutPressureDM;
ofstream foutdpdrhoDM;
ofstream foutdpdepsilonDM;
ofstream foutCentralDensityDM;
ofstream foutMassDM;
ofstream foutEnergyDM;
ofstream foutMechEnergyDM;
ofstream foutIntEnergyDM;
ofstream foutVelocityDM;

ofstream foutDensityNM;
ofstream foutEpsilonNM;
ofstream foutPressureNM;
ofstream foutVelocityNM;
ofstream foutdpdrhoNM;
ofstream foutdpdepsilonNM;
ofstream foutCentralDensityNM;
ofstream foutMassNM;
ofstream foutEnergyNM;
ofstream foutMechEnergyNM;
ofstream foutIntEnergyNM;

void openFileHydro(){
	
	if(dmFlag){
		foutDensityDM.open("./Outfile/Star_WENO_Density_DM.dat");
		foutEpsilonDM.open("./Outfile/Star_WENO_Epsilon_DM.dat");
		foutPressureDM.open("./Outfile/Star_WENO_Pressure_DM.dat");
		foutdpdrhoDM.open("./Outfile/Star_WENO_dpdrho_DM.dat");
		foutdpdepsilonDM.open("./Outfile/Star_WENO_dpdepsilon_DM.dat");
		foutCentralDensityDM.open("./Outfile/Star_WENO_CentralDensity_DM.dat");
		foutMassDM.open("./Outfile/Star_WENO_Mass_DM.dat");
		foutEnergyDM.open("./Outfile/Star_WENO_Energy_DM.dat");
		foutMechEnergyDM.open("./Outfile/Star_WENO_MechEnergy_DM.dat");
		foutIntEnergyDM.open("./Outfile/Star_WENO_IntEnergy_DM.dat");
		
		if(runDMFlag){
			foutVelocityDM.open("./Outfile/Star_WENO_Velocity_DM.dat");
		}
		
	}
	
	foutDensityNM.open("./Outfile/Star_WENO_Density_NM.dat");
	foutEpsilonNM.open("./Outfile/Star_WENO_Epsilon_NM.dat");
	foutPressureNM.open("./Outfile/Star_WENO_Pressure_NM.dat");
	foutVelocityNM.open("./Outfile/Star_WENO_Velocity_NM.dat");
	foutdpdrhoNM.open("./Outfile/Star_WENO_dpdrho_NM.dat");
	foutdpdepsilonNM.open("./Outfile/Star_WENO_dpdepsilon_NM.dat");
	foutCentralDensityNM.open("./Outfile/Star_WENO_CentralDensity_NM.dat");
	foutMassNM.open("./Outfile/Star_WENO_Mass_NM.dat");
	foutEnergyNM.open("./Outfile/Star_WENO_Energy_NM.dat");
	foutMechEnergyNM.open("./Outfile/Star_WENO_MechEnergy_NM.dat");
	foutIntEnergyNM.open("./Outfile/Star_WENO_IntEnergy_NM.dat");
	
	if(helmeosFlag){
		//do sth
	}
	
	if(outputPotentialFlag){
		//do sth
	}
	
	return;
}

void closeFileHydro(){
	
	if(dmFlag){
		foutDensityDM.close();
		foutEpsilonDM.close();
		foutPressureDM.close();
		foutdpdrhoDM.close();
		foutdpdepsilonDM.close();
		foutCentralDensityDM.close();
		foutMassDM.close();
		foutEnergyDM.close();
		foutMechEnergyDM.close();
		foutIntEnergyDM.close();
		
		if(runDMFlag){
			foutVelocityDM.close();
		}
	}
	
	foutDensityNM.close();
	foutEpsilonNM.close();
	foutPressureNM.close();
	foutVelocityNM.close();
	foutdpdrhoNM.close();
	foutdpdepsilonNM.close();
	foutCentralDensityNM.close();
	foutMassNM.close();
	foutEnergyNM.close();
	foutMechEnergyNM.close();
	foutIntEnergyNM.close();
			
	return;
}

