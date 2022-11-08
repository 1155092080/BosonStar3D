/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the profiles, such as density !
! pressure, mass, epsilon ...etc                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <iostream>
#include "output.h"
#include "./fileIO.h"
#include "../Variable/variable.h"
#include "../parameter.h"
#include "../Variable/definition.h"
using namespace std;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the central density of NM and DM by interpolation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void findCentralDensity(double &cRho1, double &cRho2){
	
	//Interpolate DM central density
	//We do this only if the users wants DM component
	if(dmFlag){
		cRho1 = (1.5*dx*1.5*dx*(*rho1)[0] - 0.5*dx*0.5*dx*(*rho1)[1])/(1.5*dx*1.5*dx - 0.5*dx*0.5*dx);
	}
	
	cRho2 = (1.5*dx*1.5*dx*(*rho2)[0] - 0.5*dx*0.5*dx*(*rho2)[1])/(1.5*dx*1.5*dx - 0.5*dx*0.5*dx);
	
	return;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine find the total energy (internal + Mechanical) !
! of the star assuming a spherical symmetric geometry	     	!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void findEnergy(double &e1, double &e2, double &me1, double &me2, double &ie1, double &ie2){
	e1 = 0.0;
	e2 = 0.0;
	me1 = 0.0;
	me2 = 0.0;
	ie1 = 0.0;
	ie2 = 0.0;
	double q1 = 0.0;
	double q2 = 0.0;
	
	if(spDimI == 2){
		//The sum of energy is the sum of internal + mechanical
		if((*rho2)[0] > rho2A){
			me2 = (4.0/3.0)*PI*(0.0+0.5)*dx*(0.0+0.5)*dx*(0.0+0.5)*dx*(*rho2)[0]*(0.5*(*vel2)[0]*(*vel2)[0]+0.5*(*phi)[0]);
			ie2 = (4.0/3.0)*PI*(0.0+0.5)*dx*(0.0+0.5)*dx*(0.0+0.5)*dx*(*rho2)[0]*(*epsilon2)[0];
			q2 = (4.0/3.0)*PI*(0.0+0.5)*dx*(0.0+0.5)*dx*(0.0+0.5)*dx*(*rho2)[0]*(*quantumPotential)[0]/mBoson;
		}
		
		for(int i=1;i<lengthStep;i++){
			if((*rho2)[i] > rho2A){
				me2 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i]*(0.5*(*vel2)[i]*(*vel2)[i]+0.5*(*phi)[i]);
				ie2 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i]*(*epsilon2)[i];
				q2 +=  4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i]*(*quantumPotential)[i]/mBoson;
			}
		}
		
		e2 = me2 + ie2 + q2;
		
		if(dmFlag && runDMFlag){
			
			//The sum of energy is the sum of internal + mechanical
			if((*rho1)[0] > rho1A){
				me1 = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho1)[0]*(0.5*(*vel1)[0]*(*vel1)[0]+0.5*(*phi)[0]);
				ie1 = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho1)[0]*(*epsilon1)[0];
			}
			
			for(int i=1;i<lengthStep;i++){
				if((*rho1)[i] > rho1A){
					me1 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho1)[i]*(0.5*(*vel1)[i]*(*vel1)[i] + 0.5*(*phi)[i]);
					ie1 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho1)[i]*(*epsilon1)[i];
				}
			}
			
			e1 = me1 + ie1;
			
		}else if(dmFlag && !runDMFlag){
			
			//The sum of energy is the sum of internal + mechanical
			if((*rho1)[0] > rho1A){
				me1 = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho1)[0]*0.5*(*phi)[0];
				ie1 = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho1)[0]*(*epsilon1)[0];
			}
			
			for(int i=1;i<lengthStep;i++){
				if((*rho1)[i] > rho1A){
					me1 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho1)[i]*0.5*(*phi)[i];
					ie1 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho1)[i]*(*epsilon1)[i];
				}
			}
			
			e1 = me1 + ie1;
			
		}
		
	}
	
	return;
}

void findMass(double &m1, double &m2){
	//Find the accumulated mass for first grid
	m2 = 0.0;
	
	if((*rho2)[0] > rho2A){
		m2 = (4.0/3.0)*PI*((0+0.5)*dx)*((0+0.5)*dx)*((0+0.5)*dx)*((*rho2)[0]);
	}
	
	//Do the same of the remaining grid
	for(int i=1;i<lengthStep;i++){
		if((*rho2)[i] > 1.01*rho2A){
			m2 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i];
		}
	}
	
	//Now do the DM case only if the users wants dm component
	if(dmFlag){
		m1 = 0.0;
		
		if((*rho1)[0] > rho1A){
			m1 = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho1)[0];
		}
		
		for(int i=1;i<lengthStep;i++){
			if((*rho1)[i] > 1.01*rho1A){
				m1 += 4.0*PI*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho1)[i];
			}
		}
		
	}
	
	return;
}

void outputProfileHydro(){
	
	//NM density
	foutDensityNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutDensityNM << dx*(j+0.5) << "\t" << (*rho2)[j] << endl;
	}
	foutDensityNM << endl;
	
	//NM epsilon
	foutEpsilonNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutEpsilonNM << dx*(j+0.5) << "\t" << (*epsilon2)[j] << endl;
	}
	foutEpsilonNM << endl;
	
	//NM dpdrho
	foutdpdrhoNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutdpdrhoNM << dx*(j+0.5) << "\t" << (*dpdrho2)[j] << endl;
	}
	foutdpdrhoNM << endl;
	
	//NM dpdepsilon
	foutdpdepsilonNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutdpdepsilonNM << dx*(j+0.5) << "\t" << (*dpdepsilon2)[j] << endl;
	}
	foutdpdepsilonNM << endl;
	
	//NM velocity
	foutVelocityNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutVelocityNM << dx*(j+0.5) << "\t" << (*vel2)[j] << endl;
	}
	foutVelocityNM << endl;
	
	//NM pressure
	foutPressureNM << "\"Time = " << globalTime << endl;
	for(int j=0;j<lengthStep;j++){
		foutPressureNM << dx*(j+0.5) << "\t" << (*p2)[j] << endl;
	}
	foutPressureNM << endl;
	
	if(helmeosFlag){
		//do sth
	}
	
	if(outputPotentialFlag){
		//do sth
	}
	
	if(dmFlag){
		
		//DM density
		foutDensityDM << "\"Time = " << globalTime << endl;
		for(int j=0;j<lengthStep;j++){
			foutDensityDM << dx*(j+0.5) << "\t" << (*rho1)[j] << endl;
		}
		foutDensityDM << endl;
		
		//DM epsilon
		foutEpsilonDM << "\"Time = " << globalTime << endl;
		for(int j=0;j<lengthStep;j++){
			foutEpsilonDM << dx*(j+0.5) << "\t" << (*epsilon1)[j] << endl;
		}
		foutEpsilonDM << endl;
		
		//DM pressure
		foutPressureDM << "\"Time = " << globalTime << endl;
		for(int j=0;j<lengthStep;j++){
			foutPressureDM << dx*(j+0.5) << "\t" << (*p1)[j] << endl;
		}
		foutPressureDM << endl;
		
		//DM dpdrho
		foutdpdrhoDM << "\"Time = " << globalTime << endl;
		for(int j=0;j<lengthStep;j++){
			foutdpdrhoDM << dx*(j+0.5) << "\t" << (*dpdrho1)[j] << endl;
		}
		foutdpdrhoDM << endl;
		
		//DM dpdepsilon
		foutdpdepsilonDM << "\"Time = " << globalTime << endl;
		for(int j=0;j<lengthStep;j++){
			foutdpdepsilonDM << dx*(j+0.5) << "\t" << (*dpdepsilon1)[j] << endl;
		}
		foutdpdepsilonDM << endl;
		
		if(runDMFlag){
			foutVelocityDM << "\"Time = " << globalTime << endl;
			for(int j=0;j<lengthStep;j++){
				foutVelocityDM << dx*(j+0.5) << "\t" << (*vel1)[j] << endl;
			}
			foutVelocityDM << endl;
		}
		
	}
	
	return;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output the discrete data at each time step !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

void outputHydro(){
	
	double centralRho1, centralRho2;
	double energy1, energy2;
	double mechEnergy1, mechEnergy2;
	double intEnergy1, intEnergy2;
	double mass1, mass2;
	
	findCentralDensity(centralRho1, centralRho2);
	
	findEnergy(energy1, energy2, mechEnergy1, mechEnergy2, intEnergy1, intEnergy2);
	
	//Since we assume spherical symmetry, the mass
	//Can only be found if you set to spherical symmetry
	//If would be nice if you add the other coordinate version
	if(spDimI == 2){
		findMass(mass1, mass2);
	}
	
	if(dmFlag){
		foutCentralDensityDM << globalTime << "\t" << centralRho1 << endl;
		foutMassDM << globalTime << "\t" << mass1 << endl;
		foutEnergyDM << globalTime << "\t" << energy1 << endl;
		foutMechEnergyDM << globalTime << "\t" << mechEnergy1 << endl;
		foutIntEnergyDM << globalTime << "\t" << intEnergy1 << endl;
	}
	
	//We output NM data
	foutCentralDensityNM << globalTime << "\t" << centralRho2 << endl;
	foutMassNM << globalTime << "\t" << mass2 << endl;
	foutEnergyNM << globalTime << "\t" << energy2 << endl;
	foutMechEnergyNM << globalTime << "\t" << mechEnergy2 << endl;
	foutIntEnergyNM << globalTime << "\t" << intEnergy2 << endl;
	
	return;
}

void outputAni(int n){
	
	//Init
	
	ofstream foutAniNM;
	ofstream foutAniDM;
	ofstream foutAniTime;
	
	//Open
	
	foutAniNM.open("./Outfile/Animation/NM_"+to_string(n)+".txt");
	
	if(dmFlag){
		foutAniDM.open("./Outfile/Animation/DM_"+to_string(n)+".txt");
	}

	foutAniTime.open("./Outfile/Animation/Time_"+to_string(n)+".txt");
	
	//Output
	
	for(int j=0;j<lengthStep;j++){
		foutAniNM << dx*(j+0.5) << "\t" << (*rho2)[j] << endl;
	}
	
	if(dmFlag){
		for(int j=0;j<lengthStep;j++){
			foutAniDM << dx*(j+0.5) << "\t" << (*rho1)[j] << endl;
		}
	}

	foutAniTime << globalTime/2.03017e5 << endl;
	
	//Close
	
	foutAniNM.close();
	
	if(dmFlag){
		foutAniDM.close();
	}

	foutAniTime.close();
	
	return;
}
