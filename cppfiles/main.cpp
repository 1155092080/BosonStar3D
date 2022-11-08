/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D Two Fluid Hydrodynamics Simulation Code	   !
! Can also be used as simulating type 1a supernova !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include "parameter.h"
#include "./Function/flagConflicts.h"
#include "./Function/initial.h"
#include "./Variable/variable.h"
#include "./Function/fileIO.h"
#include "./Function/output.h"
#include "./Function/findDt.h"
#include "./Function/rungeKutta.h"
#include "./Function/findQuantumPotential.h"
using namespace std;

int main(){
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cout << "Welcome to 1D HydroCode" << endl;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//We first check that whether there are conflicts on flags
	if(flagConflicts()){
		cout << "Simulation End" << endl;
		return -1;
	}
	
	//Allocate all physical variable
	buildVariable();
	
	//Initialize file IO
	openFileHydro();
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//Initialize all variable
	initial();
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//Output initial star profile to file
	outputHydro();
	outputProfileHydro();
	
	//Time solver
	cout << endl << "Starting time evolution" << endl;
	for(int n=0;n<totalTimeStep;n++){
		
		//Find an appropriate dt
		findDt();
		
		//Calling the 3rd order rungekutta time solver
		if(qpFlag){
			findQuantumPotential();
		}else{
			rungeKutta(n);
		}
		//rungeKutta(n);
		
		//Update the global simulation time
		globalTime = globalTime + dt;
		
		//Output simulation info
		if(n%printEvery == 0){
			cout << "Time: " << globalTime << "\t" << "dt: " << dt << "\t" << "Time step: " << n << "\t" << (*rho1)[0] << "\t" << (*rho2)[0] << endl;
			outputAni(n);
		}
		
		//Output profile to file
		if(n%outputEvery == 0){
			if(n != 0){
				outputHydro();
				outputProfileHydro();
			}
		}
		
		//Determine whether simulation ends
		if(globalTime >= totalTime){
			//Output final profile to file
			outputHydro();
			outputProfileHydro();
			break;
		}
		
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//Playing area
	
	cout << endl << "Playground:" << endl;
	
	if(dmFlag){
		for(int i=0;i<lengthStep;i++){
			if((*rho1)[i] <= rhoMinDM){
				cout << "Radius1: " << (i+0.5)*dx << endl;
				break;
			}
		}
	}
	
	for(int i=0;i<lengthStep;i++){
		if((*rho2)[i] <= rhoMinNM){
			cout << "Radius2: " << (i+0.5)*dx << endl;
			break;
		}
	}
	
	double tempMass = 0.0;
	
	if((*rho2)[0] > rho2A){
		tempMass = (4.0/3.0)*3.14159265358979*((0+0.5)*dx)*((0+0.5)*dx)*((0+0.5)*dx)*((*rho2)[0]);
	}
	
	//Do the same of the remaining grid
	for(int i=1;i<lengthStep;i++){
		if((*rho2)[i] > 1.01*rho2A){
			tempMass += 4.0*3.14159265358979*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i];
		}
	}
	
	double tempMass2 = 0.0;
	
	if((*rho2)[0] > rho2A){
		tempMass2 = (4.0/3.0)*3.14159265358979*((0+0.5)*dx)*((0+0.5)*dx)*((0+0.5)*dx)*((*rho2)[0]);
	}
	
	//Do the same of the remaining grid
	for(int i=1;i<lengthStep;i++){
		if((*rho2)[i] > 1.01*rho2A){
			tempMass2 += 4.0*3.14159265358979*dx*((i+0.5)*dx)*((i+0.5)*dx)*(*rho2)[i];
		}
		if(tempMass2 > 0.99*tempMass){
			cout << "Radius_99: " << (i+0.5)*dx << endl;
			break;
		}
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//Deallocate file IO's memory
	closeFileHydro();
	
	//Deallocate variable's memory
	destroyVariable();
	
	//End
	cout << endl << "Simulation End" << endl;
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	return 0;
}

