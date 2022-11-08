/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve the initial star model to be simulate !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <iostream>
#include <cmath>
#include "initial.h"
#include "../parameter.h"
#include "../Variable/variable.h"
#include "../Variable/definition.h"
#include "./getRho.h"
#include "./getVel.h"
#include "./getEpsilon.h"
#include "./checkRho.h"
#include "./fromxTox.h"
#include "./update.h"
#include "./runTestModel.h"
#include "./eosTable.h"
#include "./findQuantumPotential.h"
#include "./sedov.h"
using namespace std;

void getPoly(){
	//assign the polytropic index and exponent for NM
	if(relEOS){
		k2 = pow(3.0, 1.0/3.0)*pow(PI, 2.0/3.0)*hBar*pow(ye2, 4.0/3.0)/(4.0*pow(mb2, 4.0/3.0));
		gamma2 = 4.0/3.0;
	}else if(bosonFlag){
		k2 = 2.0*PI*scattering*hBar*hBar/mBoson/mBoson/mBoson;
		gamma2 = 2.0;
	}else{
		k2 = pow(3.0, 2.0/3.0)*pow(PI, 4.0/3.0)*hBar*hBar*pow(ye2, 5.0/3.0)/(5.0*me2*pow(mb2, 5.0/3.0));
		gamma2 = 5.0/3.0;
	}
	
	//assign the polytopic index and exponent for DM if it exist
	if(dmFlag){
		if(bosonDMFlag){
			k1 = 2.0*PI*scattering*hBar*hBar/mBoson/mBoson/mBoson;
			gamma1 = 2.0;
		}else{
			k1 = pow(3.0, 2.0/3.0)*pow(PI, 4.0/3.0)*hBar*hBar*pow(ye1, 5.0/3.0)/(5.0*me1*pow(mb1, 5.0/3.0));
			gamma1 = 5.0/3.0;
		}
	}
	
	return;
}

void initial(){
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//We need to get the polytropic index
	if(!helmeosFlag || dmFlag){
		getPoly();
	}
	
	//Generate EOS Table if necessary
	if(fermiEOS){
		eosTable();
	}
	
	if(bosonEOS || bosonDMEOS){
		eosTableBoson();
	}
	
	if(aprEOS){
		eosTableAPR();
	}
	
	if(newAPREOS){
		eosTableNewAPR();
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//We build variables neccessary for the levelset compositon
	if(levelSetFlagNM || levelSetFlagDM){
		//Skip it for now
	}
	
	if(!testModel){
		//We build variables neccessary for the chemical compositon
		//Initialize all the variables, arrays that related to chemical composition
		//This is needed before solving for the hydrostatic star
		//buildChem();
		
		//We build variables neccessary for the deflagration
		if(flameFlag){
			//Skip it for now
		}
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	/*
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Solve for the initial density, velocity and epsilon         !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	*/
	
	if(!testModel){
		
		cout << endl << "Solving for the initial hydrostatic star" << endl;
		
		//We determine whether the user wants a 1F or 2F star
		if(qpDMFlag && dmFlag){
			getRhoBoson2F();
		}else if(qpFlag){
			getRhoBoson1F();
		}else if(dmFlag){
			getRho2F();
		}else{
			getRho1F();
		}
		
		////////////////////////////
		
		//solve for the initial velocity and epsilon
		getVel();
		getEpsilon();
		//findQuantumPotential();
		
		//We initialize the deflagration
		if(flameFlag){
			//Skip it for now
		}
		
		//Check whether the density reached atmospheric density
		checkRho();
		
		if(xisotranFlag || flameFlag){
			//Skip it for now
		}
		
	}else{
		cout << "Running test model" << endl;
		runTestModel();
	}
	
	if(levelSetFlagNM || levelSetFlagDM){
		//Skip it for now
	}
	
	if(!testModel){
		//We check the density again because there may be some change in density
		checkRho();
	}
	
	//Sedov test
	if(!testModel){
		if(sedovFlag){
			cout << "Sedov" << endl;
			sedov();
		}
	}
	
	//Convert to conservative variables
	fromRVEToU(uNew);
	
	//Update the variables
	update(uNew, uOld, 1, 1);
	
	return;
}

