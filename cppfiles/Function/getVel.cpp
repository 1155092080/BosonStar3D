/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine assign the initial velocity to both nm and dm !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include "getVel.h"
#include "../Variable/variable.h"
#include "./boundary.h"
#include "../parameter.h"


void getVel(){
	
	int j1 = lengthStep-1;
	int j2 = lengthStep-1;
	
	//We locate the atmosphere for  NM
	for(int j=0;j<lengthStep;j++){
		if((*rho2)[j]==rho2A){
			j2 = j;
			break;
		}
	}
	
	//We apply a radial velocity perturbation
	for(int j=0;j<=j2;j++){
		(*vel2)[j] = iniVel2*(j+5.0e-1)/(j2+5.0e-1);
	}
	
	boundary1D(vel2, true);
	
	if(dmFlag && runDMFlag){
		for(int j=0;j<lengthStep;j++){
			if((*rho1)[j]==rho1A){
				j1 = j;
				break;
			}
		}
		
		for(int j=0;j<=j1;j++){
			(*vel1)[j] = iniVel1*(j+5.0e-1)/(j1+5.0e-1);
		}
		
		boundary1D(vel1, true);
	}
	
	return;
}

