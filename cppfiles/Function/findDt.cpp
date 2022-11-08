/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This subroutine find the maximum time step dt that !
!  statisfy the CFL condition of stability	      	  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <cmath>
#include <algorithm>
#include "findDt.h"
#include "../Variable/variable.h"
#include "../parameter.h"
using namespace std;

void findDt(){
	
	//Initialize by setting an arbitrarily large number
	double dtdm = 1.0e5;
	double dtnm = 1.0e5;
	
	for(int i=0;i<lengthStep;i++){
		if(dmFlag && runDMFlag){
			if((*rho1)[i] > rho1A){
				double lambdaLocdm = abs((*vel1)[i]) + sqrt((*dpdrho1)[i]+(*dpdepsilon1)[i]*(*p1)[i]/((*rho1)[i]*(*rho1)[i]));
				dtdm = min(cfl*dx/lambdaLocdm, dtdm);
			}
		}else if(dmFlag && !runDMFlag){
			if((*rho1)[i] > rho1A){
				double lambdaLocdm = sqrt((*dpdrho1)[i]+(*dpdepsilon1)[i]*(*p1)[i]/((*rho1)[i]*(*rho1)[i]));
				dtdm = min(cfl*dx/lambdaLocdm, dtdm);
			}
		}
		
		if((*rho2)[i] > rho2A){
			double lambdaLocnm = abs((*vel2)[i]) + sqrt(abs((*dpdrho2)[i]+(*dpdepsilon2)[i]*(*p2)[i]/((*rho2)[i]*(*rho2)[i])));
			dtnm = min(cfl*dx/lambdaLocnm, dtnm);
		}
	}
	
	double dtLoc = min(dtnm, dtdm);
	
	if(globalTime + dtLoc > totalTime){
		dtLoc = min(totalTime - globalTime, dtLoc);
	}
	
	dt = dtLoc;
	
	return;
}

