/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the specific energy density epsilon      !
! which is necessary for hydro dynamics evolution. We assume     !
! completely degenerate fermi gas EOS for DM and either          !
! finite temperature EOS or completely degenerate EOS for NM     !
! If you want to use your own EOS, you need to take care of this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <cmath>
#include "getEpsilon.h"
#include "../Variable/definition.h"
#include "../Variable/variable.h"
#include "../parameter.h"

void getEpsilon(){
	
	if(dmFlag){
		for(int j=0;j<lengthStep;j++){
			(*epsilon1)[j] = k1*pow((*rho1)[j],(gamma1-1.0))/(gamma1-1.0);
		}
		
		epsilon1A = (*epsilon1)[lengthStep-1];
	}
	
	if(helmeosFlag){
		//Skip first
	}else{
		if(fermiEOS){
			for(int j=-5;j<lengthStep+5;j++){
				double dlfmmo = pow(6.0*PI*PI*hBar*hBar*hBar*(*rho2)[j]*ye2/(gs2*mb2*me2*me2*me2), 1.0/3.0);
				
				if(dlfmmo <= 6.0e-2){
					(*epsilon2)[j] = (me2*ye2/mb2)*(0.3*dlfmmo*dlfmmo - 3.0/56.0 * pow(dlfmmo, 4.0) + 3.0/144.0 * pow(dlfmmo, 6.0));
				}else{
					(*epsilon2)[j] = (3.0*me2*ye2/(8.0*mb2*dlfmmo*dlfmmo*dlfmmo)) * (dlfmmo*sqrt(1.0+dlfmmo*dlfmmo)*(1.0+2.0*dlfmmo*dlfmmo)-log(dlfmmo+sqrt(1.0+dlfmmo*dlfmmo))) - (me2*ye2/mb2);
				}
				
			}
		}else if(bosonEOS){
			for(int j=0;j<lengthStep;j++){
				double K  = 8.0*PI*scattering*hBar*hBar/(4.0*mBoson*mBoson*mBoson);
				(*epsilon2)[j] = (1.0/(36.0*K)) * (2.0* (sqrt((*rho2)[j]*(12.0*K+(1.0/(*rho2)[j])))-1.0) /(*rho2)[j] - 12.0*K*log(1.0/(*rho2)[j]) + 12.0*K*log(2.0* (sqrt((*rho2)[j]*(12.0*K+(1.0/(*rho2)[j])))+1.0) /(*rho2)[j] + 12.0*K));
			}
		}else{
			for(int j=0;j<lengthStep;j++){
				(*epsilon2)[j] = k2*pow((*rho2)[j],(gamma2-1.0))/(gamma2-1.0);
			}
		}
	}
	
	epsilon2A = (*epsilon2)[lengthStep-1];
	
	return;
}

