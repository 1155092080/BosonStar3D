#include <cmath>
#include "findPressure.h"
#include "../Variable/definition.h"
#include "../Variable/variable.h"
#include "../parameter.h"

void findPressure(){
	
	if(dmFlag){
		for(int j=-5;j<lengthStep+5;j++){
			(*p1)[j] = k1*pow((*rho1)[j], gamma1);
			(*dpdrho1)[j] = k1*gamma1*pow((*rho1)[j], gamma1-1.0);
			(*dpdepsilon1)[j] = (*rho1)[j] * (gamma1-1.0);
		}
	}
	
	if(helmeosFlag){
		//do sth
	}else{
		if(!testModel){
			if(fermiEOS){
				for(int j=-5;j<lengthStep+5;j++){
					double dlfmmo = pow(6.0*PI*PI*hBar*hBar*hBar*(*rho2)[j]*ye2/(gs2*mb2*me2*me2*me2), 1.0/3.0);
					
					if(dlfmmo <= 6.0e-2){
						(*p2)[j] = (gs2*me2*me2*me2*me2)/(30.0*PI*PI*hBar*hBar*hBar) * (pow(dlfmmo, 5.0) - 5.0*pow(dlfmmo, 7.0)/14.0 + 5.0*pow(dlfmmo, 9.0)/24.0);
						(*dpdrho2)[j] = me2*ye2/mb2 * (pow(dlfmmo, 4.0)/3.0 - pow(dlfmmo, 6.0)/6.0 + pow(dlfmmo, 8.0)/8.0);
						(*dpdepsilon2)[j] = 0.0;
					}else{
						(*p2)[j] = (gs2*me2*me2*me2*me2)/(16.0*PI*PI*hBar*hBar*hBar) * (dlfmmo*sqrt(1.0+dlfmmo*dlfmmo)*(2.0*dlfmmo*dlfmmo/3.0 - 1.0) + log(dlfmmo+sqrt(1.0+dlfmmo*dlfmmo)));
						(*dpdrho2)[j] = ye2*me2/(3.0*mb2) * (dlfmmo*dlfmmo/sqrt(1.0+dlfmmo*dlfmmo));
						(*dpdepsilon2)[j] = 0.0;
					}
				}
			}else if(bosonEOS){
				for(int j=-5;j<lengthStep+5;j++){
					double K  = 8.0*PI*scattering*hBar*hBar/(4.0*mBoson*mBoson*mBoson);
					(*p2)[j] = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*(*rho2)[j])-1.0) * (sqrt(1.0+12.0*K*(*rho2)[j])-1.0);
					(*dpdrho2)[j] = 1.0/3.0 * (sqrt(1.0+12.0*K*(*rho2)[j])-1.0) / sqrt(1.0+12.0*K*(*rho2)[j]);
					(*dpdepsilon2)[j] = 0.0;
				}
			}else{
				for(int j=-5;j<lengthStep+5;j++){
					(*p2)[j] = k2*pow((*rho2)[j], gamma2);
					(*dpdrho2)[j] = k2*gamma2*pow((*rho2)[j], gamma2-1.0);
					//(*p2)[j] = (gamma2-1.0) * (*rho2)[j] * (*epsilon2)[j];
					//(*dpdrho2)[j] = (gamma2-1.0) * (*epsilon2)[j];
					(*dpdepsilon2)[j] = (*rho2)[j] * (gamma2-1.0);
				}
			}
		}else{
			for(int j=-5;j<lengthStep+5;j++){
				(*p2)[j] = (gamma2-1.0) * (*rho2)[j] * (*epsilon2)[j];
				(*dpdrho2)[j] = (gamma2-1.0) * (*epsilon2)[j];
				(*dpdepsilon2)[j] = (*rho2)[j] * (gamma2-1.0);
			}
		}
	}
	
	return;
}

