/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert the natural variables (like pressure, density...)       !
! to conservative variables in the conservative hyperbolic equation               !
! If you want to add your new physics into the hydro - code, please do it by hand !
! Write down the corresponding conservative equation, and write down the flux     !
! And source term accordingly                                                     ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include "fromxTox.h"
#include "../Variable/variable.h"
#include "../parameter.h"

void fromRVEToU(PhyQua **u){
	
	for(int j=-5;j<lengthStep+5;j++){
		if(dmFlag && runDMFlag){
			(*u[irho1])[j] = (*rho1)[j];
			(*u[ivel1])[j] = (*rho1)[j] * (*vel1)[j];
			(*u[itau1])[j] = (*rho1)[j] * (*epsilon1)[j] + 5.0e-1 * (*rho1)[j] * (*vel1)[j] * (*vel1)[j];
			
			if(levelSetFlagDM){
				//do sth
			}
		}
		
		(*u[irho2])[j] = (*rho2)[j];
		(*u[ivel2])[j] = (*rho2)[j] * (*vel2)[j];
		(*u[itau2])[j] = (*rho2)[j] * (*epsilon2)[j] + 5.0e-1 * (*rho2)[j] * (*vel2)[j] * (*vel2)[j];
		
		if(levelSetFlagNM){
			//do sth
		}
		
		if(flameFlag){
			//do sth
		}
		
		if(xisotranFlag){
			//do sth
		}
	}
	
	return;
}

void fromUToRVE(PhyQua **u){
	
	for(int j=-5;j<lengthStep+5;j++){
		
		if(dmFlag && runDMFlag){
			(*rho1)[j] = (*u[irho1])[j];
			(*vel1)[j] = (*u[ivel1])[j] / (*rho1)[j];
			(*epsilon1)[j] = (*u[itau1])[j] / (*rho1)[j] - 0.5*(*vel1)[j]*(*vel1)[j];
		}
		
		(*rho2)[j] = (*u[irho2])[j];
		(*vel2)[j] = (*u[ivel2])[j] / (*rho2)[j];
		(*epsilon2)[j] = (*u[itau2])[j] / (*rho2)[j] - 0.5*(*vel2)[j]*(*vel2)[j];
	}
	
	return;
}

