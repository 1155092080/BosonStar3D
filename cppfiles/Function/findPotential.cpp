#include <iostream>
#include <cmath>
#include "findPotential.h"
#include "../Class/phyQua.h"
#include "../Variable/definition.h"
#include "../Variable/variable.h"
#include "../parameter.h"
#include "./boundary.h"

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the total gravitating mass !
! Written by Leung Shing Chi in 2016			!
! The subroutine automatically calculates the density	!
! which is responsible for gravity			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void findGravRho(){
	
	for(int i=-5;i<lengthStep+5;i++){
		(*rho)[i] = 0.0;
	}
	
	for(int j=-5;j<lengthStep+5;j++){
		if((*rho2)[j] > rho2A){
			(*rho)[j] += (*rho2)[j];
		}
	}
	
	if(dmFlag){
		for(int j=-5;j<lengthStep+5;j++){
			if((*rho1)[j] > rho1A){
				(*rho)[j] += (*rho1)[j];
			}
		}
	}
	
	return;
}

void findPotential(){
	
	findGravRho();
	
	if(wGravityI){
		PhyQua *temp = new PhyQua(lengthStep);
		
		(*temp)[0] = (4.0/3.0)*PI*((0.0+0.5)*dx)*((0.0+0.5)*dx)*((0.0+0.5)*dx)*(*rho)[0];
		(*temp)[1] = 9.0*(*temp)[0] + (4.0/3.0)*PI*dx*((1.0+0.5)*dx)*((1.0+0.5)*dx)*(*rho)[1];
		
		for(int j=2;j<lengthStep;j++){
			(*temp)[j] = (*temp)[j-2] + 4.0 * PI * (dx/3.0) *
						 ((*rho)[j-2] * (((j-2.0)+0.5)*dx) * (((j-2.0)+0.5)*dx) +
					  	 4.0*(*rho)[j-1]*(((j-1)+0.5)*dx)*(((j-1)+0.5)*dx) +
					  	 (*rho)[j] * ((j+0.5)*dx) * ((j+0.5)*dx));
		}
		
		for(int j=0;j<lengthStep;j++){
			(*phip)[j] = (*temp)[j] / (((j+0.5)*dx)*((j+0.5)*dx));
		}
		
		boundary1D(phip, true);
		
		(*phi)[lengthStep-1] = 0.0;
		
		for(int j=lengthStep-2;j>=0;j--){
			(*phi)[j] = (*phi)[j+1] - ((*phip)[j] + (*phip)[j+1]) * dx/2.0;
		}
		
		boundary1D(phi, false);
		
		//potentialRelax
		PhyQua *phiNew = new PhyQua(lengthStep);
		PhyQua *error = new PhyQua(lengthStep);
		for(int n=0;n<relaxMax;n++){
			for(int j=0;j<lengthStep-1;j++){
				(*phiNew)[j] = 0.5 * ((*phi)[j+1] + (*phi)[j-1]) + 0.5 * dx/((j+0.5)*dx) * ((*phi)[j+1] - (*phi)[j-1]) - 
							   2.0*PI*dx*dx*(*rho)[j];
			}
			(*phiNew)[lengthStep-1] = (*phi)[lengthStep-1];
			boundary1D(phiNew, false);
			
			for(int j=0;j<lengthStep-1;j++){
				(*error)[j] = ((*phiNew)[j] - (*phi)[j])/(*phi)[j];
			}
			(*error)[lengthStep-1] = 0.0;
			boundary1D(error, false);
			
			for(int j=-5;j<lengthStep+5;j++){
				(*phi)[j] = (*phiNew)[j];
			}
			
			for(int j=0;j<lengthStep;j++){
				(*phip)[j] = (-(*phi)[j+2] + 8.0*(*phi)[j+1] - 8.0*(*phi)[j-1] + (*phi)[j-2])/(1.2e1*dx);
			}
			boundary1D(phip, true);
			
			bool needRelax = false;
			
			for(int j=0;j<lengthStep;j++){
				if(abs((*error)[j]) > tolerance){
					needRelax = true;
				}
			}
			
			if(needRelax){
				continue;
			}else{
				break;
			}
			
		}
		
		delete temp;
		delete phiNew;
		delete error;
		
	}else{
		for(int j=-5;j<lengthStep+5;j++){
			(*phi)[j] = 0.0;
			(*phip)[j] = 0.0;
		}
	}
	
	return;
}

