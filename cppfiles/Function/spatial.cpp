/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the spatial discritization of the hyperbolic equation !
! Where we have used Lax - Friedrichs flux splitting method     !
! To approximate flux at the boundary                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include "spatial.h"
#include "../Library/parallel.h"
#include "../parameter.h"
#include "../Variable/variable.h"
#include "./wenoModule.h"

#include "../Variable/definition.h"
#include "./boundary.h"
using namespace std;

void alphaSplit(){
	
	double lambda[2];
	double lambdaMax[2] = {0.0,0.0};
	
	for(int j=0;j<lengthStep;j++){
		lambda[1] = abs((*vel2)[j]) + sqrt(abs((*dpdrho2)[j] + ((*p2)[j]/((*rho2)[j]*(*rho2)[j])) * (*dpdepsilon2)[j]));
		lambdaMax[1] = max(lambda[1], lambdaMax[1]);
	}
	
	for(int j=imin2;j<=imax2;j++){
		alpha[j] = lambdaMax[1];	
	}
	
	if(dmFlag && runDMFlag){
		
		for(int j=0;j<lengthStep;j++){
			lambda[0] = abs((*vel1)[j]) + sqrt(abs((*dpdrho1)[j] + ((*p1)[j]/((*rho1)[j]*(*rho1)[j])) * (*dpdepsilon1)[j]));
			lambdaMax[0] = max(lambda[0], lambdaMax[0]);
		}
		
		for(int j=imin1;j<=imax1;j++){
			alpha[j] = lambdaMax[0];	
		}
		
	}
	
	return;
}

void spatial(PhyQua **u){
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//We then allocate the necessary array
	alpha = new double[noOfEq];
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if(movingGridFlag){
		//Skip for now
	}
	
	for(int j=-5;j<lengthStep+5;j++){
		if(dmFlag && runDMFlag){
			
			(*f[irho1])[j] = (*u[irho1])[j] * (*vel1)[j];
			(*f[ivel1])[j] = (*u[ivel1])[j] * (*vel1)[j] + (*p1)[j];
			(*f[itau1])[j] = (*vel1)[j] * ((*u[itau1])[j] + (*p1)[j]);
			
			(*sa[irho1])[j] = (*f[irho1])[j];
			(*sa[ivel1])[j] = (*u[ivel1])[j] * (*vel1)[j];
			(*sa[itau1])[j] = (*f[itau1])[j];
			
			(*sb[irho1])[j] = 0.0;
			(*sb[ivel1])[j] = 0.0;
			(*sb[itau1])[j] = 0.0;
			
			if((*rho1)[j] > rhoMinDM){
				(*sb[ivel1])[j] = (*u[irho1])[j] * (*phip)[j];
				(*sb[itau1])[j] = (*u[ivel1])[j] * (*phip)[j];
			}
			
			(*sc[irho1])[j] = 0.0;
			(*sc[ivel1])[j] = 0.0;
			(*sc[itau1])[j] = 0.0;
			
			if(j>-1 && j<lengthStep && (*rho1)[j]>rhoMinDM && qpDMFlag){
				(*sc[ivel1])[j] = (*u[irho1])[j]/mBoson * (-(*quantumPotential1)[j+2] + 8.0*(*quantumPotential1)[j+1] - 8.0*(*quantumPotential1)[j-1] + (*quantumPotential1)[j-2])/(12.0*dx);
			}
			
			/*if(j>-1 && j<lengthStep && (*rho1)[j]>rhoMinDM && qpDMFlag){
				//(*sc[ivel2])[j] = (*u[irho2])[j]/mBoson * ((*quantumPotential)[j+3] - 9.0*(*quantumPotential)[j+2] + 45.0*(*quantumPotential)[j+1] - 45.0*(*quantumPotential)[j-1] + 9.0*(*quantumPotential)[j-2] - (*quantumPotential)[j-3])/(60.0*dx);
				(*sc[ivel1])[j] = (*u[irho1])[j]/mBoson * (-(*quantumPotential1)[j+2] + 8.0*(*quantumPotential1)[j+1] - 8.0*(*quantumPotential1)[j-1] + (*quantumPotential1)[j-2])/(12.0*dx);
			}else if(j>-1 && j<lengthStep && (*rho1)[j-1]>rhoMinDM && qpDMFlag){
				(*sc[ivel1])[j-1] = (*u[irho1])[j-1]/mBoson * (147.0*(*quantumPotential1)[j-1] - 360.0*(*quantumPotential1)[j-2] + 450.0*(*quantumPotential1)[j-3] - 400.0*(*quantumPotential1)[j-4] + 225.0*(*quantumPotential1)[j-5] - 72.0*(*quantumPotential1)[j-6] + 10.0*(*quantumPotential1)[j-7])/(60.0*dx);
				(*sc[ivel1])[j-2] = (*u[irho1])[j-2]/mBoson * (147.0*(*quantumPotential1)[j-2] - 360.0*(*quantumPotential1)[j-3] + 450.0*(*quantumPotential1)[j-4] - 400.0*(*quantumPotential1)[j-5] + 225.0*(*quantumPotential1)[j-6] - 72.0*(*quantumPotential1)[j-7] + 10.0*(*quantumPotential1)[j-8])/(60.0*dx);
			}*/
			
			
			if(levelSetFlagDM){
				//Skip for now
			}
			
		}
		
		//Now we do the same thing for NM
		
		(*f[irho2])[j] = (*u[irho2])[j] * (*vel2)[j];
		(*f[ivel2])[j] = (*u[ivel2])[j] * (*vel2)[j] + (*p2)[j] + (*quantumPotential)[j];
		(*f[itau2])[j] = (*vel2)[j] * ((*u[itau2])[j] + (*p2)[j]);
		
		(*sa[irho2])[j] = (*f[irho2])[j];
		(*sa[ivel2])[j] = (*u[ivel2])[j] * (*vel2)[j] + (*quantumPotential)[j] - (*quantumPotential1)[j];
		(*sa[itau2])[j] = (*f[itau2])[j];
		
		(*sb[irho2])[j] = 0.0;
		(*sb[ivel2])[j] = 0.0;
		(*sb[itau2])[j] = 0.0;
		
		if((*rho2)[j] > rhoMinNM){
			(*sb[ivel2])[j] = (*u[irho2])[j] * (*phip)[j];
			(*sb[itau2])[j] = (*u[ivel2])[j] * (*phip)[j];
		}
		
		(*sc[irho2])[j] = 0.0;
		(*sc[ivel2])[j] = 0.0;
		(*sc[itau2])[j] = 0.0;
		
		/*Work
		if(j>-1 && j<lengthStep && (*rho2)[j] > rhoMinNM && qpFlag){
			(*sc[ivel2])[j] = (*u[irho2])[j]/mBoson * (-(*quantumPotential)[j+2] + 8.0*(*quantumPotential)[j+1] - 8.0*(*quantumPotential)[j-1] + (*quantumPotential)[j-2])/(12.0*dx);
		}*/
		
		/*if(j>-1 && j<lengthStep && (*rho2)[j] > rhoMinNM && qpFlag){
			(*sc[ivel2])[j] = (*u[irho2])[j]/mBoson * (-(*quantumPotential)[j+2] + 8.0*(*quantumPotential)[j+1] - 8.0*(*quantumPotential)[j-1] + (*quantumPotential)[j-2])/(12.0*dx);
		}else if(j>-1 && j<lengthStep && (*rho2)[j-1]> rhoMinNM && qpFlag){
			(*sc[ivel2])[j-1] = (*u[irho2])[j-1]/mBoson * (147.0*(*quantumPotential)[j-1] - 360.0*(*quantumPotential)[j-2] + 450.0*(*quantumPotential)[j-3] - 400.0*(*quantumPotential)[j-4] + 225.0*(*quantumPotential)[j-5] - 72.0*(*quantumPotential)[j-6] + 10.0*(*quantumPotential)[j-7])/(60.0*dx);
			(*sc[ivel2])[j-2] = (*u[irho2])[j-2]/mBoson * (147.0*(*quantumPotential)[j-2] - 360.0*(*quantumPotential)[j-3] + 450.0*(*quantumPotential)[j-4] - 400.0*(*quantumPotential)[j-5] + 225.0*(*quantumPotential)[j-6] - 72.0*(*quantumPotential)[j-7] + 10.0*(*quantumPotential)[j-8])/(60.0*dx);
		}*/
		
		if(levelSetFlagNM){
			//Skip for now
		}
		
		if(flameFlag){
			//Skip for now
		}
		
		if(xisotranFlag){
			//Skip for now
		}
		
	}
	
	if(movingGridFlag){
		//Skip for now
	}
	
	alphaSplit();
	
	parallel_for(0, noOfEq, [](int i, PhyQua **u){
		
		//Allocate
		PhyQua *vFirst = new PhyQua(lengthStep);
		PhyQua *vPFirst = new PhyQua(lengthStep);
		PhyQua *vMFirst = new PhyQua(lengthStep);
		
		PhyQua *vSecond = new PhyQua(lengthStep);
		PhyQua *vPSecond = new PhyQua(lengthStep);
		PhyQua *vMSecond = new PhyQua(lengthStep);
		
		PhyQua *fluxP = new PhyQua(lengthStep);
		PhyQua *fluxM = new PhyQua(lengthStep);
		
		//////////////////////
		
		for(int j=-5;j<lengthStep+5;j++){
			(*vFirst)[j] = 0.5*((*f[i])[j] + alpha[i]*(*u[i])[j]);
		}
		
		weno(vFirst, vPFirst, vMFirst);
		
		for(int j=-5;j<lengthStep+5;j++){
			(*fluxP)[j] = (*vMFirst)[j];
		}
		
		for(int j=-5;j<lengthStep+5;j++){
			(*vSecond)[j] = 0.5*((*f[i])[j] - alpha[i]*(*u[i])[j]);
		}
		
		weno(vSecond, vPSecond, vMSecond);
		
		for(int j=-5;j<lengthStep+4;j++){
			(*fluxM)[j] = (*vPSecond)[j+1];
		}
		
		//we sum up the splitted flux
		for(int j=-5;j<lengthStep+5;j++){
			(*flux[i])[j] = (*fluxP)[j] + (*fluxM)[j];
			if(i<=imax1 && runDMFlag){
				if((*rho1)[j] <= rhoMinDM && (*flux[i])[j] > 0.0){
					(*flux[i])[j] = 0.0;
				}
			}else{
				if((*rho2)[j] <= rhoMinNM && (*flux[i])[j] > 0.0){
					(*flux[i])[j] = 0.0;
				}
			}
		}
		
		//////////////////////
		
		//Deallocate
		delete vFirst;
		delete vPFirst;
		delete vMFirst;
		
		delete vSecond;
		delete vPSecond;
		delete vMSecond;
		
		delete fluxP;
		delete fluxM;
		
		for(int j=0;j<lengthStep;j++){
			(*dfdx[i])[j] = ((*flux[i])[j] - (*flux[i])[j-1])/dx;
		}
		
		if(movingGridFlag){
			//Skip for now
		}
		
		for(int j=0;j<lengthStep;j++){
			(*l[i])[j] = -(*dfdx[i])[j] - spDimI / ((j+0.5)*dx) * (*sa[i])[j] - (*sb[i])[j] - (*sc[i])[j];
		}
		
	}, u);
	
	////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	/*for(int j=0;j<lengthStep;j++){
		cout << j << "\t" << (*rho2)[j] << "\t" << -(*dfdx[ivel2])[j] - spDimI / ((j+0.5)*dx) * (*sa[ivel2])[j] - (*sb[ivel2])[j] << "\t" << -(*sc[ivel2])[j] << endl;
	}
	cout << endl;*/
	
	//We deallocate the necessary array
	delete [] alpha;
	
	return;
}

