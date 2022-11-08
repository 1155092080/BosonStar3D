/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-4 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include <iomanip>

#include <iostream>
#include <cmath>
#include "getRho.h"
#include "../parameter.h"
#include "../Variable/definition.h"
#include "../Variable/variable.h"
#include "./chemicalModule.h"
#include "./boundary.h"
#include "./output.h"

using namespace std;

//The number of hydrostatic equation and the improved accuracy
int noOfEqIni;
int more = pow(10, iniAcc);

// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1 
int binarySearch(double *arr[2], int l, int r, double x)
{ 
    if (r >= l) {
        int mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (arr[mid][0] == x) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid][0] > x) 
            return binarySearch(arr, l, mid - 1, x); 
  
        // Else the element can only be present 
        // in right subarray 
        return binarySearch(arr, mid + 1, r, x); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return l;
}

int binarySearch2(double *arr[2], int l, int r, double x)
{
    if (r >= l) {
        int mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (arr[mid][1] == x) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid][1] > x) 
            return binarySearch2(arr, l, mid - 1, x); 
  
        // Else the element can only be present 
        // in right subarray 
        return binarySearch2(arr, mid + 1, r, x); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return l;
}

int binarySearchGeneral(double *arr, int l, int r, double x)
{ 
    if (r >= l) {
        int mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (arr[mid] == x) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid] > x) 
            return binarySearchGeneral(arr, l, mid - 1, x); 
  
        // Else the element can only be present 
        // in right subarray 
        return binarySearchGeneral(arr, mid + 1, r, x); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return l;
}

void iniDer2F(double *der, double x, double *y){
	der[0] = 4.0*PI*x*x*y[2];
	if(TOVFlag){
		if(y[0]+y[3] == 0.0){
			der[1] = 0.0;
		}else{
			der[1] = -(y[0]+y[3])*y[2]/(x*x)*(1.0+y[1]/y[2])*(1.0+4.0*PI*x*x*x*(y[1]+y[4])/(y[0]+y[3]))/(1.0-2.0*(y[0]+y[3])/x);
		}
	}else{
		der[1] = -(y[0]+y[3])*y[2]/(x*x);
	}
	der[2] = 0.0;
	der[3] = 4.0*PI*x*x*y[5];
	if(TOVFlag){
		if(y[0]+y[3] == 0.0){
			der[4] = 0.0;
		}else{
			der[4] = -(y[0]+y[3])*y[5]/(x*x)*(1.0+y[4]/y[5])*(1.0+4.0*PI*x*x*x*(y[1]+y[4])/(y[0]+y[3]))/(1.0-2.0*(y[0]+y[3])/x);
		}
	}else{
		der[4] = -(y[0]+y[3])*y[5]/(x*x);
	}
	der[5] = 0.0;
	return;
}

void iniDer1F(double *der, double x, double *y){
	der[0] = 4.0*PI*x*x*y[2];
	if(TOVFlag){
		if(y[0] == 0.0){
			der[1] = 0.0;
		}else{
			der[1] = -(y[0]*y[2])/(x*x)*(1.0+y[1]/y[2])*(1.0+4.0*PI*x*x*x*y[1]/y[0])/(1.0-2.0*y[0]/x);
		}
	}else{
		der[1] = -(y[0]*y[2])/(x*x);
	}
	der[2] = 0.0;
	return;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the pressure given the density    !
! of a completely degenerate ideal fermi gas. It is used only !
! in solving for the initial star model. Do not confuse it    !
! with the subroutine findpressure, which aims to update the  !
! pressure after hydrodynamic evolution                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void getRhoEOSRToP(double &iniP, double iniRho, int type){
	if(type == 1){
		if(bosonDMEOS){
			double K  = 8.0*PI*scattering*hBar*hBar/4.0/mBoson/mBoson/mBoson;
			iniP = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*iniRho)-1.0) * (sqrt(1.0+12.0*K*iniRho)-1.0);
		}else{
			iniP = k1*pow(iniRho, gamma1);
		}
	}else if(type == 2){
		if(fermiEOS){
			double dlfmmo = pow(6.0*PI*PI*hBar*hBar*hBar*iniRho*ye2/(gs2*mb2*me2*me2*me2), 1.0/3.0);
			
			if(dlfmmo <= 6.0e-2){
				iniP = (gs2*me2*me2*me2*me2)/(30.0*PI*PI*hBar*hBar*hBar) * (pow(dlfmmo, 5.0) - 5.0*pow(dlfmmo, 7.0)/14.0 + 5.0*pow(dlfmmo, 9.0)/24.0);				
			}else{
				iniP = (gs2*me2*me2*me2*me2)/(16.0*PI*PI*hBar*hBar*hBar) * (dlfmmo*sqrt(1.0+dlfmmo*dlfmmo) * (2.0*dlfmmo*dlfmmo/3.0 - 1.0) + log(dlfmmo + sqrt(1.0+dlfmmo*dlfmmo)));
			}
			
		}else if(bosonEOS){
			double K  = 8.0*PI*scattering*hBar*hBar/4.0/mBoson/mBoson/mBoson;
			iniP = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*iniRho)-1.0) * (sqrt(1.0+12.0*K*iniRho)-1.0);
		}else if(aprEOS || newAPREOS){
			if(iniRho < eosTableList[1][1] || iniRho > eosTableList[16*eosTableWidth-2][1] || std::isnan(iniRho)){
				iniP = 0.0;
				return;
			}
			
			double x1,x2,x3,x4,y1,y2,y3,y4;
			int i = binarySearch2(eosTableList, 0, 16*eosTableWidth-1, iniRho);
			
			x1 = eosTableList[i-2][1];
			x2 = eosTableList[i-1][1];
			x3 = eosTableList[i][1];
			x4 = eosTableList[i+1][1];
			y1 = eosTableList[i-2][0];
			y2 = eosTableList[i-1][0];
			y3 = eosTableList[i][0];
			y4 = eosTableList[i+1][0];
			
			iniP = (iniRho-x2)*(iniRho-x3)*(iniRho-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
				   (iniRho-x1)*(iniRho-x3)*(iniRho-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
				   (iniRho-x1)*(iniRho-x2)*(iniRho-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
				   (iniRho-x1)*(iniRho-x2)*(iniRho-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
			
		}else{
			iniP = k2*pow(iniRho, gamma2);
		}
	}
	
	return;
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine aims at finding the density corresponds !
! to a given pressure in the case of ideal degenerate     !
! fermi gas equation of state by interpolation            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void getRhoEOSPToR(double &iniRho, double iniP, int type){
	if(type == 1){
		if(bosonDMEOS){
			if(iniP < eosTableListBoson[1][0] || iniP > eosTableListBoson[16*eosTableWidth-2][0] || std::isnan(iniP)){
				iniRho = 0.0;
				return;
			}
			
			double x1,x2,x3,x4,y1,y2,y3,y4;
			int i = binarySearch(eosTableListBoson, 0, 16*eosTableWidth-1, iniP);
			
			x1 = eosTableListBoson[i-2][0];
			x2 = eosTableListBoson[i-1][0];
			x3 = eosTableListBoson[i][0];
			x4 = eosTableListBoson[i+1][0];
			y1 = eosTableListBoson[i-2][1];
			y2 = eosTableListBoson[i-1][1];
			y3 = eosTableListBoson[i][1];
			y4 = eosTableListBoson[i+1][1];
			
			iniRho = (iniP-x2)*(iniP-x3)*(iniP-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
					 (iniP-x1)*(iniP-x3)*(iniP-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
					 (iniP-x1)*(iniP-x2)*(iniP-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
					 (iniP-x1)*(iniP-x2)*(iniP-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
		}else{
			iniRho = pow((iniP/k1), (1.0/gamma1));
		}
	}else if(type == 2){
		if(fermiEOS || aprEOS || newAPREOS){
			if(iniP < eosTableList[1][0] || iniP > eosTableList[16*eosTableWidth-2][0] || std::isnan(iniP)){
				iniRho = 0.0;
				return;
			}
			
			double x1,x2,x3,x4,y1,y2,y3,y4;
			int i = binarySearch(eosTableList, 0, 16*eosTableWidth-1, iniP);
			
			x1 = eosTableList[i-2][0];
			x2 = eosTableList[i-1][0];
			x3 = eosTableList[i][0];
			x4 = eosTableList[i+1][0];
			y1 = eosTableList[i-2][1];
			y2 = eosTableList[i-1][1];
			y3 = eosTableList[i][1];
			y4 = eosTableList[i+1][1];
			
			iniRho = (iniP-x2)*(iniP-x3)*(iniP-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
					 (iniP-x1)*(iniP-x3)*(iniP-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
					 (iniP-x1)*(iniP-x2)*(iniP-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
					 (iniP-x1)*(iniP-x2)*(iniP-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
			
		}else if(bosonEOS){
			if(iniP < eosTableListBoson[1][0] || iniP > eosTableListBoson[16*eosTableWidth-2][0] || std::isnan(iniP)){
				iniRho = 0.0;
				return;
			}
			
			double x1,x2,x3,x4,y1,y2,y3,y4;
			int i = binarySearch(eosTableListBoson, 0, 16*eosTableWidth-1, iniP);
			
			x1 = eosTableListBoson[i-2][0];
			x2 = eosTableListBoson[i-1][0];
			x3 = eosTableListBoson[i][0];
			x4 = eosTableListBoson[i+1][0];
			y1 = eosTableListBoson[i-2][1];
			y2 = eosTableListBoson[i-1][1];
			y3 = eosTableListBoson[i][1];
			y4 = eosTableListBoson[i+1][1];
			
			iniRho = (iniP-x2)*(iniP-x3)*(iniP-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
					 (iniP-x1)*(iniP-x3)*(iniP-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
					 (iniP-x1)*(iniP-x2)*(iniP-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
					 (iniP-x1)*(iniP-x2)*(iniP-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
		}else{
			iniRho = pow((iniP/k2), (1.0/gamma2));
		}
	}
	
	return;
}

void getRho2F(){
	
	noOfEqIni = 6;
	
	//The extra array arising from the extra accuracy
	int lengthMoreStep = lengthStep*more;
	//The smaller length step arising from the extra accuracy
	double dxMore = dxOrg/more;
	
	//Quantity related to chemical composition if
	//you are using variable compositions of star
	double *xIsoIni = new double[totalIon];
	
	//This is necessary for any finite temperature EOS since the deriative of pressure plays a role
	double *dpdrho2Temp = new double[lengthMoreStep];
	
	//Variables essential in the RK-4 method
	double *yZero = new double[noOfEqIni];
	double *yOne = new double[noOfEqIni];
	double *yTwo = new double[noOfEqIni];
	double *yThree = new double[noOfEqIni];
	double *der = new double[noOfEqIni];
	
	double **y = new double *[noOfEqIni];
	for(int i=0;i<noOfEqIni;i++){
		y[i] = new double[lengthMoreStep];
	}
	
	double p2C, p2A, p1C, p1A;
	
	//We assign initial chemical composition accordingly
	xIsoIni[cc12] = xc12Ini;
	xIsoIni[cni56] = xni56Ini;
	
	//At the atmosphere, there should not be any composition
	double xISOA = 0.0;
	
	//Now convert the composition into mean atomic and mass number
	azBar(xIsoIni, aBar2Ini, zBar2Ini);
	
	//We choose which table is needed to read according to the chosen EOS
	if(helmeosFlag){
		//do sth
	}else{
		getRhoEOSRToP(p2C, rho2C, 2);
		getRhoEOSRToP(p2A, rho2A, 2);
	}
	getRhoEOSRToP(p1C, rho1C, 1);
	getRhoEOSRToP(p1A, rho1A, 1);
	
	//We assign the value of y at center
	y[0][0] = 0.0;
	y[1][0] = p1C;
	y[2][0] = rho1C;
	y[3][0] = 0.0;
	y[4][0] = p2C;
	y[5][0] = rho2C;
	
	//this is the main structure of RK-4 method
	
	for(int j=0;j<lengthMoreStep-1;j++){
		
		double iniP2, iniRho2, iniP1, iniRho1;
		
		//k1
		for(int i=0;i<noOfEqIni;i++){
			yZero[i] = y[i][j];
		}
		
		double x = j*dxMore;
		
		iniDer2F(der, x, yZero);
		
		//We artificially make the value to zero to avoid singularity at center
		if(j==0){
			der[1] = 0.0;
			der[4] = 0.0;
		}
		
		double k1[6];
		
		for(int i=0;i<noOfEqIni;i++){
			k1[i] = dxMore*der[i];
			yOne[i] = y[i][j]+k1[i]/2.0;
		}
		
		iniP1 = yOne[1];
		iniP2 = yOne[4];
		if(iniP1 <= p1A){
			iniP1 = p1A;
		}
		if(iniP2 <= p2A){
			iniP2 = p2A;
		}
		getRhoEOSPToR(iniRho1, iniP1, 1);
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yOne[2] = iniRho1;
		yOne[5] = iniRho2;
		
		//k2
		
		x = j*dxMore+(1.0/2.0)*dxMore;
		
		iniDer2F(der, x, yOne);
		
		double k2[6];
		
		for(int i=0;i<noOfEqIni;i++){
			k2[i] = dxMore*der[i];
			yTwo[i] = y[i][j]+k2[i]/2.0;
		}
		
		iniP1 = yTwo[1];
		iniP2 = yTwo[4];
		if(iniP1 <= p1A){
			iniP1 = p1A;
		}
		if(iniP2 <= p2A){
			iniP2 = p2A;
		}
		getRhoEOSPToR(iniRho1, iniP1, 1);
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yTwo[2] = iniRho1;
		yTwo[5] = iniRho2;
		
		//k3
		
		x = j*dxMore+(1.0/2.0)*dxMore;
		
		iniDer2F(der, x, yTwo);
		
		double k3[6];
		for(int i=0;i<noOfEqIni;i++){
			k3[i] = dxMore*der[i];
			yThree[i] = y[i][j]+k3[i];
		}
		
		iniP1 = yThree[1];
		iniP2 = yThree[4];
		if(iniP1 <= p1A){
			iniP1 = p1A;
		}
		if(iniP2 <= p2A){
			iniP2 = p2A;
		}
		getRhoEOSPToR(iniRho1, iniP1, 1);
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yThree[2] = iniRho1;
		yThree[5] = iniRho2;
		
		//k4
		
		x = j*dxMore+dxMore;
		
		iniDer2F(der, x, yThree);
		
		double k4[6];
		for(int i=0;i<noOfEqIni;i++){
			k4[i] = dxMore*der[i];
			y[i][j+1] = y[i][j] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}
		
		iniP1 = y[1][j+1];
		iniP2 = y[4][j+1];
		if(iniP1 <= p1A){
			iniP1 = p1A;
			y[1][j+1] = p1A;
		}
		if(iniP2 <= p2A){
			iniP2 = p2A;
			y[4][j+1] = p2A;
		}
		getRhoEOSPToR(iniRho1, iniP1, 1);
		getRhoEOSPToR(iniRho2, iniP2, 2);
		
		y[2][j+1] = iniRho1;
		y[5][j+1] = iniRho2;
		
	}
	
	for(int j=0;j<lengthStep;j++){
		(*rho1)[j] = y[2][(j+1)*more - more/2];
		(*rho2)[j] = y[5][(j+1)*more - more/2];
		
		if((*rho1)[j] < rho1A || std::isnan((*rho1)[j])){
			(*rho1)[j] = rho1A;
		}
		
		if((*rho2)[j] < rho2A || std::isnan((*rho2)[j])){
			(*rho2)[j] = rho2A;
		}
	}
	
	boundary1D(rho1, false);
	boundary1D(rho2, false);
	
	//We first assign the appropriate value to the arrays
	for(int i=0;i<lengthStep;i++){
		if((*rho2)[i] > 0.0){
			(*aBar2)[i] = aBar2Ini;
			(*zBar2)[i] = zBar2Ini;
			for(int j=0;j<totalIon;j++){
				(*xIso[j])[i] = xIsoIni[j];
			}
			(*temp2)[i] = temp2Ini;
		}
	}
	
	boundary2DX(xIso);
	boundary1D(aBar2, false);
	boundary1D(zBar2, false);
	boundary1D(temp2, false);
	
	delete [] xIsoIni;
	delete [] dpdrho2Temp;
	delete [] yZero;
	delete [] yOne;
	delete [] yTwo;
	delete [] yThree;
	delete [] der;
	
	for(int i=0;i<noOfEqIni;i++){
		delete [] y[i];
	}
	delete [] y;
	
	return;
}

void getRho1F(){
	
	noOfEqIni = 3;
	
	//The extra array arising from the extra accuracy
	int lengthMoreStep = lengthStep*more;
	//The smaller length step arising from the extra accuracy
	double dxMore = dxOrg/more;
	
	//Quantity related to chemical composition if
	//you are using variable compositions of star
	double *xIsoIni = new double[totalIon];
	
	//This is necessary for any finite temperature EOS since the deriative of pressure plays a role
	double *dpdrho2Temp = new double[lengthMoreStep];
	
	//Variables essential in the RK-4 method
	double *yZero = new double[noOfEqIni];
	double *yOne = new double[noOfEqIni];
	double *yTwo = new double[noOfEqIni];
	double *yThree = new double[noOfEqIni];
	double *der = new double[noOfEqIni];
	
	double **y = new double *[noOfEqIni];
	for(int i=0;i<noOfEqIni;i++){
		y[i] = new double[lengthMoreStep];
	}
	
	double p2C, p2A;
	
	//We assign initial chemical composition accordingly
	xIsoIni[cc12] = xc12Ini;
	xIsoIni[cni56] = xni56Ini;
	
	//At the atmosphere, there should not be any composition
	double xISOA = 0.0;
	
	//Now convert the composition into mean atomic and mass number
	azBar(xIsoIni, aBar2Ini, zBar2Ini);
	
	//We choose which table is needed to read according to the chosen EOS
	if(helmeosFlag){
		//do sth
	}else{
		getRhoEOSRToP(p2C, rho2C, 2);
		getRhoEOSRToP(p2A, rho2A, 2);
	}
	
	//We assign the value of y at center
	y[0][0] = 0.0;
	y[1][0] = p2C;
	y[2][0] = rho2C;
	
	//this is the main structure of RK-4 method
	
	for(int j=0;j<lengthMoreStep-1;j++){
		
		double iniP2, iniRho2;
		
		//k1
		for(int i=0;i<noOfEqIni;i++){
			yZero[i] = y[i][j];
		}
		
		double x = j*dxMore;
		
		iniDer1F(der, x, yZero);
		
		//We artificially make the value to zero to avoid singularity at center
		if(j==0){
			der[1] = 0.0;
		}
		
		double k1[3];
		
		for(int i=0;i<noOfEqIni;i++){
			k1[i] = dxMore*der[i];
			yOne[i] = y[i][j]+k1[i]/2.0;
		}
		
		iniP2 = yOne[1];
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yOne[2] = iniRho2;
		
		//k2
		
		x = j*dxMore+(1.0/2.0)*dxMore;
		
		iniDer1F(der, x, yOne);
		
		double k2[3];
		
		for(int i=0;i<noOfEqIni;i++){
			k2[i] = dxMore*der[i];
			yTwo[i] = y[i][j]+k2[i]/2.0;
		}
		
		iniP2 = yTwo[1];
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yTwo[2] = iniRho2;
		
		//k3
		
		x = j*dxMore+(1.0/2.0)*dxMore;
		
		iniDer1F(der, x, yTwo);
		
		double k3[3];
		for(int i=0;i<noOfEqIni;i++){
			k3[i] = dxMore*der[i];
			yThree[i] = y[i][j]+k3[i];
		}
		
		iniP2 = yThree[1];
		getRhoEOSPToR(iniRho2, iniP2, 2);
		yThree[2] = iniRho2;
		
		//k4
		
		x = j*dxMore+dxMore;
		
		iniDer1F(der, x, yThree);
		
		double k4[3];
		for(int i=0;i<noOfEqIni;i++){
			k4[i] = dxMore*der[i];
			y[i][j+1] = y[i][j] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}
		
		iniP2 = y[1][j+1];
		if(iniP2 <= p2A){
			iniP2 = p2A;
			y[1][j+1] = p2A;
		}
		getRhoEOSPToR(iniRho2, iniP2, 2);
		y[2][j+1] = iniRho2;
		
	}
	
	for(int j=0;j<lengthStep;j++){
		(*rho2)[j] = y[2][(j+1)*more - more/2];
		
		if((*rho2)[j] < rho2A || std::isnan((*rho2)[j])){
			(*rho2)[j] = rho2A;
		}
	}
	
	boundary1D(rho2, false);
	
	//We first assign the appropriate value to the arrays
	for(int i=0;i<lengthStep;i++){
		if((*rho2)[i] > 0.0){
			(*aBar2)[i] = aBar2Ini;
			(*zBar2)[i] = zBar2Ini;
			for(int j=0;j<totalIon;j++){
				(*xIso[j])[i] = xIsoIni[j];
			}
			(*temp2)[i] = temp2Ini;
		}
	}
	
	boundary2DX(xIso);
	boundary1D(aBar2, false);
	boundary1D(zBar2, false);
	boundary1D(temp2, false);
	
	delete [] xIsoIni;
	delete [] dpdrho2Temp;
	delete [] yZero;
	delete [] yOne;
	delete [] yTwo;
	delete [] yThree;
	delete [] der;
	
	for(int i=0;i<noOfEqIni;i++){
		delete [] y[i];
	}
	delete [] y;
	
	return;
}

void iniDerBoson(double *der, double x, double *y){
	der[0] = x*x*y[1];
	der[1] = y[2];
	der[2] = y[3];
	if(TOVFlag){
		if(y[0] == 0.0){
			der[3] = 0.0;
		}else{
			der[3] = 2.0*mu*y[1]*y[2] + 2.0*y[0]*y[1]/x/x * (1.0+mu*mu*y[1]/psi) * (1.0+mu*mu*y[1]*y[1]*x*x*x/y[0]/psi) / (1.0 - 4.0*mu*y[0]/x/psi) + 2.0/y[1]*y[2]*y[3] - 1.0/y[1]/y[1]*y[2]*y[2]*y[2] - 2.0/x*y[3] + 2.0/x/y[1]*y[2]*y[2] + 2.0/x/x*y[2];
			//der[3] = 2.0*mu*y[1]*y[2] + 2.0*y[0]*y[1]/x/x * (1.0+mu*mu*y[1]/psi+mu/y[1]/psi*(y[2]*y[2]/y[1]-y[3]-2.0/x*y[2])) * (1.0+mu*mu*y[1]*y[1]*x*x*x/y[0]/psi+mu*x*x*x/y[0]/psi*(y[2]*y[2]/y[1]-y[3]-2.0/x*y[2])) / (1.0 - 4.0*mu*y[0]/x/psi) + 2.0/y[1]*y[2]*y[3] - 1.0/y[1]/y[1]*y[2]*y[2]*y[2] - 2.0/x*y[3] + 2.0/x/y[1]*y[2]*y[2] + 2.0/x/x*y[2];
		}
	}else{
		der[3] = 2.0*mu*y[1]*y[2] + 2.0*y[0]*y[1]/x/x + 2.0/y[1]*y[2]*y[3] - 1.0/y[1]/y[1]*y[2]*y[2]*y[2] - 2.0/x*y[3] + 2.0/x/y[1]*y[2]*y[2] + 2.0/x/x*y[2];
	}
	
	return;
}

void getRhoBoson1F(){
	
	noOfEqIni = 4;
	
	psi = 4.0*scattering/mBoson;
	
	//Step Number
	int stepNum = 100000;
	double dxrk = 30.0/(stepNum-1.0);
	
	//Init
	double *yZero = new double[noOfEqIni];
	double *yOne = new double[noOfEqIni];
	double *yTwo = new double[noOfEqIni];
	double *yThree = new double[noOfEqIni];
	double *der = new double[noOfEqIni];

	double **y = new double *[noOfEqIni];
	for(int i=0;i<noOfEqIni;i++){
		y[i] = new double[stepNum];
	}
	
	double *shell = new double[stepNum];
	double *tempRho = new double [stepNum];
	double *r = new double [stepNum];
	
	int farPoint = -1;
	double minValue = 1.01;
	
	/////////////////////
	
	double start = -15.0;
	double end = -0.0001;
	
	while(abs(start-end) > 1e-15){
		
		double A[7];
		double dA = (end-start)/6.0;
		for(int i=0;i<7;i++){
			A[i] = dA*i + start;
		}
		
		double maxFarPoint = -1.0;
		int maxFarPointAA = -1;
		
		for(int aa=0;aa<7;aa++){
			
			y[0][0] = 0.0;
			y[1][0] = 1.0;
			y[2][0] = 0.0;
			y[3][0] = A[aa];
			
			for(int j=0;j<stepNum;j++){
				
				//k1
				for(int i=0;i<noOfEqIni;i++){
					yZero[i] = y[i][j];
				}
				
				double x = j*dxrk;
				
				iniDerBoson(der, x, yZero);
				
				if(j==0){
					der[3] = 0.0;
				}
				
				double k1[4];
				
				for(int i=0;i<noOfEqIni;i++){
					k1[i] = dxrk*der[i];
					yOne[i] = y[i][j]+k1[i]/2.0;
				}
				
				//k2
				
				x = j*dxrk+(1.0/2.0)*dxrk;
				
				iniDerBoson(der, x, yOne);
				
				double k2[4];
				
				for(int i=0;i<noOfEqIni;i++){
					k2[i] = dxrk*der[i];
					yTwo[i] = y[i][j]+k2[i]/2.0;
				}
				
				//k3
				
				x = j*dxrk+(1.0/2.0)*dxrk;
				
				iniDerBoson(der, x, yTwo);
				
				double k3[4];
				for(int i=0;i<noOfEqIni;i++){
					k3[i] = dxrk*der[i];
					yThree[i] = y[i][j]+k3[i];
				}
				
				//k4
				
				x = j*dxrk+dxrk;
				
				iniDerBoson(der, x, yThree);
				
				double k4[4];
				for(int i=0;i<noOfEqIni;i++){
					k4[i] = dxrk*der[i];
					y[i][j+1] = y[i][j] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
				}
				
			}
			
			farPoint = -1;
			minValue = 1.01;
			for(int qq=0;qq<stepNum;qq++){
				if(y[1][qq] < minValue){
					farPoint = qq;
					minValue = y[1][qq];
				}else{
					break;
				}
			}
			
			if(farPoint > maxFarPoint){
				maxFarPoint = farPoint;
				maxFarPointAA = aa;
			}
			
		}
		
		start = A[max(maxFarPointAA-1, 0)];
		end = A[min(maxFarPointAA+1, 6)];
	}
	
	//[Post processing]
	
	for(int qq=0;qq<stepNum;qq++){
		shell[qq] = y[1][qq]*qq*dxrk*qq*dxrk;
	}
	
	double volume = shell[0] + shell[farPoint-100];
	for(int i=1;i<farPoint-100;i++){
		if(i%3 == 0){
			volume += 2.0*shell[i];
		}else{
			volume += 3.0*shell[i];
		}
	}
	volume = 3.0*dxrk*volume/8.0;
	
	double n0 = pow(1.0/volume, 4.0);
	double chi = mu/sqrt(n0);
	double mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
	double b = hBar*hBar/2.0/mStar/mBoson/mBoson;
	double rho0 = mStar*n0/4.0/PI/b/b/b;
	
	double scale = b/pow(n0, 1.0/4.0);
	
	cout << endl;
	cout << "Results of the fourth order ODE solver:" << endl;
	cout << "farPoint: " << (double)farPoint/(double)stepNum << endl;
	cout << "A2: " << start << endl;
	cout << "Length scale: " << scale << endl;
	/*if(totalLength/scale < 5.0 || totalLength/scale > 30.0){
		cout << "Too large/small length scale" << endl;
		exit(-1);
	}*/
	
	///////////////
	//[Convert back to original coordinates]
	
	double ratio = rho2A/rho2C;
	rho2C = rho0;
	rho2A = rho2C * ratio;
	rhoMinNM = 1.1 * rho2A;
	
	for(int i=0;i<stepNum;i++){
		
		if(y[1][i] < ratio || i > farPoint || std::isnan(y[1][i])){
			y[1][i] = ratio;
		}
		
		r[i] = i*dxrk * volume * b;
		tempRho[i] = rho0*y[1][i];
	}
	
	for(int i=0;i<lengthStep;i++){
		double xPos = (i+0.5)*dx;
		for(int j=0;j<stepNum;j++){
			if(xPos < r[j]){
				double x1 = r[j-2];
				double x2 = r[j-1];
				double x3 = r[j];
				double x4 = r[j+1];
				double y1 = tempRho[j-2];
				double y2 = tempRho[j-1];
				double y3 = tempRho[j];
				double y4 = tempRho[j+1];
				
				(*rho2)[i] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;

				break;
			}
		}
	}

	double tmp;
	findMass(tmp, M);
	double omegaCU = abs(mBoson*mBoson/scattering/hBar);
	double a0 = sqrt(hBar/mBoson/omegaCU);
	double N = M*1.9891e30/(mBoson/5.02788e-34*0.001);
	delta = 4.0*PI*scattering*N/a0;
	for(int i=0;i<2*lengthStep;i++){
		for(int j=0;j<2*lengthStep;j++){
			for(int k=0;k<2*lengthStep;k++){
				double xPos = sqrt((i+0.5-lengthStep)*dx*(i+0.5-lengthStep)*dx+(j+0.5-lengthStep)*dx*(j+0.5-lengthStep)*dx+(k+0.5-lengthStep)*dx*(k+0.5-lengthStep)*dx);
				int l = binarySearchGeneral(r, 0, stepNum-1, xPos);

				if(l>stepNum-2){
					l = stepNum-2;
				}

				if(l<2){
					l = 2;
				}

				double x1 = r[l-2];
				double x2 = r[l-1];
				double x3 = r[l];
				double x4 = r[l+1];
				double y1 = tempRho[l-2];
				double y2 = tempRho[l-1];
				double y3 = tempRho[l];
				double y4 = tempRho[l+1];
				
				psiSave[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real((xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4);

				psiSave[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real(sqrt(psiSave[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real()/1.619e-18*1000.0/(M*1.9891e30))*pow(Ra, 3.0/2.0));
				psiSave[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag(0.0);
			}
		}
	}

	for(int j=lengthStep;j<2*lengthStep;j++){
		double tmp = abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+j])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+j]);
		tmp /= (Ra*Ra*Ra);
		tmp *= 1.619e-18;
		tmp /= 1000.0;
		tmp *= (M*1.9891e30);
		(*rho2)[j-lengthStep] = tmp;
	}

	//Hot fix
	rhoMinDM = 1.1 * rho1A;
	rhoMinNM = 1.1 * rho2A;
	
	//Hot fix 2
	for(int i=0;i<lengthStep-1;i++){
		if((*rho2)[i+1] >= (*rho2)[i] || (*rho2)[i+1] < 0.0){
			(*rho2)[i+1] = rho2A;
		}
	}
	
	boundary1D(rho2, false);

	//[Deallocate]
	
	delete [] yZero;
	delete [] yOne;
	delete [] yTwo;
	delete [] yThree;
	delete [] der;

	for(int i=0;i<noOfEqIni;i++){
		delete [] y[i];
	}
	delete [] y;
	
	delete [] r;
	delete [] tempRho;
	delete [] shell;
	
	return;
}

/////////////////////////////////////

//RK4 eq.

//g'=h (f'')
double j(double h){
	return h;
}

//h'=i (f''')
double k(double i){
	return i;
}

//f'=g (f')
double l_rk(double g){
	return g;
}

//i'=... (f'''')
double m(double x, double f, double f2, double g, double h, double i){
	return 2.0*f*f + 2.0*f*f2 - 4.0/x *i + 10.0*g*h/f/x - 6.0*g*g*g/f/f/x + 3.0*i*g/f + 2.0*h*h/f - 7.0*g*g*h/f/f + 3.0*g*g*g*g/f/f/f + 2.0*f*mu*(h+2.0*g/x);
}

double der1(double r, double y3){
	return 4.0*PI*r*r*y3;
}

double der2(double r, double y1, double y3, double y4){
	return -((y1+y4)*y3)/(r*r);
}

double der3(){
	return 0.0;
}

void loop(double *soll, double x, double f, double g, double h, double i, double dx){
	double k1 = l_rk(g);
	double l1 = k(i);
	double m1 = j(h);
	double n1 = m(x,f,0.0,g,h,i);
	
	double k2 = l_rk(g+dx/2.0*m1);
	double l2 = k(i+dx/2.0*n1);
	double m2 = j(h+dx/2.0*l1);
	double n2 = m(x+dx/2.0, f+dx/2.0*k1, 0.0, g+dx/2.0*m1, h+dx/2.0*l1, i+dx/2.0*n1);
	
	double k3 = l_rk(g+dx/2.0*m2);
	double l3 = k(i+dx/2.0*n2);
	double m3 = j(h+dx/2.0*l2);
	double n3 = m(x+dx/2.0, f+dx/2.0*k2, 0.0, g+dx/2.0*m2, h+dx/2.0*l2, i+dx/2.0*n2);
	
	double k4 = l_rk(g+dx*m3);
	double l4 = k(i+dx*n3);
	double m4 = j(h+dx*l3);
	double n4 = m(x+dx, f+dx*k3, 0.0, g+dx*m3, h+dx*l3, i+dx*n3);
	
    soll[0] = f+dx/6.0*(k1+2.0*k2+2.0*k3+k4);
    soll[1] = h+dx/6.0*(l1+2.0*l2+2.0*l3+l4);
    soll[2] = g+dx/6.0*(m1+2.0*m2+2.0*m3+m4);
    soll[3] = i+dx/6.0*(n1+2.0*n2+2.0*n3+n4);
    
    return;
}

void loop2(double *soll, double x, double f, double f2, double g, double h, double i, double dx){
	double k1 = l_rk(g);
	double l1 = k(i);
	double m1 = j(h);
	double n1 = m(x,f,f2,g,h,i);
	
	double k2 = l_rk(g+dx/2.0*m1);
	double l2 = k(i+dx/2.0*n1);
	double m2 = j(h+dx/2.0*l1);
	double n2 = m(x+dx/2.0, f+dx/2.0*k1, f2, g+dx/2.0*m1, h+dx/2.0*l1, i+dx/2.0*n1);
	
	double k3 = l_rk(g+dx/2.0*m2);
	double l3 = k(i+dx/2.0*n2);
	double m3 = j(h+dx/2.0*l2);
	double n3 = m(x+dx/2.0, f+dx/2.0*k2, f2, g+dx/2.0*m2, h+dx/2.0*l2, i+dx/2.0*n2);
	
	double k4 = l_rk(g+dx*m3);
	double l4 = k(i+dx*n3);
	double m4 = j(h+dx*l3);
	double n4 = m(x+dx, f+dx*k3, f2, g+dx*m3, h+dx*l3, i+dx*n3);
	
    soll[0] = f+dx/6.0*(k1+2.0*k2+2.0*k3+k4);
    soll[1] = h+dx/6.0*(l1+2.0*l2+2.0*l3+l4);
    soll[2] = g+dx/6.0*(m1+2.0*m2+2.0*m3+m4);
    soll[3] = i+dx/6.0*(n1+2.0*n2+2.0*n3+n4);
    
    return;
}

void loop3(double *sol, double r, double y1, double y2, double y3, double y4, double y5, double dr){
	
	double iniP2, iniRho2;
	
	double k1 = der1(r, y3);
	double l1 = der2(r, y1, y3, y4);
	double m1 = der3();
	
	iniP2 = y2+dr/2.0*l1;
	getRhoEOSPToR(iniRho2, iniP2, 2);
	y3 = iniRho2;
	
	double k2 = der1(r+dr/2.0, y3+dr/2.0*m1);
	double l2 = der2(r+dr/2.0, y1+dr/2.0*k1, y3+dr/2.0*m1, y4);
	double m2 = der3();
	
	iniP2 = y2+dr/2.0*l2;
	getRhoEOSPToR(iniRho2, iniP2, 2);
	y3 = iniRho2;
	
	double k3 = der1(r+dr/2.0, y3+dr/2.0*m2);
	double l3 = der2(r+dr/2.0, y1+dr/2.0*k2, y3+dr/2.0*m2, y4);
	double m3 = der3();
	
	iniP2 = y2+dr*l3;
	getRhoEOSPToR(iniRho2, iniP2, 2);
	y3 = iniRho2;
	
	double k4 = der1(r+dr, y3+dr*m3);
	double l4 = der2(r+dr, y1+dr*k3, y3+dr*m3, y4);
	double m4 = der3();
	
	iniP2 = y2+dr/6.0*(l1+2.0*l2+2.0*l3+l4);
	getRhoEOSPToR(iniRho2, iniP2, 2);
	y3 = iniRho2;
	
    sol[0] = y1+dr/6.0*(k1+2.0*k2+2.0*k3+k4);
    sol[1] = y2+dr/6.0*(l1+2.0*l2+2.0*l3+l4);
    sol[2] = y3+dr/6.0*(m1+2.0*m2+2.0*m3+m4);
    sol[3] = y4+4.0*PI*r*r*y5*dr;
	
	return;
}

void getRhoBoson2F(){
	
	//Step Number
	int stepNum = 500000;
	
	///////////////
	
	//Init
	double *x = new double [stepNum]; //x coordinate
	double *f = new double [stepNum]; //f1
	double *g = new double [stepNum]; //df1/dx
	double *h = new double [stepNum]; //d2f1/dx2
	double *i = new double [stepNum]; //d3f1/dx3
	double *r = new double [stepNum]; //r coordinate
	
	double *f2 = new double [stepNum]; //f2
	double *temp1 = new double [stepNum]; //M2
	double *temp2 = new double [stepNum]; //P2
	double *temp3 = new double [stepNum]; //rho2
	double *temp4 = new double [stepNum]; //M1
	double *temp5 = new double [stepNum]; //rho1
	
	double *tempRho = new double [stepNum]; //useless now
	double *shell = new double[stepNum]; //shell
	double sol[4];
	double sol2[4];
	
	double dxrk = (30.0-1e-15)/(stepNum-1.0); //dx
	for(int i=0;i<stepNum;i++){
		x[i] = dxrk*i + 1e-15;
	}
	
	//Part 1: Solve for f2=0
	double start = -15.0;
	double end = -0.0001;
	
	while(abs(start-end) > 1e-15){
		
		double A[7];
		double dA = (end-start)/6.0;
		for(int i=0;i<7;i++){
			A[i] = dA*i + start;
		}
		
		double maxFarPoint = -1.0;
		int maxFarPointAA = -1;
		
		for(int aa=0;aa<7;aa++){
			f[0] = 1.0;
			h[0] = A[aa];
			g[0] = 0.0;
			i[0] = 0.0;
			
			for(int qq=0;qq<stepNum-1;qq++){
				loop(sol,x[qq],f[qq],g[qq],h[qq],i[qq], dxrk);
				f[qq+1] = sol[0];
				h[qq+1] = sol[1];
				g[qq+1] = sol[2];
				i[qq+1] = sol[3];
			}
			
			int farPoint = -1;
			double minValue = 1.01;
			for(int qq=0;qq<stepNum;qq++){
				if(f[qq] < minValue && f[qq] > 0){
					farPoint = qq;
					minValue = f[qq];
				}else{
					break;
				}
			}
			
			if(farPoint > maxFarPoint){
				maxFarPoint = farPoint;
				maxFarPointAA = aa;
			}
		}
		
		start = A[max(maxFarPointAA-1, 0)];
		end = A[min(maxFarPointAA+1, 6)];
		
	}
	
	//Best profile
	f[0] = 1.0;
	h[0] = start;
	g[0] = 0.0;
	i[0] = 0.0;
	
	for(int qq=0;qq<stepNum-1;qq++){
		loop(sol,x[qq],f[qq],g[qq],h[qq],i[qq], dxrk);
		f[qq+1] = sol[0];
		h[qq+1] = sol[1];
		g[qq+1] = sol[2];
		i[qq+1] = sol[3];
	}
	
	int farPoint = -1;
	double minValue = 1.01;
	for(int qq=0;qq<stepNum;qq++){
		if(f[qq] < minValue && f[qq] > 0){
			farPoint = qq;
			minValue = f[qq];
		}else{
			break;
		}
	}
	
	for(int qq=0;qq<stepNum;qq++){
		shell[qq] = f[qq]*x[qq]*x[qq];
	}
	
	double volume = shell[0] + shell[farPoint];
	for(int i=1;i<farPoint;i++){
		if(i%3 == 0){
			volume += 2.0*shell[i];
		}else{
			volume += 3.0*shell[i];
		}
	}
	volume = 3.0*dxrk*volume/8.0;
	
	///////////////
	
	double n0 = pow(1.0/volume, 4.0);
	double chi = mu/sqrt(n0);
	double mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
	double b = hBar*hBar/2.0/mStar/mBoson/mBoson;
	double rho0 = mStar*n0/4.0/PI/b/b/b;
	for(int i=0;i<stepNum;i++){
		r[i] = x[i] * volume * b;
		if(i > farPoint){
			tempRho[i] = 1e-6 * rho0;
		}else{
			tempRho[i] = f[i] * rho0;
		}
	}
	double drrk = r[1] - r[0];
	
	cout << "Part1: " << x[farPoint] << "\t" << mStar << "\t" << n0 << endl;
	
	//Part 2: Solve for f2 != 0
	
	//init
	double error = 999.9;
	
	while(abs(error) > 1e-4){
		
		start = -15.0;
		end = -0.0001;
		double iniF2 = rho2C*4.0*PI*b*b*b/mStar/n0;
		
		double p2C, p2A;
		getRhoEOSRToP(p2C, rho2C, 2);
		getRhoEOSRToP(p2A, rho2A, 2);
		
		while(abs(start-end) > 1e-15){
			
			double A[7];
			double dA = (end-start)/6.0;
			for(int i=0;i<7;i++){
				A[i] = dA*i + start;
			}
			
			double maxFarPoint = -1.0;
			int maxFarPointAA = -1;
			
			for(int aa=0;aa<7;aa++){
				
				temp1[0] = 0.0;
				temp2[0] = p2C;
				temp3[0] = rho2C;
				temp4[0] = 0.0;
				temp5[0] = rho0;
				
				f[0] = 1.0;
				f2[0] = iniF2;
				h[0] = A[aa];
				g[0] = 0.0;
				i[0] = 0.0;
				
				for(int qq=0;qq<stepNum-1;qq++){
					
					loop2(sol,x[qq],f[qq],f2[qq],g[qq],h[qq],i[qq],dxrk);
					loop3(sol2,r[qq],temp1[qq],temp2[qq],temp3[qq],temp4[qq],temp5[qq],drrk);
					f[qq+1] = sol[0];
					h[qq+1] = sol[1];
					g[qq+1] = sol[2];
					i[qq+1] = sol[3];
					
					temp1[qq+1] = sol2[0];
					temp2[qq+1] = sol2[1];
					temp3[qq+1] = sol2[2];
					temp4[qq+1] = sol2[3];
					
					if(f[qq+1] >= f[qq] || f[qq+1] < 0.0){
						f[qq+1] = f[qq];
					}
					
					if(std::isnan(temp3[qq+1]) || temp3[qq+1] < rho2A || temp3[qq+1] >= temp3[qq]){
						temp3[qq+1] = rho2A;
						temp2[qq+1] = p2A;
					}
					
					f2[qq+1] = temp3[qq+1]*4.0*PI*b*b*b/mStar/n0;
					temp5[qq+1] = f[qq+1]*rho0;
				}
				
				int farPoint = -1;
				double minValue = 1.01;
				for(int qq=0;qq<stepNum;qq++){
					if(f[qq] < minValue && f[qq] > 0){
						farPoint = qq;
						minValue = f[qq];
					}else{
						break;
					}
				}
				
				if(farPoint > maxFarPoint){
					maxFarPoint = farPoint;
					maxFarPointAA = aa;
				}
				
			}
			
			start = A[max(maxFarPointAA-1, 0)];
			end = A[min(maxFarPointAA+1, 6)];
			
		}
		
		//Best profile
		temp1[0] = 0.0;
		temp2[0] = p2C;
		temp3[0] = rho2C;
		temp4[0] = 0.0;
		temp5[0] = rho0;
		
		f[0] = 1.0;
		f2[0] = iniF2;
		h[0] = start;
		g[0] = 0.0;
		i[0] = 0.0;
		
		for(int qq=0;qq<stepNum-1;qq++){
			loop2(sol,x[qq],f[qq],f2[qq],g[qq],h[qq],i[qq],dxrk);
			loop3(sol2,r[qq],temp1[qq],temp2[qq],temp3[qq],temp4[qq],temp5[qq],drrk);
			f[qq+1] = sol[0];
			h[qq+1] = sol[1];
			g[qq+1] = sol[2];
			i[qq+1] = sol[3];
			
			temp1[qq+1] = sol2[0];
			temp2[qq+1] = sol2[1];
			temp3[qq+1] = sol2[2];
			temp4[qq+1] = sol2[3];
			
			if(f[qq+1] >= f[qq] || f[qq+1] < 0.0){
				f[qq+1] = f[qq];
			}
			
			if(std::isnan(temp3[qq+1]) || temp3[qq+1] < rho2A || temp3[qq+1] >= temp3[qq]){
				temp3[qq+1] = rho2A;
			}
			
			f2[qq+1] = temp3[qq+1]*4.0*PI*b*b*b/mStar/n0;
			temp5[qq+1] = f[qq+1]*rho0;
			if(f[qq+1] < 1e-6 || f[qq+1] > 1.0 || std::isnan(f[qq+1])){
				temp5[qq+1] = 1e-6*rho0;
			}
		}
		
		int farPoint = -1;
		double minValue = 1.01;
		for(int qq=0;qq<stepNum;qq++){
			if(f[qq] < minValue && f[qq] > 0){
				farPoint = qq;
				minValue = f[qq];
			}else{
				break;
			}
		}
		
		for(int qq=0;qq<stepNum;qq++){
			shell[qq] = f[qq]*x[qq]*x[qq];
		}
		
		double volume = shell[0] + shell[farPoint];
		for(int i=1;i<farPoint;i++){
			if(i%3 == 0){
				volume += 2.0*shell[i];
			}else{
				volume += 3.0*shell[i];
			}
		}
		volume = 3.0*dxrk*volume/8.0;
		
		double newN0 = pow(1.0/volume, 4.0);
		error = (newN0 - n0)/n0;
		n0 = newN0;
		chi = mu/sqrt(n0);
		mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
		b = hBar*hBar/2.0/mStar/mBoson/mBoson;
		rho0 = mStar*n0/4.0/PI/b/b/b;
		
		for(int i=0;i<stepNum;i++){
			r[i] = x[i] * volume * b;
			if(i > farPoint){
				tempRho[i] = 1e-6 * rho0;
			}else{
				tempRho[i] = f[i] * rho0;
			}
		}
		drrk = r[1] - r[0];
		
		rho1C = rho0;
		rho1A = rho1C * 1e-6;
		double scale = b/pow(n0, 1.0/4.0);
		
		cout << "Part2: " << x[farPoint] << "\t" << mStar << "\t" << scale << "\t" << error << endl;
		
	}
	
	///////////////
	
	double **tt = new double *[stepNum];
	for(int qq=0;qq<stepNum;qq++){
		tt[qq] = new double[3];
	}
	
	for(int qq=0;qq<stepNum;qq++){
		tt[qq][0] = r[qq]; //r
		tt[qq][1] = tempRho[qq]; //Boson
		tt[qq][2] = temp3[qq]; //NM
	}
	
	///////////////
	
	//Deallocate
	
	delete [] x;
	delete [] f;
	delete [] g;
	delete [] h;
	delete [] i;
	delete [] shell;
	delete [] r;
	delete [] tempRho;
	
	delete [] f2;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] temp4;
	delete [] temp5;
	
	///////////////
	
	for(int i=0;i<lengthStep;i++){
		double xPos = (i+0.5)*dx;
		for(int j=0;j<stepNum;j++){
			if(xPos < tt[j][0]){
				double x1 = tt[j-2][0];
				double x2 = tt[j-1][0];
				double x3 = tt[j][0];
				double x4 = tt[j+1][0];
				double y1 = tt[j-2][1];
				double y2 = tt[j-1][1];
				double y3 = tt[j][1];
				double y4 = tt[j+1][1];
				double z1 = tt[j-2][2];
				double z2 = tt[j-1][2];
				double z3 = tt[j][2];
				double z4 = tt[j+1][2];
				
				(*rho1)[i] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
						 
				(*rho2)[i] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * z1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * z2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * z3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * z4;
				
				if((*rho1)[i] < rho1A || std::isnan((*rho1)[i])){
					(*rho1)[i] = rho1A;
				}
				
				if((*rho2)[i] < rho2A || std::isnan((*rho2)[i])){
					(*rho2)[i] = rho2A;
				}
				
				break;
			}
		}
	}
	
	//Hot fix
	rhoMinDM = 1.1 * rho1A;
	rhoMinNM = 1.1 * rho2A;
	
	//Hot fix 2
	for(int i=0;i<lengthStep-1;i++){
		if((*rho2)[i+1] >= (*rho2)[i] || (*rho2)[i+1] < 0.0){
			(*rho2)[i+1] = rho2A;
		}
		
		if((*rho1)[i+1] >= (*rho1)[i] || (*rho1)[i+1] < 0.0){
			(*rho1)[i+1] = rho1A;
		}
	}
	
	boundary1D(rho1, false);
	boundary1D(rho2, false);
	
	for(int i=0;i<stepNum;i++){
		delete [] tt[i];
	}
	delete [] tt;
	
	return;
}

