#include "wenoModule.h"
#include "../parameter.h"
#include "../Variable/variable.h"

void constc(double c[][3]){
	c[0][0] = 1.1e1/6.0;
	c[0][1] = -7.0/6.0;
	c[0][2] = 1.0/3.0;
	
	c[1][0] = 1.0/3.0;
	c[1][1] = 5.0/6.0;
	c[1][2] = -1.0/6.0;
	
	c[2][0] = -1.0/6.0;
	c[2][1] = 5.0/6.0;
	c[2][2] = 1.0/3.0;
	
	c[3][0] = 1.0/3.0;
	c[3][1] = -7.0/6.0;
	c[3][2] = 1.1e1/6.0;
	
	return;
}

void constd(double d[]){
	d[0] = 3.0/1.0e1;
	d[1] = 3.0/5.0;
	d[2] = 1.0/1.0e1;
	
	return;
}

void weno(PhyQua *v, PhyQua *vPlus, PhyQua *vMinus){
	
	double c[4][3];
	double d[3];
	double td[3];
	double vrhs[3];
	double vlhs[3];
	double alphaWENO[3];
	double beta[3];
	double omega[3];
	double talpha[3];
	double tomega[3];
	double tmp;
	
	constc(c);
	constd(d);
	
	for(int i=0;i<=2;i++){
		td[i] = d[2-i];
	}
	
	for(int i=-3;i<lengthStep+3;i++){
		for(int r=0;r<=2;r++){
			vrhs[r] = 0.0;
			for(int j=0;j<=2;j++){
				vrhs[r] += c[r+1][j] * (*v)[i-r+j];
			}
		}
		
		for(int r=0;r<=2;r++){
			vlhs[r] = 0.0;
			for(int j=0;j<=2;j++){
				vlhs[r] += c[r][j] * (*v)[i-r+j];
			}
		}
		
		beta[0] = (1.3e1/1.2e1)*((*v)[i]-2.0*(*v)[i+1]+(*v)[i+2])*((*v)[i]-2.0*(*v)[i+1]+(*v)[i+2])+
				  (1.0/4.0)*(3.0*(*v)[i]-4.0*(*v)[i+1]+(*v)[i+2])*(3.0*(*v)[i]-4.0*(*v)[i+1]+(*v)[i+2]);
		
		beta[1] = (1.3e1/1.2e1)*((*v)[i-1]-2.0*(*v)[i]+(*v)[i+1])*((*v)[i-1]-2.0*(*v)[i]+(*v)[i+1])+
				  (1.0/4.0)*((*v)[i-1]-(*v)[i+1])*((*v)[i-1]-(*v)[i+1]);
		
		beta[2] = (1.3e1/1.2e1)*((*v)[i-2]-2.0*(*v)[i-1]+(*v)[i])*((*v)[i-2]-2.0*(*v)[i-1]+(*v)[i])+
				  (1.0/4.0)*((*v)[i-2]-4.0*(*v)[i-1]+3.0*(*v)[i])*((*v)[i-2]-4.0*(*v)[i-1]+3.0*(*v)[i]);
		
		for(int r=0;r<=2;r++){
			alphaWENO[r] = d[r]/((beta[r]+smallpara)*(beta[r]+smallpara));
		}
		
		tmp = 0.0;
		
		for(int r=0;r<=2;r++){
			tmp += alphaWENO[r];
		}
		
		for(int r=0;r<=2;r++){
			omega[r] = alphaWENO[r]/tmp;
		}
		
		for(int r=0;r<=2;r++){
			talpha[r] = td[r]/((beta[r]+smallpara)*(beta[r]+smallpara));
		}
		
		tmp = 0.0;
		
		for(int r=0;r<=2;r++){
			tmp += talpha[r];
		}
		
		for(int r=0;r<=2;r++){
			tomega[r] = talpha[r]/tmp;
		}
		
		(*vMinus)[i] = 0.0;
		
		for(int r=0;r<=2;r++){
			(*vMinus)[i] += omega[r]*vrhs[r];
		}
		
		(*vPlus)[i] = 0.0;
		
		for(int r=0;r<=2;r++){
			(*vPlus)[i] += tomega[r]*vlhs[r];
		}
		
	}
	
	return;
}

