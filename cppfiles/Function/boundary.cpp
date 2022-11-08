#include <iostream>
#include "boundary.h"
#include "../Class/phyQua.h"
#include "../Variable/variable.h"
#include "../parameter.h"

void boundary1D(PhyQua *v, bool sign){
	
	int fac;
	
	if(!sign){
		fac = 1;
	}else{
		fac = -1;
	}
	
	//Inner boundary
	if(boundaryFlag[0] == 0){
		for(int i=0;i<5;i++){
			(*v)[-1-i] = (*v)[lengthStep-1-i];
		}
	}else if(boundaryFlag[0] == 1){
		for(int i=0;i<5;i++){
			(*v)[-1-i] = fac * (*v)[i];
		}
	}else{
		for(int i=0;i<5;i++){
			(*v)[-1-i] = (*v)[0];
		}
	}
	
	//Outer boundary
	if(boundaryFlag[1] == 0){
		for(int i=0;i<5;i++){
			(*v)[lengthStep+i] = (*v)[i];
		}
	}else if(boundaryFlag[1] == 1){
		for(int i=0;i<5;i++){
			(*v)[lengthStep+i] = fac * (*v)[lengthStep-1-i];
		}
	}else{
		for(int i=0;i<5;i++){
			(*v)[lengthStep+i] = (*v)[lengthStep-1];
		}
	}
	
	return;
}

void boundary2DX(PhyQua **x){
	
	//Inner boundary
	if(boundaryFlag[0] == 0){
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[-1-j] = (*x[i])[lengthStep-1-j];
			}
		}
	}else if(boundaryFlag[0] == 1){
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[-1-j] = (*x[i])[j];
			}
		}
	}else{
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[-1-j] = (*x[i])[0];
			}
		}
	}
	
	//Outer boundary
	if(boundaryFlag[1] == 0){
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[lengthStep+j] = (*x[i])[j];
			}
		}
	}else if(boundaryFlag[1] == 1){
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[lengthStep+j] = (*x[i])[lengthStep-1-j];
			}
		}
	}else{
		for(int i=0;i<totalIon;i++){
			for(int j=0;j<5;j++){
				(*x[i])[lengthStep+j] = (*x[i])[lengthStep-1];
			}
		}
	}
	
	return;
}

void boundary2D(PhyQua **u){
	
	//Inner boundary
	if(boundaryFlag[0] == 0){
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[-1-j] = (*u[i])[lengthStep-1-j];
			}
		}
	}else if(boundaryFlag[0] == 1){
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[-1-j] = bFac[i]*(*u[i])[j];
			}
		}
	}else{
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[-1-j] = (*u[i])[0];
			}
		}
	}
	
	//Outer boundary
	if(boundaryFlag[1] == 0){
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[lengthStep+j] = (*u[i])[j];
			}
		}
	}else if(boundaryFlag[1] == 1){
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[lengthStep+j] = bFac[i]*(*u[i])[lengthStep-1-j];
			}
		}
	}else{
		for(int i=0;i<=imax2;i++){
			for(int j=0;j<5;j++){
				(*u[i])[lengthStep+j] = (*u[i])[lengthStep-1];
			}
		}
	}
	
	return;
	
}

