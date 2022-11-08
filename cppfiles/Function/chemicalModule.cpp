#include "chemicalModule.h"
#include "../Variable/variable.h"
#include <vector>
using namespace std;

/*
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                             
 ! This subroutine is obtained from Timmes' nuclear network !
 ! program to calculate the abar and zbar for a given       !
 ! chemical composition in a flexible structure             !
 ! Merged by Leung Shing Chi in 2016                        !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */
 
double sum(vector<double> &a){
	double tmp = 0.0;
	for(int i=0;i<a.size();i++){
		tmp += a[i];
	}
	return tmp;
}

void azBar(double *xMass, double &aBar, double &zBar){
	
	vector<double> yMass(totalIon, 0.0);
	
	for(int i=0;i<totalIon;i++){
		yMass[i] = xMass[i]/wion[i];
	}
	
	double wBar = 1.0/sum(yMass);
	
	vector<double> tmp(totalIon, 0.0);
	for(int i=0;i<totalIon;i++){
		tmp[i] = aion[i]*yMass[i];
	}
	double sum1 = sum(tmp);
	
	aBar = wBar * sum1;
	
	for(int i=0;i<totalIon;i++){
		tmp[i] = zion[i]*yMass[i];
	}
	
	double ye = sum(tmp);
	zBar = wBar*ye;
	
	return;
}

