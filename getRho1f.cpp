#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "../parameter.h"
using namespace std;

//double sum(vecotr<double>);
//void azBar(double *xMass, double &aBar, double &zBar);


//The number of hydrostatic equation and the improved accuracy
int noOfEqIni;
int more = pow(10, iniAcc);
int totalIon = 1;
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
	
	// xMass is the composition of each isotope, adding up together is 1
    // wion is the atomic weight of each isotope, wchich can be approximated by atomic number
    // aion is number of nucleons
    // z is number of protons (electrons)
    vector<double> yMass(totalIon, 0.0); // mean atomic mass
	
	for(int i=0;i<totalIon;i++){
		yMass[i] = xMass[i]/wion[i];
	}
	
	double wBar = 1.0/sum(yMass);
	
	vector<double> tmp(totalIon, 0.0);
	for(int i=0;i<totalIon;i++){
		tmp[i] = aion[i]*yMass[i];
	}
	double sum1 = sum(tmp); // mean fraction of nucleons
	
	aBar = wBar * sum1; // mean number of nucleons
	
	for(int i=0;i<totalIon;i++){
		tmp[i] = zion[i]*yMass[i]; // mean fraction of protons
	}
	
	double ye = sum(tmp);
	zBar = wBar*ye; // mean number of protons/electrons
	
	return;
}


void iniDer1F(double *der, double x, double *y){
	der[0] = 4.0*PI*x*x*y[2]; //dm = 4pi r^2 rho
	if(TOVFlag){
		if(y[0] == 0.0){
			der[1] = 0.0;
		}else{
			der[1] = -(y[0]*y[2])/(x*x)*(1.0+y[1]/y[2])*(1.0+4.0*PI*x*x*x*y[1]/y[0])/(1.0-2.0*y[0]/x);
		}
	}else{
		der[1] = -(y[0]*y[2])/(x*x);  // dm = -G*m*rho/r^2
	}
	der[2] = 0.0;
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
    // pressure at center and atmonsphere
    double p2C, p2A;
	
	//We assign initial chemical composition accordingly
    int cc12 = 0;// label type of isotope, here we define xc12Ini=1 and xni56Ini=0
    int cni56 = 1;
	xIsoIni[cc12] = xc12Ini;
	xIsoIni[cni56] = xni56Ini;

    //Now convert the composition into mean atomic and mass number
	azBar(xIsoIni, aBar2Ini, zBar2Ini);
	
	//At the atmosphere, there should not be any composition
	double xISOA = 0.0;

    //We choose which table is needed to read according to the chosen EOS
	if(helmeosFlag){
		//do sth
	}else{
		getRhoEOSRToP(p2C, rho2C, 2);
		getRhoEOSRToP(p2A, rho2A, 2);
	}

    //We assign the initial value of y at center, in our definition, y[0], y[1], y[2] are mass, p and rho
	y[0][0] = 0.0;
	y[1][0] = p2C;
	y[2][0] = rho2C;

    //this is the main structure of RK-4 method
	
	for(int j=0;j<lengthMoreStep-1;j++){
		
		double iniP2, iniRho2;
        //y0 = y(x0)
        //x0 = x0
        //k1 = dx*f'(x0, y0)
		for(int i=0;i<noOfEqIni;i++){
			yZero[i] = y[i][j];
		}

        double x = j*dxMore;
		
		iniDer1F(der, x, yZero); // calculate dm/dr, dp/dr and drho/dr. But rho is obtained by EOS instead of drho/dr
		
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

        //y1 = y0 + k1/2
        //x1 = x0 + dx/2
        //k2 = dx*f'(x1, y1)
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

        //y2 = y0 + k2/2
        //x2 = x0 + dx/2
        //k3 = dx*f'(x2, y2)

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

        //y3 = y0 + k3
        //x3 = x0 + dx
        //k3 = dx*f'(x3, y3)
        // y_new = y_old + k1/6 + k2/3 + k3/3 + k4/6
        x = j*dxMore+dxMore;
		
		iniDer1F(der, x, yThree);
		
		double k4[3];
		for(int i=0;i<noOfEqIni;i++){
			k4[i] = dxMore*der[i];
			y[i][j+1] = y[i][j] + k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0;
		}

        iniP2 = y[1][j+1];
        // test if updated pressure is still more than atmonsphere pressure, if not, assign it atmonspheric pressure
		if(iniP2 <= p2A){
			iniP2 = p2A;
			y[1][j+1] = p2A;
		}
        getRhoEOSPToR(iniRho2, iniP2, 2); // Calculate corresponding rho

		y[2][j+1] = iniRho2;
    }
    
    //Question: Why j+1/2 more ?
    for(int j=0;j<lengthStep;j++){
		(*rho2)[j] = y[2][(j+1)*more - more/2];
		
		if((*rho2)[j] < rho2A || std::isnan((*rho2)[j])){
			(*rho2)[j] = rho2A;
		}
	}
}
