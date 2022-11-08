#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <fftw3.h>
#include "../parameter.h"
#include "../Class/phyQua.h"
#include "./definition.h"
using namespace std;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Previously initial.h

double k2, gamma2, k1, gamma1;

double **eosTableList = NULL;
double **eosTableListBoson = NULL;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Previously definition.h

double globalTime;
double dx;
double dt;
int lengthStep;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PhyQua *rho = NULL;
PhyQua *phi = NULL;
PhyQua *phip = NULL;

PhyQua *rho2 = NULL;
PhyQua *p2 = NULL;
PhyQua *epsilon2 = NULL;
PhyQua *dpdrho2 = NULL;
PhyQua *dpdepsilon2 = NULL;
PhyQua *vel2 = NULL;
PhyQua *temp2 = NULL;
PhyQua *quantumPotential = NULL;

PhyQua *rho1 = NULL;
PhyQua *p1 = NULL;
PhyQua *epsilon1 = NULL;
PhyQua *dpdrho1 = NULL;
PhyQua *dpdepsilon1 = NULL;
PhyQua *vel1 = NULL;
PhyQua *quantumPotential1 = NULL;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PhyQua **uOld = NULL;
PhyQua **uNew = NULL;
PhyQua **uOne = NULL;
PhyQua **uTwo = NULL;
PhyQua **uThree = NULL;
PhyQua **uFour = NULL;
PhyQua **uTemp = NULL;
PhyQua **f = NULL;
PhyQua **sa = NULL;
PhyQua **sb = NULL;
PhyQua **sc = NULL;
PhyQua **fP = NULL;
PhyQua **fM = NULL;
PhyQua **flux = NULL;
PhyQua **dfdx = NULL;
PhyQua *v = NULL;
PhyQua *vP = NULL;
PhyQua *vM = NULL;
PhyQua *fluxP = NULL;
PhyQua *fluxM = NULL;
PhyQua **l = NULL;

int *bFac = NULL;

double rhoMinDM;
double rhoMinNM;

int noOfEq;

int imin1, irho1, imax1, ivel1, itau1, imin2, imax2, irho2, ivel2, itau2;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//Number of isotope in the star
int totalIon = 2;

int cc12, cni56;

double aBar2Ini, zBar2Ini;

vector<double> wion(totalIon, 0.0);
vector<double> aion(totalIon, 0.0);
vector<double> zion(totalIon, 0.0);

PhyQua *aBar2 = NULL;
PhyQua *zBar2 = NULL;

PhyQua **xIso = NULL;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double epsilon1A;
double epsilon2A;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double *alpha = NULL;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

complex<double> *psiSave = NULL;
complex<double> *psiStar = NULL;
complex<double> *psiStarFourier = NULL;
fftw_complex *psiStarFourierFFTW = NULL;
double *expo = NULL;
double *phi3D = NULL;

double Ra;
double epsilon;
double delta;
double M;
double omega;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void buildVariable(){
	
	//Initialize globalTime
	globalTime = 0.0;
	
	//Initialize dx
	dx = dxOrg;
	
	//Initialize dt
	dt = cfl*dx;
	
	//The total number of array stored by each variables
	lengthStep = (int)(totalLength/dxOrg);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	eosTableList = new double *[16*eosTableWidth];
	for(int i=0;i<16*eosTableWidth;i++){
		eosTableList[i] = new double[2];
	}
	
	eosTableListBoson = new double *[16*eosTableWidth];
	for(int i=0;i<16*eosTableWidth;i++){
		eosTableListBoson[i] = new double[2];
	}
	
	rho = new PhyQua(lengthStep);
	phi = new PhyQua(lengthStep);
	phip = new PhyQua(lengthStep);
	
	rho2 = new PhyQua(lengthStep);
	p2 = new PhyQua(lengthStep);
	epsilon2 = new PhyQua(lengthStep);
	dpdrho2 = new PhyQua(lengthStep);
	dpdepsilon2 = new PhyQua(lengthStep);
	vel2 = new PhyQua(lengthStep);
	temp2 = new PhyQua(lengthStep);
	quantumPotential = new PhyQua(lengthStep);
	
	rho1 = new PhyQua(lengthStep);
	p1 = new PhyQua(lengthStep);
	epsilon1 = new PhyQua(lengthStep);
	dpdrho1 = new PhyQua(lengthStep);
	dpdepsilon1 = new PhyQua(lengthStep);
	quantumPotential1 = new PhyQua(lengthStep);
	
	if(runDMFlag){
		vel1 = new PhyQua(lengthStep);
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	bFac = new int[50];
	
	/*
	We now find out how many conservative equation is needed to solve
	We first do the DM part
	We do the DM part only if the DM is precense and movable
	*/
	
	if(dmFlag && runDMFlag){
		
		//DM density conservative equation is needed to be solved
		//increase the maximum number of equation needed to be solved
		imin1 = noOfEq;
		
		//assign the order of conservative equation accordingly
		irho1 = noOfEq;
		
		//boundary factor of density related to boundary condition
		bFac[noOfEq] = 1;
		
		//increase no of eq
		noOfEq++;
		
		//DM velocity conservative equation is needed to be solved
		imax1 = noOfEq;
		ivel1 = noOfEq;
		bFac[noOfEq] = -1;
		noOfEq++;
		
		//DM epsilon conservative equation is needed to be solved
		imax1 = noOfEq;
		itau1 = noOfEq;
		bFac[noOfEq] = 1;
		noOfEq++;
		
		//Now we check if we need to solve extra conservative equation
		if(levelSetFlagDM){
			//Skip for now
		}
		
	}
	
	//Now we turn on to the NM part
	
	//NM density conservative equation is needed to be solved
	imin2 = noOfEq;
	irho2 = noOfEq;
	bFac[noOfEq] = 1;
	noOfEq++;
	
	//NM velocity conservative equation is needed to be solved
	imax2 = noOfEq;
	ivel2 = noOfEq;
	bFac[noOfEq] = -1;
	noOfEq++;
	
	//NM epsilon conservative equation is needed to be solved
	imax2 = noOfEq;
	itau2 = noOfEq;
	bFac[noOfEq] = 1;
	noOfEq++;
	
	//Now we check if we need to solve extra conservative equation
	if(levelSetFlagNM){
		//Skip for now
	} 
	
	//We check whether there is a deflagration
	if(flameFlag){
		//Skip for now
	}
	
	/*
	! Set up the conservative equation for chemical composition !
	! CAUTION : If you change the number of elements !
	! Please also change this section to include the elements you want !
	*/
	
	if(xisotranFlag){
		//Skip for now
	}
	
	//Now we allocate the u, and l accordingly
	uOld = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uOld[i] = new PhyQua(lengthStep);
	}
	
	uNew = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uNew[i] = new PhyQua(lengthStep);
	}
	
	l = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		l[i] = new PhyQua(lengthStep);
	}
	
	uOne = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uOne[i] = new PhyQua(lengthStep);
	}
	
	uTwo = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uTwo[i] = new PhyQua(lengthStep);
	}
	
	uThree = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uThree[i] = new PhyQua(lengthStep);
	}
	
	uFour = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uFour[i] = new PhyQua(lengthStep);
	}
	
	uTemp = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		uTemp[i] = new PhyQua(lengthStep);
	}
	
	f = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		f[i] = new PhyQua(lengthStep);
	}
	
	sa = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		sa[i] = new PhyQua(lengthStep);
	}
	
	sb = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		sb[i] = new PhyQua(lengthStep);
	}
	
	sc = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		sc[i] = new PhyQua(lengthStep);
	}
	
	flux = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		flux[i] = new PhyQua(lengthStep);
	}
	
	dfdx = new PhyQua *[noOfEq];
	for(int i=0;i<noOfEq;i++){
		dfdx[i] = new PhyQua(lengthStep);
	}
	
	rhoMinNM = 1.1 * rho2A;
	rhoMinDM = 1.1 * rho1A;
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	//allocate mean atomic number and mean atomic mass
	aBar2 = new PhyQua(lengthStep);
	zBar2 = new PhyQua(lengthStep);
	
	vector<string> ioName(totalIon, "");
	vector<double> bion(totalIon, 0.0);
	vector<double> nion(totalIon, 0.0);
	vector<double> mion(totalIon, 0.0);
	
	//allocate mean atomic composition
	xIso = new PhyQua *[totalIon];
	for(int i=0;i<totalIon;i++){
		xIso[i] = new PhyQua(lengthStep);
	}
	
	/*
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 ! This subroutine calculate all the nuclear variables for the elements ! 
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	*/	
	
	//this is the conversion between different units of energy
	//Include MeV to erg and MeV to gram 
	const double ev2erg = 1.60217648740e-12;
	
	//speed of light
	const double clight = 2.99792458e10;
	
	//avogardo number
	const double avo = 6.0221367e23;
	
	//Of course, the conversion is fixed by physical constant
	const double mev2erg = ev2erg*1.0e6;
	
	const double mev2gr = mev2erg/(clight*clight);
	
	//set the id numbers of the elements
	//If you add more elements, this needed to be changed
	cc12 = 0;
	cni56 = 1;
	
	//set the names of the elements
	ioName[cc12] = "c12";
	ioName[cni56] = "ni56"; 
	
	//set the number of nucleons in the element
	aion[cc12] = 12.0;
	aion[cni56] = 56.0;
	
	//set the number of protons in the element
	zion[cc12] = 6.0;
	zion[cni56] = 28.0;
	
	//set the binding energy of (MeV) the element
	bion[cc12] = 92.16294;
	bion[cni56] = 484.00300;
	
	for(int i=0;i<totalIon;i++){
		//set the number of neutrons and mass
		nion[i] = aion[i] - zion[i];
		//mass in gram of each isotope
		//the calculation is total neutron mass + total proton mass
		//minus the binding energy 
		mion[i] = nion[i]*1.67492721184e-24 + zion[i]*1.67262163783e-24 - bion[i]*mev2gr;
		
		//Why do this twice??
		//molar mass
		wion[i] = avo * mion[i];
		//a common approximation
		wion[i] = aion[i];
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	psiSave = new complex<double>[2*lengthStep*2*lengthStep*2*lengthStep];
	psiStar = new complex<double>[2*lengthStep*2*lengthStep*2*lengthStep];
	psiStarFourier = new complex<double>[2*lengthStep*2*lengthStep*2*lengthStep];
	psiStarFourierFFTW = (fftw_complex*)fftw_malloc(2*lengthStep*2*lengthStep*2*lengthStep*sizeof(fftw_complex));
	expo = new double[2*lengthStep*2*lengthStep*2*lengthStep];
	phi3D = new double[2*lengthStep*2*lengthStep*2*lengthStep];

	Ra = sqrt(abs(scattering*hBar*hBar/mBoson/mBoson/mBoson))/6.77193e-6/100;
	epsilon = 1.0;
	omega = mBoson*mBoson/scattering/hBar*2.03017e5;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	return;
}

void destroyVariable(){
	
	for(int i=0;i<16*eosTableWidth;i++){
		delete eosTableList[i];
	}
	delete [] eosTableList;
	
	for(int i=0;i<16*eosTableWidth;i++){
		delete eosTableListBoson[i];
	}
	delete [] eosTableListBoson;
	
	delete rho;
	delete phi;
	delete phip;
	
	delete rho2;
	delete p2;
	delete epsilon2;
	delete dpdrho2;
	delete dpdepsilon2;
	delete vel2;
	delete temp2;
	delete quantumPotential;
	
	delete rho1;
	delete p1;
	delete epsilon1;
	delete dpdrho1;
	delete dpdepsilon1;
	delete quantumPotential1;

	delete psiSave;
	delete psiStar;
	delete psiStarFourier;
	fftw_free(psiStarFourierFFTW);
	delete expo;
	delete phi3D;
	
	if(runDMFlag){
		delete vel1;
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	for(int i=0;i<noOfEq;i++){
		delete uNew[i];
	}
	delete [] uNew;
	
	for(int i=0;i<noOfEq;i++){
		delete uOld[i];
	}
	delete [] uOld;
	
	for(int i=0;i<noOfEq;i++){
		delete l[i];
	}
	delete [] l;
	
	for(int i=0;i<noOfEq;i++){
		delete uOne[i];
	}
	delete [] uOne;
	
	for(int i=0;i<noOfEq;i++){
		delete uTwo[i];
	}
	delete [] uTwo;
	
	for(int i=0;i<noOfEq;i++){
		delete uThree[i];
	}
	delete [] uThree;
	
	for(int i=0;i<noOfEq;i++){
		delete uFour[i];
	}
	delete [] uFour;
	
	for(int i=0;i<noOfEq;i++){
		delete uTemp[i];
	}
	delete [] uTemp;
	
	
	for(int i=0;i<noOfEq;i++){
		delete f[i];
	}
	delete [] f;
	
	for(int i=0;i<noOfEq;i++){
		delete sa[i];
	}
	delete [] sa;
	
	for(int i=0;i<noOfEq;i++){
		delete sb[i];
	}
	delete [] sb;
	
	for(int i=0;i<noOfEq;i++){
		delete sc[i];
	}
	delete [] sc;
	
	for(int i=0;i<noOfEq;i++){
		delete flux[i];
	}
	delete [] flux;
	
	for(int i=0;i<noOfEq;i++){
		delete dfdx[i];
	}
	delete [] dfdx;
	
	delete [] bFac;
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	delete aBar2;
	delete zBar2;
	
	for(int i=0;i<totalIon;i++){
		delete xIso[i];
	}
	delete [] xIso;
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	return;
}
