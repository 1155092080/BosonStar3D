#ifndef __VARIABLE_H__
#define __VARIABLE_H__

#include "../Class/phyQua.h"
#include <string> 
#include <vector>
#include <complex>
#include <fftw3.h>

void buildVariable();
void destroyVariable();

//Previously initial
extern double k2;
extern double gamma2;
extern double k1;
extern double gamma1;

extern double **eosTableList;
extern double **eosTableListBoson;

//Previously definition
extern double globalTime;
extern double dx;
extern double dt;

//Previously buildHydro
extern PhyQua *rho;
extern PhyQua *phi;
extern PhyQua *phip;

extern PhyQua *rho2;
extern PhyQua *p2;
extern PhyQua *epsilon2;
extern PhyQua *dpdrho2;
extern PhyQua *dpdepsilon2;
extern PhyQua *vel2;
extern PhyQua *temp2;

extern PhyQua *rho1;
extern PhyQua *p1;
extern PhyQua *epsilon1;
extern PhyQua *dpdrho1;
extern PhyQua *dpdepsilon1;
extern PhyQua *vel1;
extern PhyQua *quantumPotential;
extern PhyQua *quantumPotential1;

//Previously buildWENO
extern PhyQua **uOld;
extern PhyQua **uNew;
extern PhyQua **uOne;
extern PhyQua **uTwo;
extern PhyQua **uThree;
extern PhyQua **uFour;
extern PhyQua **uTemp;
extern PhyQua **l;
extern PhyQua **f;
extern PhyQua **sa;
extern PhyQua **sb;
extern PhyQua **sc;
extern PhyQua **flux;
extern PhyQua **dfdx;
extern double rhoMinDM;
extern double rhoMinNM;

extern int *bFac;

extern int noOfEq;

extern int imin1;
extern int imax1;
extern int irho1;
extern int ivel1;
extern int itau1; 

extern int imin2;
extern int imax2;
extern int irho2;
extern int ivel2;
extern int itau2;

//Previously chemicalModule
extern PhyQua *aBar2;
extern PhyQua *zBar2;

extern PhyQua **xIso;

extern int totalIon;
extern std::vector<std::string> ioName;
extern std::vector<double> aion;
extern std::vector<double> zion;
extern std::vector<double> bion;
extern std::vector<double> nion;
extern std::vector<double> mion;
extern std::vector<double> wion;

extern int cc12;
extern int cni56;

extern double aBar2Ini, zBar2Ini;

//Previously getEpsilon
extern double epsilon1A;
extern double epsilon2A;

//Previously spatial
extern double *alpha;
extern double rhoMinDM;
extern double rhoMinNM;

//GP equation
extern std::complex<double> *psiSave;
extern std::complex<double> *psiStar;
extern std::complex<double> *psiStarFourier;
extern fftw_complex *psiStarFourierFFTW;
extern double *expo;
extern double *phi3D;

extern double Ra;
extern double epsilon;
extern double delta;
extern double M;
extern double omega;

#endif
