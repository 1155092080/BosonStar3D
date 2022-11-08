#ifndef __PARAMETER_H__
#define __PARAMETER_H__

//dx, cfl
extern double dxOrg;
extern double cfl;
extern int totalTimeStep;
extern double totalTime;

//Flag
extern bool wGravityI;
extern bool dmFlag;
extern bool runDMFlag;
extern bool levelSetFlagDM;
extern bool levelSetFlagNM;
extern bool sedovFlag;
extern bool helmeosFlag;
extern bool movingGridFlag;
extern bool flameFlag;
extern bool xisotranFlag;
extern bool outputPotentialFlag;

extern double gs2;
extern double ye2;
extern double me2;
extern double mb2;

extern double ye1;
extern double me1;
extern double mb1;

extern double totalLength;
extern int lengthStep;

extern int testModel;

extern int iniAcc;

extern double rho1C;
extern double rho1A;
extern double rho2C;
extern double rho2A;

extern double xc12Ini;
extern double xni56Ini;

extern double temp2Ini;
extern double temp2A; 

extern double iniVel1;
extern double iniVel2;

extern double vel1A;
extern double vel2A;

extern int boundaryFlag[2];

extern int spDimI;

extern double smallpara;

extern int relaxMax;		
extern double tolerance;

extern int eosTableWidth;

extern int printEvery;
extern int outputEvery;

extern double scattering;
extern double mBoson;
extern bool bosonFlag;
extern bool relEOS; 
extern bool fermiEOS;
extern bool bosonEOS;
extern bool qpFlag;

extern double mu;
extern double psi;
extern bool bosonDMFlag;
extern bool bosonDMEOS;
extern bool qpDMFlag;

extern bool aprEOS;
extern bool newAPREOS;
extern bool TOVFlag;

#endif
