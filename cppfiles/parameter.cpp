/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! User-defined parameters	   					   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include "parameter.h"

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is related to the basic parameters governing the simulation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

//! Boundary condition for the simulation box !
//! The first one is inner boundary and second is outer !
//! 0 = periodic, 1 = reflecting, 2 = outgoing !
int boundaryFlag[2] = {1,2};

//Spacial structure: 0 = 1-D Cartesian; 1 = cylindrical symmetric; 2 = spherical symmetric
int spDimI = 2;

//Physical length (dimensionless) of the simulation box
double totalLength = 10000.0;

//Value of spatial grid size dx
double dxOrg = 40.0;

//Value of CFL number that govern the stability condition (dt / dx)
double cfl = 0.1;

//Physical time (dimensionless) that stop the program
double totalTime = 2e5;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Test model														 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
*/

//Hydro Test Model. Choose from 1, 2, 3 ,4, 5, 6. 0 Means not running hydro test
int testModel = 0;

//Sedov model
bool sedovFlag = false;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Output setting													 ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
*/

//Print to console every
int printEvery = 1;

//Output to file every
int outputEvery = 10000;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is some logical flag to determine whether some extra feature should be turned on or not !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

//Presence of gravity
bool wGravityI = true;

//The presence of dark matter component
bool dmFlag = false;

//Whether dark matter component are movable
bool runDMFlag = false;

//! 1 = use Fermi Gas EOS For NM !
bool fermiEOS = true;

////Use APR for NM
bool aprEOS = false;
bool newAPREOS = false;

//! 0 = use newtonian polytropic EOS For NM !
bool relEOS = true;

//...(Skip it please)
bool levelSetFlagDM = false;
bool levelSetFlagNM = false;
bool helmeosFlag = false;
bool movingGridFlag = false;
bool outputPotentialFlag = false;
bool flameFlag = false;
bool xisotranFlag = false;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section governs the initial condition of the hydrostatic star ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
*/

//Initial central density of dark matter component, do not set to zero
double rho1C = 1.0;//0.005 * 1.619e-9;
//double rho1C = 0.0;

//Initial central density of normal matter component, do not set to zero
double rho2C = 2*1e10 * 1.619e-18;
//double rho2C = 0.0046 * 1.619e-9;

//Density of the dark matter at the atmosphere, do not set to zero
double rho1A = rho1C * 1.0e-6;

//Density of the normal matter at the atmosphere, do not set to zero
double rho2A = rho2C * 1.0e-8;

//Initial velocity of normal matter, it can be set to zero
double iniVel2 = 0.0;

//Initial velocity of dark matter, it can be set to zero
double iniVel1 = 0.0;

//Velocity of the dark matter at the atmosphere, it can be set to zero
double vel1A = 0.0;

//Velocity of the normal matter at the atmosphere, it can be set to zero
double vel2A = 0.0;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This section governs the physics of the EOS of NM or DM    !
! CAUTION : We assumed an ideal completely degenerate fermi  !
! gas EOS. To change the EOS, you need to input the required !
! parameters by yourself. For example, temperature           !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
*/

//Fermionic mass (dark matter particles) for dark matter
double me1 = 8.965e-59;
//double me1 = 4.5771e-61;

//Fermionic mass (electrons) for normal matter
double me2 = 4.5771e-61;

//Baryonic mass for dark matter
double mb1 = 8.965e-59;
//double mb1 = 8.4158e-58;

//Baryonic mass for normal matter
double mb2 = 8.4158e-58;

//Dark matter fraction for dark matter
double ye1 = 1.0;
//double ye1 = 0.5;

//Electron fraction for normal matter
double ye2 = 0.5;

//Multiplicity factor for normal matter
double gs2 = 2.0;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section governs the chemical composition the hydrostatic star !
! The composition Xi is defined by rhoi/rho, the fraction of density !
! occupied by that elements for each desity at a radial distance     !
! Caution : The composition should sum up to 1, also, once you       !
! assume variables composition, you can only use helmeos	    	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

//Composition of carbon 12
double xc12Ini = 1.0;

//Composition of Nikle 56
double xni56Ini = 0.0;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is the extra physics related to the !
! Helmeos finite temperature EOS. If you hate it   !
! Please feel free to delete them.                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

//Inital temperature (assume isothermal) of the star
double temp2Ini = 1.0e-1;

//atmospheric temperature of the star
double temp2A = 1.0e-1;

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some unimportant parameters that should not be changed				   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
//The maximum number of iteration for Rungekutta method
int totalTimeStep = 0;

//Accuracy of the initial date if it is generated
//(Do NOT change this unless you know what you are doing)
int iniAcc = 2;

//The small parameter in WENO scheme that avoid the coefficient to be dividing by zero
//(Do NOT change this unless you know what you are doing)
double smallpara = 1.0e-80;

//! Maximum run time in relaxation of the potential !
//! (Do NOT change this unless you know what you are doing) !
int relaxMax = 100000;

//! Tolerance in relaxation of the potential !
//! (Do NOT change this unless you know what you are doing) !			
double tolerance = 3.0e-8;

int eosTableWidth = 2400;

/*
Test Boson Feature
*/

//Use TOV correction
bool TOVFlag = false;

//Use polytropic Boson EOS for NM
bool bosonFlag = false;
//Use rel. Boson EOS for NM
bool bosonEOS = false;
//Use 4-th order solver for Boson NM
bool qpFlag = false;

//Use polytropic Boson EOS for DM
bool bosonDMFlag = false;
//Use rel. Boson EOS for DM
bool bosonDMEOS = false;
//Use 4-th order solver for Boson DM
bool qpDMFlag = false;

//Properties of Boson
double scattering = 5e-074;
double mBoson = 2e-077;
double mu = 3;//0.5;
double psi;
