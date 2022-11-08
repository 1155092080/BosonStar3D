/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the 3rd order rungekutta time evolution of the hydro - equation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

#include "rungeKutta.h"
#include "../Library/parallel.h"
#include "../Class/phyQua.h"
#include "../parameter.h"
#include "../Variable/variable.h"
#include "./fromxTox.h"
#include "./spatial.h"
#include "./boundary.h"
#include "./checkRho.h"
#include "./update.h"
#include "./findQuantumPotential.h"

void rungeKutta(int n){
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if(movingGridFlag){
		//do sth
	}
	
	spatial(uOld);
	
	for(int i=0;i<noOfEq;i++){
		for(int j=0;j<lengthStep;j++){
			(*uOne[i])[j] = (*uOld[i])[j] + (1.0/7.0)*dt*(*l[i])[j];
		}
	}
	
	if(movingGridFlag){
		//Skip for now
	}
	
	boundary2D(uOne);
	fromUToRVE(uOne);
	
	if(testModel == 0){
		
		checkRho();
		
		//spongeFlag(Skip)
		
		fromRVEToU(uOne);
		
		//xisotranFlag(Skip)
	}
	
	update(uOne, uTemp, 0, 0);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	spatial(uOne);
	
	for(int i=0;i<noOfEq;i++){
		for(int j=0;j<lengthStep;j++){
			(*uTwo[i])[j] = (*uOld[i])[j] + (3.0/1.6e1)*dt*(*l[i])[j];
		}
	}
	
	if(movingGridFlag){
		//Skip for now
	}
	
	boundary2D(uTwo);
	fromUToRVE(uTwo);
	
	if(testModel == 0){
		checkRho();
		
		//spongeFlag(Skip)
		
		fromRVEToU(uTwo);
		
		//xisotranFlag(Skip)
	}
	
	update(uTwo, uTemp, 0, 0);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	spatial(uTwo);
	
	for(int i=0;i<noOfEq;i++){
		for(int j=0;j<lengthStep;j++){
			(*uThree[i])[j] = (*uOld[i])[j] + (1.0/3.0)*dt*(*l[i])[j];
		}
	}
	
	if(movingGridFlag){
		//do sth
	}
	
	boundary2D(uThree);
	fromUToRVE(uThree);
	
	if(testModel == 0){
		checkRho();
		
		//spongeFlag(Skip)
		
		fromRVEToU(uThree);
		
		//xisotranFlag(Skip)
	}
	
	update(uThree, uTemp, 0, 0);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	spatial(uThree);
	
	for(int i=0;i<noOfEq;i++){
		for(int j=0;j<lengthStep;j++){
			(*uFour[i])[j] = (*uOld[i])[j] + (2.0/3.0)*dt*(*l[i])[j];
		}
	}
	
	if(movingGridFlag){
		//do sth
	}
	
	boundary2D(uFour);
	fromUToRVE(uFour);
	
	if(testModel == 0){
		checkRho();
		
		//spongeFlag(Skip)
		
		fromRVEToU(uFour);
		
		//xisotranFlag(Skip)
	}
	
	update(uFour, uTemp, 0, 0);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	spatial(uFour);
	
	for(int i=0;i<noOfEq;i++){
		for(int j=0;j<lengthStep;j++){
			(*uNew[i])[j] = -(3.0/4.0)*(*uOld[i])[j] + (7.0/4.0)*(*uOne[i])[j]+
						 (3.0/4.0)*dt*(*l[i])[j];
		}
	}
	
	if(movingGridFlag){
		//do sth
	}
	
	boundary2D(uNew);
	fromUToRVE(uNew);
	
	if(testModel == 0){
		//xisotranFlag
		
		//helmeosFlag
		
		//levelSetFlag
		
		//flameFlag
	}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if(testModel == 0){
		//flameFlag (Skip for now)
	}
	
	////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if(testModel == 0){
		//movingGridFlag
		
		checkRho();
		
		//spongeFlag
		
		fromRVEToU(uNew);
	}
	
	/*if(qpDMFlag){
		findQuantumPotential();
	}*/
	
	//if(n%2 == 0){
	update(uNew, uOld, 1, 1);
	//}else{
	//	update(uNew, uOld, 0, 1);
	//}
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	return;
}

