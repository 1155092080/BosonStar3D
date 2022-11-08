#include "update.h"
#include "../parameter.h"
#include "../Variable/variable.h"
#include "./findPotential.h"
#include "./findPressure.h"
#include "./findQuantumPotential.h"

#include <iostream>
using namespace std;

void update(PhyQua **v, PhyQua **u, int p, int gate){
	
	for(int i=0;i<noOfEq;i++){
		for(int j=-5;j<lengthStep+5;j++){
			(*u[i])[j] = (*v[i])[j];
		}
	}
	
	//Since we assume spherical symmetry, the potential
	//Can only be found if you set to spherical symmetry
	//If would be nice if you add the other coordinate version
	if(spDimI == 2){
		if(p == 1){
			//findQuantumPotential();
			findPotential();
		}
	}
	
	if(helmeosFlag){
		//do sth
	}
	
	findPressure();
	
	return;
}

