#include <cmath> 
#include "runTestModel.h"
#include "../parameter.h"
#include "../Variable/variable.h"
#include "../Variable/definition.h"
#include "./boundary.h"

void runTestModel(){

	//Test cases

	if(testModel == 1){
		//Test 1
		gamma2 = 1.4;
		for(int i=0;i<lengthStep*3/10;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 0.75;
			(*p2)[i] = 1.0;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*3/10;i<lengthStep;i++){
			(*rho2)[i] = 0.125;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 0.1;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 2){
		//Test2
		gamma2 = 1.4;
		
		for(int i=0;i<lengthStep*5/10;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = -2.0;
			(*p2)[i] = 0.4;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*5/10;i<lengthStep;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 2.0;
			(*p2)[i] = 0.4;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 3){
		//Test3
		gamma2 = 1.4;
		
		for(int i=0;i<lengthStep*5/10;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 1000.0;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*5/10;i<lengthStep;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 0.01;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 4){
		//Test4
		gamma2 = 1.4;
		
		for(int i=0;i<lengthStep*4/10;i++){
			(*rho2)[i] = 5.99924;
			(*vel2)[i] = 19.5975;
			(*p2)[i] = 460.894;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*4/10;i<lengthStep;i++){
			(*rho2)[i] = 5.99242;
			(*vel2)[i] = -6.19633;
			(*p2)[i] = 46.0950;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 5){
		//Test5
		gamma2 = 1.4;
		
		for(int i=0;i<lengthStep*8/10;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = -19.5975;
			(*p2)[i] = 1000.0;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*8/10;i<lengthStep;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = -19.59745;
			(*p2)[i] = 0.01;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 6){
		//Test 6 (Diffusion test)
		for(int i=0;i<lengthStep*5/10;i++){
			(*rho2)[i] = 1.0 - 0.1;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 0.4 - 0.0001;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		double mm = 0.6571;
		double tt = 20.0;
		
		for(int j=lengthStep/2;j<(int)(mm*lengthStep);j++){
			double test = (*rho2)[0] * 1.0/2.0 * (sin(tt*((double)j/lengthStep + PI/(2*tt) - 0.5)) + 1.0);
			double test2 = (*p2)[0] * 1.0/2.0 * (sin(tt*((double)j/lengthStep + PI/(2*tt) - 0.5)) + 1.0);
			(*rho2)[j] = test;
			(*vel2)[j] = 0.0;
			(*p2)[j] = test2;
			(*epsilon2)[j] = (*p2)[j] / (*rho2)[j] / (gamma2 - 1.0);
		}
		
		for(int j=(int)(mm*lengthStep);j<lengthStep;j++){
			(*rho2)[j] = 0.0;
			(*vel2)[j] = 0.0;
			(*p2)[j] = 0.0;
			(*epsilon2)[j] = (*p2)[j] / (*rho2)[j] / (gamma2 - 1.0);
		}
		
		for(int i=0;i<lengthStep;i++){
			(*rho2)[i] += 0.001;
			(*p2)[i] += 0.0001;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	if(testModel == 7){
		//Sedov
		gamma2 = 1.4;
		
		for(int i=0;i<lengthStep*0.01;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 95492.96;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		for(int i=lengthStep*0.01;i<lengthStep;i++){
			(*rho2)[i] = 1.0;
			(*vel2)[i] = 0.0;
			(*p2)[i] = 1e-3;
			(*epsilon2)[i] = (*p2)[i] / (*rho2)[i] / (gamma2 - 1.0);
		}
		
		boundary1D(rho2, false);
		boundary1D(vel2, true);
		boundary1D(p2, false);
		boundary1D(epsilon2, false);
	}
	
	return;
}
