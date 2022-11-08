#include "checkRho.h"
#include "../Variable/variable.h"
#include "../parameter.h"
#include "./boundary.h"

void checkRho(){
	
	//Assign the appropriate values of hydro variables
	for(int j=0;j<lengthStep;j++){
		if((*rho2)[j]<=rho2A){
			(*rho2)[j] = rho2A;
			(*vel2)[j] = vel2A;
			(*epsilon2)[j] = epsilon2A;
			(*temp2)[j] = temp2A;
		}
	}
	
	//Artifical 
	/*double vLimit = 0.1;
	for(int j=0;j<lengthStep;j++){
		if((*vel2)[j] > vLimit || (*vel2)[j] < -vLimit){
			//if((j+0.5)*dx > 2500){
				if((*vel2)[j] > 0.0){
					(*vel2)[j] = vLimit;
				}else{
					(*vel2)[j] = -vLimit;
				}
			//}
		}
	}*/
	
	boundary1D(rho2, false);
	boundary1D(vel2, true);
	boundary1D(epsilon2, false);
	boundary1D(temp2, false);
	
	if(xisotranFlag || flameFlag){
		//do sth
	}
	
	if(dmFlag){
		
		for(int j=0;j<lengthStep;j++){
			if((*rho1)[j]<=rho1A){
				(*rho1)[j] = rho1A;
				(*epsilon1)[j] = epsilon1A;
				if(runDMFlag){
					(*vel1)[j] = vel1A;
				}
			}
		}
		
		boundary1D(rho1, false);
		boundary1D(epsilon1, false);
		if(runDMFlag){
			boundary1D(vel1, true);
		}
		
	}
	
	return;
}

