#include <cmath>
#include <fstream>
#include "eosTable.h"
#include "../Variable/variable.h"
#include "../Variable/definition.h"
#include "../parameter.h"

using namespace std;

void eosTable(){
	
	int start = -20;
	int end = -5;
	
	int k = 0;
	
	for(int i=start;i<=end;i++){
		for(int j=0;j<eosTableWidth;j++){
			double den = (pow(10.0, i+1) - pow(10.0, i))*j/(double)eosTableWidth + pow(10.0, i);
			
			//fermimo
			double dlfmmo = pow(6.0*PI*PI*hBar*hBar*hBar*den*ye2/(gs2*mb2*me2*me2*me2), 1.0/3.0);
			
			double pree2;
			if(dlfmmo <= 6.0e-2){
				pree2 = (gs2*me2*me2*me2*me2)/(30.0*PI*PI*hBar*hBar*hBar) * (pow(dlfmmo, 5.0) - 0.5*pow(dlfmmo, 7.0)/14.0 + 5.0*pow(dlfmmo, 9.0)/24.0);				
			}else{
				pree2 = (gs2*me2*me2*me2*me2)/(16.0*PI*PI*hBar*hBar*hBar) * (dlfmmo*pow(1.0+dlfmmo*dlfmmo, 1.0/2.0) * (2.0*dlfmmo*dlfmmo/3.0 - 1.0) + log(dlfmmo + sqrt(1.0+dlfmmo*dlfmmo)));
			}
			
			eosTableList[k*eosTableWidth+j][0] = pree2;
			eosTableList[k*eosTableWidth+j][1] = den;
			
		}
		
		k++;
	}
	
	return;
}

void eosTableBoson(){
	
	int start = -15;
	int end = 0;
	
	int k = 0;
	
	double K  = 8.0*PI*scattering*hBar*hBar/4.0/mBoson/mBoson/mBoson;
	
	for(int i=start;i<=end;i++){
		for(int j=0;j<eosTableWidth;j++){
			double den = (pow(10.0, i+1) - pow(10.0, i))*j/(double)eosTableWidth + pow(10.0, i);
			double pree2 = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*den)-1.0) * (sqrt(1.0+12.0*K*den)-1.0);
			
			eosTableListBoson[k*eosTableWidth+j][0] = pree2;
			eosTableListBoson[k*eosTableWidth+j][1] = den;
			
		}
		
		k++;
	}
	
	return;
}

void eosTableAPR(){
	
	ifstream fin("./Library/aprEOS.table");
	
	double tmp;
	
	for(int i=0;i<38400;i++){
		fin >> tmp >> eosTableList[i][1] >> tmp >> eosTableList[i][0] >> tmp;
		eosTableList[i][0] *= 1.80139e-39/6.2415e-34;
		eosTableList[i][1] *= tmp*1.79e-30*1000*1e45/1e6 * 1.619e-18;
	}
	
	fin.close();
	
	return;
}

void eosTableNewAPR(){
	
	ifstream fin("./Library/newaprEOS.table");
	
	double tmp;
	
	for(int i=0;i<38400;i++){
		fin >> tmp >> eosTableList[i][1] >> tmp >> eosTableList[i][0] >> tmp;
		eosTableList[i][0] *= 1.80139e-39/6.2415e-34;
		eosTableList[i][1] *= tmp*1.79e-30*1000*1e45/1e6 * 1.619e-18;
	}
	
	fin.close();
	
	return;
}
