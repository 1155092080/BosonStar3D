#include <iostream> 
#include "sedov.h"
#include "../parameter.h"
#include "../Variable/variable.h"

using namespace std;

void sedov(){
	
	int endPt = lengthStep/(totalLength/500.0);
	
	for(int i=0;i<endPt;i++){
		(*epsilon2)[i] += 1.0 * ((double)(endPt-i)/(double)endPt);
	}
	
	return;
}
