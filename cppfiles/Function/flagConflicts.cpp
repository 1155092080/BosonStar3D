#include <iostream>
#include "flagConflicts.h"
#include "../parameter.h"
using namespace std;

//Check whether there are conflicts on flags
bool flagConflicts(){
	bool conflicts = false;
	
	if(!dmFlag && runDMFlag){
		cout << "Error, the dmFlag is " << dmFlag;
		cout << " and the runDMFlag is " << runDMFlag;
		cout << " which is inconsistent" << endl;
		conflicts = true;
	}
	
	if(!dmFlag && levelSetFlagDM){
		cout << "Error, the dmFlag is " << dmFlag;
		cout << " and the levelSetFlagDm is " << levelSetFlagDM;
		cout << " which is inconsistent" << endl;
		conflicts = true;
	}
	
	if(!runDMFlag && levelSetFlagDM){
		cout << "Error, the runDMFlag is " << runDMFlag;
		cout << " and the levelSetFlagDM is " << levelSetFlagDM;
		cout << " which is inconsistent" << endl;
		conflicts = true;
	}
	
	/*if(sedovFlag && !helmeosFlag){
		cout << "Error, the sedovFlag is " << sedovFlag;
		cout << " and the helmeosFlag is " << helmeosFlag;
		cout << " which is inconsistent" << endl;
		conflicts = true;
	}*/
	
	if(movingGridFlag && runDMFlag){
		cout << "Error, the movingGridFlag is " << movingGridFlag;
		cout << " and the wGravityI is " << wGravityI;
		cout << " which is inconsistent" << endl;
		conflicts = true;
	}
	
	if(testModel != 0 && wGravityI == 1){
		cout << "Error, the testModel is " << testModel;
		cout << " and the wGravityI is " << wGravityI;
		cout << " which is inconsistent, end simulation" << endl;
		conflicts = true;
	}
	
	return conflicts;
}
