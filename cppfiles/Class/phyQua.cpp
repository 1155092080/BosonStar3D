#include <iostream>
#include <cmath>
#include "phyQua.h"
#include "../parameter.h"

PhyQua::PhyQua(int size){
	mid = new double[size];
	left = new double[5];
	right = new double[5];
}

PhyQua::~PhyQua(){
	delete [] mid;
	delete [] left;
	delete [] right;
}

double& PhyQua::operator[](int index){
	if(index < -5 || index > lengthStep+4){
		std::cout << "Warning: Array index out of bound" << std::endl;
		exit(-1);
	}else if(index < 0){
		return left[abs(index)-1];
	}else if(index > lengthStep-1){
		return right[index-lengthStep];
	}else{
		return mid[index];
	}
}
