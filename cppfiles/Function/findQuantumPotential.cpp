#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <complex>
#include <fftw3.h>
#include "findQuantumPotential.h"
#include "./boundary.h"
#include "./checkRho.h"
#include "./findPotential.h"
#include "../Variable/variable.h"
#include "../Variable/definition.h"
#include "../parameter.h"

#include <iostream>
using namespace std;

int binarySearchTest(double *arr, int l, int r, double x)
{ 
    if (r >= l) {
        int mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (arr[mid] == x) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid] > x) 
            return binarySearchTest(arr, l, mid - 1, x); 
  
        // Else the element can only be present 
        // in right subarray 
        return binarySearchTest(arr, mid + 1, r, x); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return l;
}

/*void findPotential2(){
	
	if(wGravityI){
		PhyQua *temp = new PhyQua(lengthStep);

		double dr = dx/(6.77193e-6 * 100)/Ra;
		(*temp)[0] = (4.0/3.0)*PI*((0.0+0.5)*dr)*((0.0+0.5)*dr)*((0.0+0.5)*dr)*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+0])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+0]);
		(*temp)[1] = 9.0*(*temp)[0] + (4.0/3.0)*PI*dr*((1.0+0.5)*dr)*((1.0+0.5)*dr)*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+1])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+1]);
		
		for(int j=2;j<lengthStep;j++){
			(*temp)[j] = (*temp)[j-2] + 4.0 * PI * (dr/3.0) *
						 (abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j-2])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j-2]) * (((j-2.0)+0.5)*dr) * (((j-2.0)+0.5)*dr) +
					  	 4.0*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j-1])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j-1])*(((j-1)+0.5)*dr)*(((j-1)+0.5)*dr) +
					  	 abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j]) * ((j+0.5)*dr) * ((j+0.5)*dr));
		}
		
		for(int j=0;j<lengthStep;j++){
			(*phip)[j] = (*temp)[j] / (((j+0.5)*dr)*((j+0.5)*dr));
		}
		
		boundary1D(phip, true);
		
		(*phi)[lengthStep-1] = 0.0;
		
		for(int j=lengthStep-2;j>=0;j--){
			(*phi)[j] = (*phi)[j+1] - ((*phip)[j] + (*phip)[j+1]) * dr/2.0;
		}
		
		boundary1D(phi, false);
		
		//potentialRelax
		PhyQua *phiNew = new PhyQua(lengthStep);
		PhyQua *error = new PhyQua(lengthStep);
		for(int n=0;n<relaxMax;n++){
			for(int j=0;j<lengthStep-1;j++){
				(*phiNew)[j] = 0.5 * ((*phi)[j+1] + (*phi)[j-1]) + 0.5 * dr/((j+0.5)*dr) * ((*phi)[j+1] - (*phi)[j-1]) - 
							   2.0*PI*dr*dr*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+lengthStep+j]);
			}
			(*phiNew)[lengthStep-1] = (*phi)[lengthStep-1];
			boundary1D(phiNew, false);
			
			for(int j=0;j<lengthStep-1;j++){
				(*error)[j] = ((*phiNew)[j] - (*phi)[j])/(*phi)[j];
			}
			(*error)[lengthStep-1] = 0.0;
			boundary1D(error, false);
			
			for(int j=-5;j<lengthStep+5;j++){
				(*phi)[j] = (*phiNew)[j];
			}
			
			for(int j=0;j<lengthStep;j++){
				(*phip)[j] = (-(*phi)[j+2] + 8.0*(*phi)[j+1] - 8.0*(*phi)[j-1] + (*phi)[j-2])/(1.2e1*dr);
			}
			boundary1D(phip, true);
			
			bool needRelax = false;
			
			for(int j=0;j<lengthStep;j++){
				if(abs((*error)[j]) > tolerance){
					needRelax = true;
				}
			}
			
			if(needRelax){
				continue;
			}else{
				break;
			}
			
		}
		
		delete temp;
		delete phiNew;
		delete error;
		
	}else{
		for(int j=-5;j<lengthStep+5;j++){
			(*phi)[j] = 0.0;
			(*phip)[j] = 0.0;
		}
	}
	
	return;
}*/

void findQuantumPotential(){
	
	if(qpFlag){

		//double Gdim = 1.65151;
		//double dr = dx/(6.77193e-6 * 100)/Ra;
		//findPotential2();
		findPotential();

		double x[lengthStep];
		for(int i=0;i<lengthStep;i++){
			x[i] = (i+0.5)*dx;
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					double xPos = sqrt((i+0.5-lengthStep)*dx*(i+0.5-lengthStep)*dx+(j+0.5-lengthStep)*dx*(j+0.5-lengthStep)*dx+(k+0.5-lengthStep)*dx*(k+0.5-lengthStep)*dx);
					int l = binarySearchTest(x, 0, lengthStep-1, xPos);

					if(l>lengthStep-2){
						l = lengthStep-2;
					}

					if(l<2){
						l = 2;
					}

					double x1 = x[l-2];
					double x2 = x[l-1];
					double x3 = x[l];
					double x4 = x[l+1];
					double y1 = (*phi)[l-2];
					double y2 = (*phi)[l-1];
					double y3 = (*phi)[l];
					double y4 = (*phi)[l+1];
					
					phi3D[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
							 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
							 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
							 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;

					phi3D[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k] = phi3D[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k]*(1e-7/5.59429e-55)*(5.02788e-34*1000.0)*pow(Ra, -2.0)*pow(omega, -2.0);
				}
			}
		}

		complex<double> I(0.0,1.0);

		double b = (lengthStep-0.5)*dx;
		b /= 6.77193e-6 * 100;
		b /= Ra;
		double a = -b;

		double dtDimensionless = dt;
		dtDimensionless /= 2.03017e5;
		dtDimensionless *= omega;

		//Action 1
		for(int j=0;j<2*lengthStep*2*lengthStep*2*lengthStep;j++){
			expo[j] = phi3D[j]*dtDimensionless/2.0/epsilon + delta*dtDimensionless/2.0/epsilon*abs(psiSave[j])*abs(psiSave[j]);
			psiSave[j] = (cos(expo[j])-sin(expo[j])*I)*psiSave[j];

			psiStarFourierFFTW[j][0] = psiSave[j].real();
			psiStarFourierFFTW[j][1] = psiSave[j].imag();
		}

		//Action 2
		fftw_plan plan = fftw_plan_dft_3d(2*lengthStep, 2*lengthStep, 2*lengthStep, psiStarFourierFFTW, psiStarFourierFFTW, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		fftw_cleanup();

		//Swap
		for(int i=0;i<lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourier[(i+lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[(i+lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=lengthStep;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourier[(i-lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[(i-lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+k][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+k][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=lengthStep;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+(j-lengthStep)*2*lengthStep+k][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+(j-lengthStep)*2*lengthStep+k][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<lengthStep;k++){
					psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k+lengthStep)].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k+lengthStep)].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=lengthStep;k<2*lengthStep;k++){
					psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k-lengthStep)].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k-lengthStep)].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=-lengthStep;i<lengthStep;i++){
			for(int j=-lengthStep;j<lengthStep;j++){
				for(int k=-lengthStep;k<lengthStep;k++){
					double muLr2 = (i*2.0*PI/(b-a))*(i*2.0*PI/(b-a))+(j*2.0*PI/(b-a))*(j*2.0*PI/(b-a))+(k*2.0*PI/(b-a))*(k*2.0*PI/(b-a));
					psiStarFourier[(i+lengthStep)*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+(k+lengthStep)] = (cos(epsilon*dtDimensionless*muLr2/2.0) - sin(epsilon*dtDimensionless*muLr2/2.0)*I)*psiStarFourier[(i+lengthStep)*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+(k+lengthStep)];
				}
			}
		}

		//Swap
		for(int i=0;i<lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourierFFTW[(i+lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[(i+lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}

		for(int i=lengthStep;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourierFFTW[(i-lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[(i-lengthStep)*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourier[i*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+k].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[i*2*lengthStep*2*lengthStep+(j+lengthStep)*2*lengthStep+k].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=lengthStep;j<2*lengthStep;j++){
				for(int k=0;k<2*lengthStep;k++){
					psiStarFourier[i*2*lengthStep*2*lengthStep+(j-lengthStep)*2*lengthStep+k].real(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][0]);
					psiStarFourier[i*2*lengthStep*2*lengthStep+(j-lengthStep)*2*lengthStep+k].imag(psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k][1]);
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=0;k<lengthStep;k++){
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k+lengthStep)][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k+lengthStep)][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}

		for(int i=0;i<2*lengthStep;i++){
			for(int j=0;j<2*lengthStep;j++){
				for(int k=lengthStep;k<2*lengthStep;k++){
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k-lengthStep)][0] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].real();
					psiStarFourierFFTW[i*2*lengthStep*2*lengthStep+j*2*lengthStep+(k-lengthStep)][1] = psiStarFourier[i*2*lengthStep*2*lengthStep+j*2*lengthStep+k].imag();
				}
			}
		}


		fftw_plan plan2 = fftw_plan_dft_3d(2*lengthStep, 2*lengthStep, 2*lengthStep, psiStarFourierFFTW, psiStarFourierFFTW, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan2);
		fftw_destroy_plan(plan2);
		fftw_cleanup();

		//Action 3
		for(int j=0;j<2*lengthStep*2*lengthStep*2*lengthStep;j++){

			psiSave[j].real(psiStarFourierFFTW[j][0]/(2*lengthStep*2*lengthStep*2*lengthStep));
			psiSave[j].imag(psiStarFourierFFTW[j][1]/(2*lengthStep*2*lengthStep*2*lengthStep));

			expo[j] = phi3D[j]*dtDimensionless/2.0/epsilon + delta*dtDimensionless/2.0/epsilon*abs(psiSave[j])*abs(psiSave[j]);
			psiSave[j] = (cos(expo[j])-sin(expo[j])*I)*psiSave[j];
		}

		for(int j=lengthStep;j<2*lengthStep;j++){
			double tmp = abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+j])*abs(psiSave[lengthStep*2*lengthStep*2*lengthStep+lengthStep*2*lengthStep+j]);
			tmp /= (Ra*Ra*Ra);
			tmp *= 1.619e-18;
			tmp /= 1000.0;
			tmp *= (M*1.9891e30);
			(*rho2)[j-lengthStep] = tmp;
		}

		for(int i=-5;i<lengthStep+5;i++){
			(*quantumPotential)[i] = 0.0;
		}

	}else{
		for(int i=-5;i<lengthStep+5;i++){
			(*quantumPotential)[i] = 0.0;
		}
	}
	
	if(qpDMFlag){
		for(int i=-5;i<lengthStep+5;i++){
			(*quantumPotential1)[i] = 0.0;
		}
	}else{
		for(int i=-5;i<lengthStep+5;i++){
			(*quantumPotential1)[i] = 0.0;
		}
	}
	
	return;
}
