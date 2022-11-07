
#include <iostream>


using namespace std;


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the pressure given the density    !
! of a completely degenerate ideal fermi gas. It is used only !
! in solving for the initial star model. Do not confuse it    !
! with the subroutine findpressure, which aims to update the  !
! pressure after hydrodynamic evolution                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void getRhoEOSRToP(double &iniP, double iniRho, int type){
	if(type == 1){
		if(bosonDMEOS){
			double K  = 8.0*PI*scattering*hBar*hBar/4.0/mBoson/mBoson/mBoson;
			iniP = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*iniRho)-1.0) * (sqrt(1.0+12.0*K*iniRho)-1.0);
		}else{
			iniP = k1*pow(iniRho, gamma1);
		}
	}else if(type == 2){
		if(fermiEOS){
			double dlfmmo = pow(6.0*PI*PI*hBar*hBar*hBar*iniRho*ye2/(gs2*mb2*me2*me2*me2), 1.0/3.0);
			
			if(dlfmmo <= 6.0e-2){
				iniP = (gs2*me2*me2*me2*me2)/(30.0*PI*PI*hBar*hBar*hBar) * (pow(dlfmmo, 5.0) - 5.0*pow(dlfmmo, 7.0)/14.0 + 5.0*pow(dlfmmo, 9.0)/24.0);				
			}else{
				iniP = (gs2*me2*me2*me2*me2)/(16.0*PI*PI*hBar*hBar*hBar) * (dlfmmo*sqrt(1.0+dlfmmo*dlfmmo) * (2.0*dlfmmo*dlfmmo/3.0 - 1.0) + log(dlfmmo + sqrt(1.0+dlfmmo*dlfmmo)));
			}
			
		}else if(bosonEOS){
			double K  = 8.0*PI*scattering*hBar*hBar/4.0/mBoson/mBoson/mBoson;
			iniP = 1.0/(36.0*K) * (sqrt(1.0+12.0*K*iniRho)-1.0) * (sqrt(1.0+12.0*K*iniRho)-1.0);
		}else if(aprEOS || newAPREOS){
			if(iniRho < eosTableList[1][1] || iniRho > eosTableList[16*eosTableWidth-2][1] || std::isnan(iniRho)){
				iniP = 0.0;
				return;
			}
			
			double x1,x2,x3,x4,y1,y2,y3,y4;
			int i = binarySearch2(eosTableList, 0, 16*eosTableWidth-1, iniRho);
			
			x1 = eosTableList[i-2][1];
			x2 = eosTableList[i-1][1];
			x3 = eosTableList[i][1];
			x4 = eosTableList[i+1][1];
			y1 = eosTableList[i-2][0];
			y2 = eosTableList[i-1][0];
			y3 = eosTableList[i][0];
			y4 = eosTableList[i+1][0];
			
			iniP = (iniRho-x2)*(iniRho-x3)*(iniRho-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
				   (iniRho-x1)*(iniRho-x3)*(iniRho-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
				   (iniRho-x1)*(iniRho-x2)*(iniRho-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
				   (iniRho-x1)*(iniRho-x2)*(iniRho-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
			
		}else{
			iniP = k2*pow(iniRho, gamma2);
		}
	}
	
	return;
}
//RK4 eq.

//g'=h (f'')
double j(double h){
	return h;
}

//h'=i (f''')
double k(double i){
	return i;
}

//f'=g (f')
double l_rk(double g){
	return g;
}

//i'=... (f''''), following page 17. Generally, mu=n0^(1/2)x
double m(double x, double f, double f2, double g, double h, double i){
	return 2.0*f*f + 2.0*f*f2 - 4.0/x *i + 10.0*g*h/f/x - 6.0*g*g*g/f/f/x + 3.0*i*g/f + 2.0*h*h/f - 7.0*g*g*h/f/f + 3.0*g*g*g*g/f/f/f + 2.0*f*mu*(h+2.0*g/x);
}

// dM_NM/dr=4*pi*r^2(rho_NM)
double der1(double r, double y3){
	return 4.0*PI*r*r*y3;
}
// dp_NM/dr=-G*M_total*rho_NM/r^2
double der2(double r, double y1, double y3, double y4){
	return -((y1+y4)*y3)/(r*r);
}

double der3(){
	return 0.0;
}

/*
This loop function evolves all f, f', f'' and f''' by RK4 method, and saves in soll, while f2=0
*/
void loop(double *soll, double x, double f, double g, double h, double i, double dx){
	double k1 = l_rk(g);
	double l1 = k(i);
	double m1 = j(h);
	double n1 = m(x,f,0.0,g,h,i);
	
	double k2 = l_rk(g+dx/2.0*m1);
	double l2 = k(i+dx/2.0*n1);
	double m2 = j(h+dx/2.0*l1);
	double n2 = m(x+dx/2.0, f+dx/2.0*k1, 0.0, g+dx/2.0*m1, h+dx/2.0*l1, i+dx/2.0*n1);
	
	double k3 = l_rk(g+dx/2.0*m2);
	double l3 = k(i+dx/2.0*n2);
	double m3 = j(h+dx/2.0*l2);
	double n3 = m(x+dx/2.0, f+dx/2.0*k2, 0.0, g+dx/2.0*m2, h+dx/2.0*l2, i+dx/2.0*n2);
	
	double k4 = l_rk(g+dx*m3);
	double l4 = k(i+dx*n3);
	double m4 = j(h+dx*l3);
	double n4 = m(x+dx, f+dx*k3, 0.0, g+dx*m3, h+dx*l3, i+dx*n3);
	
    soll[0] = f+dx/6.0*(k1+2.0*k2+2.0*k3+k4); // f = f + dx/6*(some g)
    soll[1] = h+dx/6.0*(l1+2.0*l2+2.0*l3+l4); // h = h + dx/6*(some i)
    soll[2] = g+dx/6.0*(m1+2.0*m2+2.0*m3+m4); // g = g + dx/6*(some h)
    soll[3] = i+dx/6.0*(n1+2.0*n2+2.0*n3+n4); // i = i + dx/6*(some f'''')
    
    return;
}
/*
This loop function evolves all f, f', f'' and f''' by RK4 method, and saves in soll, while f2 is const
*/
void loop2(double *soll, double x, double f, double f2, double g, double h, double i, double dx){
	double k1 = l_rk(g);
	double l1 = k(i);
	double m1 = j(h);
	double n1 = m(x,f,f2,g,h,i);
	
	double k2 = l_rk(g+dx/2.0*m1);
	double l2 = k(i+dx/2.0*n1);
	double m2 = j(h+dx/2.0*l1);
	double n2 = m(x+dx/2.0, f+dx/2.0*k1, f2, g+dx/2.0*m1, h+dx/2.0*l1, i+dx/2.0*n1);
	
	double k3 = l_rk(g+dx/2.0*m2);
	double l3 = k(i+dx/2.0*n2);
	double m3 = j(h+dx/2.0*l2);
	double n3 = m(x+dx/2.0, f+dx/2.0*k2, f2, g+dx/2.0*m2, h+dx/2.0*l2, i+dx/2.0*n2);
	
	double k4 = l_rk(g+dx*m3);
	double l4 = k(i+dx*n3);
	double m4 = j(h+dx*l3);
	double n4 = m(x+dx, f+dx*k3, f2, g+dx*m3, h+dx*l3, i+dx*n3);
	
    soll[0] = f+dx/6.0*(k1+2.0*k2+2.0*k3+k4);
    soll[1] = h+dx/6.0*(l1+2.0*l2+2.0*l3+l4);
    soll[2] = g+dx/6.0*(m1+2.0*m2+2.0*m3+m4);
    soll[3] = i+dx/6.0*(n1+2.0*n2+2.0*n3+n4);
    
    return;
}
/*
This loop evolves M_NM, P_NM, Pho_NM and M_DM using newtonian formulation and considering EOS. But M_DM is considered const when calculating NM and updated by RK1 at last
*/
void loop3(double *sol, double r, double y1, double y2, double y3, double y4, double y5, double dr){
	
	// Step1:
	double iniP2, iniRho2;
	
	double k1 = der1(r, y3); // dM_NM/dr
	double l1 = der2(r, y1, y3, y4); // dP_NM/dr
	double m1 = der3(); // dRho_NM/dr, but it is 0 since we don't update it by derivative while we use EOS
	
	// Step2:
	iniP2 = y2+dr/2.0*l1; // update half step of P_NM
	getRhoEOSPToR(iniRho2, iniP2, 2); // Obtain corresponding Rho_NM from P_NM
	y3 = iniRho2; // update half step of Rho_NM
	
	double k2 = der1(r+dr/2.0, y3+dr/2.0*m1); // dM_NM/dr at r+dr/2 with updated values. Here we set y3 by adding dr/2*m1, but y3 is already updated in last step
	double l2 = der2(r+dr/2.0, y1+dr/2.0*k1, y3+dr/2.0*m1, y4); // dP_NM/dr at r+dr/2 with updated detivatives
	double m2 = der3(); // dRho_NM/dr = 0 by the same reason.
	
	// Step3:
	iniP2 = y2+dr/2.0*l2; // update half step of P_NM
	getRhoEOSPToR(iniRho2, iniP2, 2); // Obtain corresponding Rho_NM from P_NM
	y3 = iniRho2; // update half step of Rho_NM
	
	// Same as above "step 2" but we use new derivatives calculated above
	double k3 = der1(r+dr/2.0, y3+dr/2.0*m2); // dM_NM/dr at r+dr/2 with updated values. Here we set y3 by adding dr/2*m1, but y3 is already updated in last step
	double l3 = der2(r+dr/2.0, y1+dr/2.0*k2, y3+dr/2.0*m2, y4); // dP_NM/dr at r+dr/2 with updated detivatives
	double m3 = der3(); // dRho_NM/dr = 0 by the same reason.
	
	// Step4:
	iniP2 = y2+dr*l3; // update full step of P_NM
	getRhoEOSPToR(iniRho2, iniP2, 2); // Obtain corresponding Rho_NM from P_NM
	y3 = iniRho2; // update full step of Rho_NM
	
	double k4 = der1(r+dr, y3+dr*m3); // dM_NM/dr at r+dr with updated values.
	double l4 = der2(r+dr, y1+dr*k3, y3+dr*m3, y4); // dP_NM/dr at r+dr with updated detivatives
	double m4 = der3(); // dRho_NM/dr = 0 by the same reason.
	
	// Summarize: only M_NM and P_NM is updated with derivatives. Rho_NM is from EOS and M_DM is by RK1
	iniP2 = y2+dr/6.0*(l1+2.0*l2+2.0*l3+l4);
	getRhoEOSPToR(iniRho2, iniP2, 2);
	y3 = iniRho2;
	
    sol[0] = y1+dr/6.0*(k1+2.0*k2+2.0*k3+k4); // Updated M_NM
    sol[1] = y2+dr/6.0*(l1+2.0*l2+2.0*l3+l4); // Updated P_NM
    sol[2] = y3+dr/6.0*(m1+2.0*m2+2.0*m3+m4); // Updated Rho_NM
    sol[3] = y4+4.0*PI*r*r*y5*dr;             // Updated M_DM
	
	return;
}


void getRhoBoson2F(){
	//Initialize
	double *x = new double [stepNum]; //x coordinate
	double *f = new double [stepNum]; //f1
	double *g = new double [stepNum]; //df1/dx
	double *h = new double [stepNum]; //d2f1/dx2
	double *i = new double [stepNum]; //d3f1/dx3
	double *r = new double [stepNum]; //r coordinate
	
	double *f2 = new double [stepNum]; //f2
	double *temp1 = new double [stepNum]; //M2
	double *temp2 = new double [stepNum]; //P2
	double *temp3 = new double [stepNum]; //rho2
	double *temp4 = new double [stepNum]; //M1
	double *temp5 = new double [stepNum]; //rho1
	
	double *tempRho = new double [stepNum]; //useless now
	double *shell = new double[stepNum]; //shell
	double sol[4];
	double sol2[4];
	
	/*Question: Why 30.0-1e-15?*/
	double dxrk = (30.0-1e-15)/(stepNum-1.0); //dx from 30-1e-15 to 30 with stepNum
	for(int i=0;i<stepNum;i++){
		x[i] = dxrk*i + 1e-15;
	}
	//Part 1: solve for f2=0 (double derivative of f is 0)
	/*Question: Why these values?*/
	double start = -15.0;
	double end = -0.0001;

	while(abs(start-end) > 1e-15){
		
		/*Question: What does these A stand for?*/ 
		double A[7];
		double dA = (end-start)/6.0;
		for(int i=0;i<7;i++){
			A[i] = dA*i + start;
		}
		double maxFarPoint = -1.0;
		int maxFarPointAA = -1;

		//Initial condition for the shooting method. With f, f', f'''=0 and f''=A2
		for(int aa=0;aa<7;aa++){
			f[0] = 1.0;
			h[0] = A[aa];
			g[0] = 0.0;
			i[0] = 0.0;

			// For each trial of A2, evolve f, f', f'', f''' through spatial space by RK4 method
			for(int qq=0;qq<stepNum-1;qq++){
				loop(sol,x[qq],f[qq],g[qq],h[qq],i[qq], dxrk);
				f[qq+1] = sol[0];
				h[qq+1] = sol[1];
				g[qq+1] = sol[2];
				i[qq+1] = sol[3];
			}

			/* Question: What is it ??? Sol: Find the minimum value and corresponding index*/
			int farPoint = -1;
			double minValue = 1.01;
			for(int qq=0;qq<stepNum;qq++){
				if(f[qq] < minValue && f[qq] > 0){
					farPoint = qq;
					minValue = f[qq];
				}else{
					break;
				}
			}
			if(farPoint > maxFarPoint){
				maxFarPoint = farPoint;
				maxFarPointAA = aa;
			}
		}
		// Here, we aim to find the best range of A2 value
		start = A[max(maxFarPointAA-1, 0)];
		end = A[min(maxFarPointAA+1, 6)];
	}
	// After above loop, the range of A would be less than 1e-15 as set in condition of while loop
	//Best profile
	f[0] = 1.0;
	h[0] = start;
	g[0] = 0.0;
	i[0] = 0.0;
	
	for(int qq=0;qq<stepNum-1;qq++){
		loop(sol,x[qq],f[qq],g[qq],h[qq],i[qq], dxrk);
		f[qq+1] = sol[0];
		h[qq+1] = sol[1];
		g[qq+1] = sol[2];
		i[qq+1] = sol[3];
	}
	
	int farPoint = -1;
	double minValue = 1.01;
	for(int qq=0;qq<stepNum;qq++){
		if(f[qq] < minValue && f[qq] > 0){
			farPoint = qq;
			minValue = f[qq];
		}else{
			break;
		}
	}

	for(int qq=0;qq<stepNum;qq++){
		shell[qq] = f[qq]*x[qq]*x[qq];
	}
	/*Question: Why we can calculate n in this way?*/
	double volume = shell[0] + shell[farPoint];
	for(int i=1;i<farPoint;i++){
		if(i%3 == 0){
			volume += 2.0*shell[i];
		}else{
			volume += 3.0*shell[i];
		}
	}
	volume = 3.0*dxrk*volume/8.0;

	// Canvert all value back to original unit, as P17
	double n0 = pow(1.0/volume, 4.0);
	double chi = mu/sqrt(n0);
	double mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
	double b = hBar*hBar/2.0/mStar/mBoson/mBoson;
	double rho0 = mStar*n0/4.0/PI/b/b/b;
	for(int i=0;i<stepNum;i++){
		r[i] = x[i] * volume * b;
		if(i > farPoint){
			tempRho[i] = 1e-6 * rho0;
		}else{
			tempRho[i] = f[i] * rho0;
		}
	}
	double drrk = r[1] - r[0];
	
	cout << "Part1: " << x[farPoint] << "\t" << mStar << "\t" << n0 << endl;


	//Part 2: SOlve for f2!=0, now we need to do iterations until error less than a value
	
	//init error
	double error = 999.9;

	while(abs(error) > 1e-4){
		
		start = -15.0;
		end = -0.0001;
		// initial f2 for NM center density
		double iniF2 = rho2C*4.0*PI*b*b*b/mStar/n0;
		
		/*Question: Dencity of NM at center and atmosphere, but why?*/
		double p2C, p2A;
		getRhoEOSRToP(p2C, rho2C, 2);
		getRhoEOSRToP(p2A, rho2A, 2);
		
		while(abs(start-end) > 1e-15){
			
			double A[7];
			double dA = (end-start)/6.0;
			for(int i=0;i<7;i++){
				A[i] = dA*i + start;
			}
			
			double maxFarPoint = -1.0;
			int maxFarPointAA = -1;
			
			for(int aa=0;aa<7;aa++){
				
				// Recall 1,2,3,4,5 are m2, p2, rho2, m1, rho1
				temp1[0] = 0.0;
				temp2[0] = p2C;
				temp3[0] = rho2C;
				temp4[0] = 0.0;
				temp5[0] = rho0;
				
				f[0] = 1.0;
				f2[0] = iniF2;
				h[0] = A[aa];
				g[0] = 0.0;
				i[0] = 0.0;
				
				for(int qq=0;qq<stepNum-1;qq++){
					
					loop2(sol,x[qq],f[qq],f2[qq],g[qq],h[qq],i[qq],dxrk);
					loop3(sol2,r[qq],temp1[qq],temp2[qq],temp3[qq],temp4[qq],temp5[qq],drrk);
					f[qq+1] = sol[0];
					h[qq+1] = sol[1];
					g[qq+1] = sol[2];
					i[qq+1] = sol[3];
					
					temp1[qq+1] = sol2[0];
					temp2[qq+1] = sol2[1];
					temp3[qq+1] = sol2[2];
					temp4[qq+1] = sol2[3];
					
					//Because only when negative slop means stable states, and density should be larger than 0 
					if(f[qq+1] >= f[qq] || f[qq+1] < 0.0){
						f[qq+1] = f[qq];
					}
					//Again, only negative slope is acceptable and the center density should be larger than atmosphere.
					if(std::isnan(temp3[qq+1]) || temp3[qq+1] < rho2A || temp3[qq+1] >= temp3[qq]){
						temp3[qq+1] = rho2A;
						temp2[qq+1] = p2A;
					}
					
					f2[qq+1] = temp3[qq+1]*4.0*PI*b*b*b/mStar/n0;
					temp5[qq+1] = f[qq+1]*rho0;
				}
				
				int farPoint = -1;
				double minValue = 1.01;
				for(int qq=0;qq<stepNum;qq++){
					if(f[qq] < minValue && f[qq] > 0){
						farPoint = qq;
						minValue = f[qq];
					}else{
						break;
					}
				}
				
				if(farPoint > maxFarPoint){
					maxFarPoint = farPoint;
					maxFarPointAA = aa;
				}
				
			}
			
			start = A[max(maxFarPointAA-1, 0)];
			end = A[min(maxFarPointAA+1, 6)];
			
		}
		
		//Best profile: After reducing A range, do one more iteration for final profile
		temp1[0] = 0.0;
		temp2[0] = p2C;
		temp3[0] = rho2C;
		temp4[0] = 0.0;
		temp5[0] = rho0;
		
		f[0] = 1.0;
		f2[0] = iniF2;
		h[0] = start;
		g[0] = 0.0;
		i[0] = 0.0;
		
		for(int qq=0;qq<stepNum-1;qq++){
			// By RK4, evolve DM part first with f2 a constant.
			loop2(sol,x[qq],f[qq],f2[qq],g[qq],h[qq],i[qq],dxrk);
			// By RK4, evolve NM part first with M_DM a constant. And then evolve M_DM with RK1
			loop3(sol2,r[qq],temp1[qq],temp2[qq],temp3[qq],temp4[qq],temp5[qq],drrk);
			f[qq+1] = sol[0];
			h[qq+1] = sol[1];
			g[qq+1] = sol[2];
			i[qq+1] = sol[3];
			
			temp1[qq+1] = sol2[0];
			temp2[qq+1] = sol2[1];
			temp3[qq+1] = sol2[2];
			temp4[qq+1] = sol2[3];
			
			if(f[qq+1] >= f[qq] || f[qq+1] < 0.0){
				f[qq+1] = f[qq];
			}
			
			if(std::isnan(temp3[qq+1]) || temp3[qq+1] < rho2A || temp3[qq+1] >= temp3[qq]){
				temp3[qq+1] = rho2A;
			}
			
			f2[qq+1] = temp3[qq+1]*4.0*PI*b*b*b/mStar/n0;
			temp5[qq+1] = f[qq+1]*rho0;
			if(f[qq+1] < 1e-6 || f[qq+1] > 1.0 || std::isnan(f[qq+1])){
				temp5[qq+1] = 1e-6*rho0;
			}
		}
		
		int farPoint = -1;
		double minValue = 1.01;
		for(int qq=0;qq<stepNum;qq++){
			if(f[qq] < minValue && f[qq] > 0){
				farPoint = qq;
				minValue = f[qq];
			}else{
				break;
			}
		}
		
		for(int qq=0;qq<stepNum;qq++){
			shell[qq] = f[qq]*x[qq]*x[qq];
		}
		
		double volume = shell[0] + shell[farPoint];
		for(int i=1;i<farPoint;i++){
			if(i%3 == 0){
				volume += 2.0*shell[i];
			}else{
				volume += 3.0*shell[i];
			}
		}
		volume = 3.0*dxrk*volume/8.0;
		
		double newN0 = pow(1.0/volume, 4.0);
		error = (newN0 - n0)/n0;
		n0 = newN0;
		chi = mu/sqrt(n0);
		mStar = sqrt(chi*hBar*hBar/4.0/scattering/mBoson);
		b = hBar*hBar/2.0/mStar/mBoson/mBoson;
		rho0 = mStar*n0/4.0/PI/b/b/b;
		/*Question: How to calculate n0, rho0 and the relationship to volume*/
		for(int i=0;i<stepNum;i++){
			r[i] = x[i] * volume * b;
			if(i > farPoint){
				tempRho[i] = 1e-6 * rho0;
			}else{
				tempRho[i] = f[i] * rho0;
			}
		}
		drrk = r[1] - r[0];
		
		rho1C = rho0;
		rho1A = rho1C * 1e-6;
		double scale = b/pow(n0, 1.0/4.0);
		
		cout << "Part2: " << x[farPoint] << "\t" << mStar << "\t" << scale << "\t" << error << endl;
		
	}
	///////////////
	
	double **tt = new double *[stepNum];
	for(int qq=0;qq<stepNum;qq++){
		tt[qq] = new double[3];
	}
	
	for(int qq=0;qq<stepNum;qq++){
		tt[qq][0] = r[qq]; //r
		tt[qq][1] = tempRho[qq]; //Boson
		tt[qq][2] = temp3[qq]; //NM
	}
	
	///////////////
	
	//Deallocate
	
	delete [] x;
	delete [] f;
	delete [] g;
	delete [] h;
	delete [] i;
	delete [] shell;
	delete [] r;
	delete [] tempRho;
	
	delete [] f2;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] temp4;
	delete [] temp5;
	
	//Here we do a polynomic interpolation. Check http://www.appstate.edu/~grayro/comphys/lecture4_11.pdf
	// We want to remap what we calculated to a uniform grid
	for(int i=0;i<lengthStep;i++){
		double xPos = (i+0.5)*dx;
		for(int j=0;j<stepNum;j++){
			if(xPos < tt[j][0]){
				double x1 = tt[j-2][0];
				double x2 = tt[j-1][0];
				double x3 = tt[j][0];
				double x4 = tt[j+1][0];
				double y1 = tt[j-2][1];
				double y2 = tt[j-1][1];
				double y3 = tt[j][1];
				double y4 = tt[j+1][1];
				double z1 = tt[j-2][2];
				double z2 = tt[j-1][2];
				double z3 = tt[j][2];
				double z4 = tt[j+1][2];
				
				(*rho1)[i] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * y1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * y2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * y3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * y4;
						 
				(*rho2)[i] = (xPos-x2)*(xPos-x3)*(xPos-x4)/((x1-x2)*(x1-x3)*(x1-x4)) * z1 +
						 (xPos-x1)*(xPos-x3)*(xPos-x4)/((x2-x1)*(x2-x3)*(x2-x4)) * z2 +
						 (xPos-x1)*(xPos-x2)*(xPos-x4)/((x3-x1)*(x3-x2)*(x3-x4)) * z3 +
						 (xPos-x1)*(xPos-x2)*(xPos-x3)/((x4-x1)*(x4-x2)*(x4-x3)) * z4;
				
				if((*rho1)[i] < rho1A || std::isnan((*rho1)[i])){
					(*rho1)[i] = rho1A;
				}
				
				if((*rho2)[i] < rho2A || std::isnan((*rho2)[i])){
					(*rho2)[i] = rho2A;
				}
				
				break;
			}
		}
	}
	
	//Hot fix, or a cutoff of the minimum DM and NM density according to the value at atmosphere
	rhoMinDM = 1.1 * rho1A;
	rhoMinNM = 1.1 * rho2A;
	
	//Hot fix 2, make sure the remapping doesn't create unstable profile
	for(int i=0;i<lengthStep-1;i++){
		if((*rho2)[i+1] >= (*rho2)[i] || (*rho2)[i+1] < 0.0){
			(*rho2)[i+1] = rho2A;
		}
		
		if((*rho1)[i+1] >= (*rho1)[i] || (*rho1)[i+1] < 0.0){
			(*rho1)[i+1] = rho1A;
		}
	}
	
	boundary1D(rho1, false);
	boundary1D(rho2, false);
	
	for(int i=0;i<stepNum;i++){
		delete [] tt[i];
	}
	delete [] tt;
	
	return;
}