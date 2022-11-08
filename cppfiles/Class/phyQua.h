#ifndef __PHYQUA_H__
#define __PHYQUA_H__

class PhyQua{
	
	public:
		PhyQua(int);
		~PhyQua();
		
		double& operator[] (int);
	
	private:
		double *mid;
		double *left;
		double *right;
	
};

#endif

