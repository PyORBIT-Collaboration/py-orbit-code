//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implement 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef QUADEMFIELD_HH_
#define QUADEMFIELD_HH_

#include "Python.h"

#include "BaseFieldSource.hh"

namespace OrbitUtils{
	
	class  QuadEMfield: public BaseFieldSource
	{
		public:
		
			QuadEMfield();
			~QuadEMfield();
		
			void getElectricMagneticField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);

			
		 	void Ex0(double val);
		 	void Ey0(double val);
		 	void Ez0(double val);
			
		 	void Bx0(double val);
		 	void By0(double val);
		 	void Bz0(double val);
			
		 	void x0(double val);
		 	void y0(double val);
		 	void z0(double val);
			
		 	void cxEx(double val);
		 	void cyEx(double val);
		 	void czEx(double val);
			
		 	void cxEy(double val);
		 	void cyEy(double val);
		 	void czEy(double val);
			
		 	void cxEz(double val);
		 	void cyEz(double val);
		 	void czEz(double val);
			
		 	void cxBx(double val);
		 	void cyBx(double val);
		 	void czBx(double val);
			
		 	void cxBy(double val);
		 	void cyBy(double val);
		 	void czBy(double val);
			
		 	void cxBz(double val);
		 	void cyBz(double val);
		 	void czBz(double val);	
			
			
			
		private:
			


		
			double Ex0v;
			double Ey0v;
			double Ez0v;
			
			double Bx0v;
			double By0v;
			double Bz0v;
			
			double x0v;
			double y0v;
			double z0v;
			
			double cxExv;
			double cyExv;
			double czExv;
			
			double cxEyv;
			double cyEyv;
			double czEyv;
			
			double cxEzv;
			double cyEzv;
			double czEzv;
			
			double cxBxv;
			double cyBxv;
			double czBxv;
			
			double cxByv;
			double cyByv;
			double czByv;
			
			double cxBzv;
			double cyBzv;
			double czBzv;
			
			
	};
};




#endif /*QUADEMFIELD_HH_*/



