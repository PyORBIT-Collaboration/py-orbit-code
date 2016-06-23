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
#ifndef CONE_LF_MODE_HH_
#define CONE_LF_MODE_HH_

#include "Python.h"

#include "BaseLaserFieldSource.hh"
#include "FieldOrientation.hh"



namespace OrbitUtils{
	
	class  ConeLFmode: public BaseLaserFieldSource{
		public:
		
			ConeLFmode(double Cnm,double wx,double wy,double f_x,double f_y,double la);
			~ConeLFmode();
		
			void getLaserElectricMagneticField(double x, double y, double z, double t, 
					tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
					tcomplex& H_x, tcomplex& H_y, tcomplex& H_z);

					
			tcomplex getNonOrientedU( double x, double y, double z, double t);
		
			
			double getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t);
			
			void setLaserFieldOrientation(double x_0, double y_0, double z_0,
										double k_x, double k_y, double k_z,
										double m_x, double m_y, double m_z,
										double n_Ex, double n_Ey, double n_Ez);
			
			

			
		private:
			
			

			double HermiteGaussianOmega(double m, double x, double y, double z, double px, double py, double pz, double t);

			
			

			
			  //Laser parameters
			  tcomplex U;
			  double Laser_lambda;
			  double k;
			  double rx;
			  double ry;
			  double ax;
			  double ay;
			
			  //parameters of orientation of laser field
			  double nHx, nHy, nHz, nEx, nEy, nEz;	  
			  FieldOrientation*   orient; 
			

	};
};



#endif /*CONE_LF_MODE_HH_*/
