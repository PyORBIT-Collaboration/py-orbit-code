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
#ifndef FROISSART_STORA_LF_HH_
#define FROISSART_STORA_LF_HH_

#include "Python.h"

#include "BaseLaserFieldSource.hh"



namespace OrbitUtils{
	
	class  FroissartStoraLF: public BaseLaserFieldSource{
		public:
		
			FroissartStoraLF(double Om,double G,double El);
			~FroissartStoraLF();
		
			tcomplex getLaserEMField(double x, double y, double z, double t, 
					tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
					tcomplex& H_x, tcomplex& H_y, tcomplex& H_z);

					
			tcomplex getNonOrientedU(double t);
		
			
			double getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t);
			
			void setLaserFieldPolarization(double n_Ex, double n_Ey, double n_Ez);
			


			
		private:
			
			
//			double FroissartStoraTestOmega(double x, double y, double z, double px, double py, double pz, double t);
//			double HermiteGaussianOmega(double m, double x, double y, double z, double px, double py, double pz, double t);
//			tcomplex E_FroissartStora(double x, double y, double z, double t);
			
			

			
			  //Laser parameters

			  double Omega;
			  double Gamma;
			  double amp_ELas;


			
			  //parameters of polarization of laser field
			  double nEx,nEy,nEz;

			  

			

	};
};



#endif /*FROISSART_STORA_LF_HH_*/
