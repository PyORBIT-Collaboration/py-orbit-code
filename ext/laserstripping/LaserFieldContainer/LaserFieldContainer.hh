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
#ifndef LASERFIELDCONTAINER_HH_
#define LASERFIELDCONTAINER_HH_



#include "Python.h"

#include "BaseLaserFieldSource.hh"
#include <vector>


namespace OrbitUtils{
	
	class  LaserFieldContainer: public BaseLaserFieldSource{
		public:
		
			LaserFieldContainer();
			~LaserFieldContainer();
		
			void getLaserElectricField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
			void getLaserMagneticField(double x, double y, double z, double t, tcomplex& f_x, tcomplex& f_y, tcomplex& f_z);
	
			double getFrequencyOmega(double m,double x, double y, double z, double px, double py, double pz, double t);
		
			
		private:
			
			std::vector<BaseLaserFieldSource*>	ref;
			

	};
};



#endif /*LASERFIELDCONTAINER_HH_*/



