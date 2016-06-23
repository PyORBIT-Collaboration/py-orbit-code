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

#ifndef SNSSTRIPPINGMAGNET_HH_
#define SNSSTRIPPINGMAGNET_HH_

#include "Python.h"
#include "BaseFieldSource.hh"

#include <string>


namespace OrbitUtils{
	
	class  SNSstrippingMagnet: public BaseFieldSource

	{
		public:
		
			SNSstrippingMagnet(std::string data_file);
			~SNSstrippingMagnet();			

			void getElectricMagneticField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);
			
			double getBx(double zi);

			void setA(double A);
			void setPhase(double phase);
                        void setPosition(double pos);
  
		
		  private:
	  
                        double *z;
                        double *Bx;

                        
                        double dz, l, pos, A, phase;
                        int len;				
			
	};

};



#endif /*SNSSTRIPPINGMAGNET_HH_*/



