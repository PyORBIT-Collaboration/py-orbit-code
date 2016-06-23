//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for field source. It should be sub-classed. The units
//    for E and B are unknown at this point. They should be defined in 
//    subclasses.
//
///////////////////////////////////////////////////////////////////////////
#ifndef BASELASERFIELDSOURCE_HH_
#define BASELASERFIELDSOURCE_HH_

#include "CppPyWrapper.hh"
#include "tcomplex.hh"

#define J tcomplex(0.,1.)


namespace OrbitUtils{
	class  BaseLaserFieldSource: public CppPyWrapper
	{
	public:
		
		BaseLaserFieldSource();
		virtual ~BaseLaserFieldSource();
		
		virtual tcomplex getLaserEMField(double x, double y, double z, double t, 
				tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
				tcomplex& H_x, tcomplex& H_y, tcomplex& H_z);
		
		virtual bool region(double x, double y, double z);
		

		virtual double getFrequencyOmega(double m,double x, double y, double z, double px, double py, double pz, double t);

		
		

	};
};

#endif /*BASELASERFIELDSOURCE_HH_*/
