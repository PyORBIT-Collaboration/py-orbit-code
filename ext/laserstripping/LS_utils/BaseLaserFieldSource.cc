//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    BaseFieldSource.cc
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


#include "BaseLaserFieldSource.hh"

using namespace OrbitUtils;

BaseLaserFieldSource::BaseLaserFieldSource(): CppPyWrapper(NULL)
{
}

BaseLaserFieldSource::~BaseLaserFieldSource()
{
}

tcomplex BaseLaserFieldSource::getLaserEMField(double x, double y, double z, double t, 
		tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
		tcomplex& H_x, tcomplex& H_y, tcomplex& H_z)
{	
  E_x = 0.0; E_y = 0.0; E_z = 0.0; 
  H_x = 0.0; H_y = 0.0; H_z = 0.0; 
  
  return tcomplex(0.,0.);
}

bool BaseLaserFieldSource::region(double x, double y, double z){
	
	return true;
}


double BaseLaserFieldSource::getFrequencyOmega(double m,double x, double y, double z,double px, double py, double pz, double t)
{	
		return 0.;
}



