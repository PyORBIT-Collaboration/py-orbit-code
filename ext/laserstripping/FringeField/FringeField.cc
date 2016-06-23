//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppBaseFieldSource.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implements 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//
///////////////////////////////////////////////////////////////////////////
#include "FringeField.hh"

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cmath>
#include <fstream>
#include <iostream>

using namespace OrbitUtils;

FringeField::FringeField(double d,double F, int derivative)
{
	
	gap = d;
	Bx = F;
	der = derivative;

}




FringeField::~FringeField()
{ 


}







void FringeField::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{
	double coeff = -der*z/gap;
	
	double coeff_2 = coeff;
	double coeff_3 = coeff*coeff_2;
	double coeff_4 = coeff*coeff_3;
	double coeff_5 = coeff*coeff_4;
	double coeff_6 = coeff*coeff_5;
	
	
	E_x=0;
	E_y=0;
	E_z=0;
	
	H_y=0;
	H_z=0;
	
	H_x=Bx/(1+exp(0.478959+1.911289*coeff_2-1.185953*coeff_3+1.630554*coeff_4-1.082657*coeff_5+0.318111*coeff_6));


return;

}




