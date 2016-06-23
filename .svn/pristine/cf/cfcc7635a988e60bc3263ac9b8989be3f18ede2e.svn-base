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
#include "ConstEMfield.hh"

#include "orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

ConstEMfield::ConstEMfield()
{
	
	ExConst=0;
	EyConst=0;
	EzConst=0;
	
	BxConst=0;
	ByConst=0;
	BzConst=0;
	
}

ConstEMfield::ConstEMfield(double E_x,double E_y,double E_z,double B_x,double B_y,double B_z)
{
	
	ExConst=E_x;
	EyConst=E_y;
	EzConst=E_z;
	
	BxConst=B_x;
	ByConst=B_y;
	BzConst=B_z;
	

}


ConstEMfield::~ConstEMfield()
{

}

void ConstEMfield::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& B_x, double& B_y, double& B_z)

{
	
	 E_x = ExConst; E_y = EyConst; E_z = EzConst; 
	 B_x = BxConst; B_y = ByConst; B_z = BzConst; 

}


void ConstEMfield::Ex(double val)	{ExConst=val;}
void ConstEMfield::Ey(double val)	{EyConst=val;}
void ConstEMfield::Ez(double val)	{EzConst=val;}

void ConstEMfield::Bx(double val)	{BxConst=val;}
void ConstEMfield::By(double val)	{ByConst=val;}
void ConstEMfield::Bz(double val)	{BzConst=val;}






