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
#include "QuadEMfield.hh"

#include "orbit_mpi.hh"
#include <iostream>

using namespace OrbitUtils;

QuadEMfield::QuadEMfield()
{
	 Ex0v=0;
	 Ey0v=0;
	 Ez0v=0;
	
	 Bx0v=0;
	 By0v=0;
	 Bz0v=0;
	
	 x0v=0;
	 y0v=0;
	 z0v=0;
	
	 cxExv=0;
	 cyExv=0;
	 czExv=0;
	
	 cxEyv=0;
	 cyEyv=0;
	 czEyv=0;
	
	 cxEzv=0;
	 cyEzv=0;
	 czEzv=0;
	
	 cxBxv=0;
	 cyBxv=0;
	 czBxv=0;
	
	 cxByv=0;
	 cyByv=0;
	 czByv=0;
	
	 cxBzv=0;
	 cyBzv=0;
	 czBzv=0;
}


QuadEMfield::~QuadEMfield()
{
}



void QuadEMfield::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& B_x, double& B_y, double& B_z)

{	
	x-=x0v;
	y-=y0v;
	z-=z0v;
	
	 E_x = Ex0v+cxExv*x+cyExv*y+czExv*z; 
	 E_y = Ey0v+cxEyv*x+cyEyv*y+czEyv*z; 
	 E_z = Ez0v+cxEzv*x+cyEzv*y+czEzv*z; 
	 
	 B_x = Bx0v+cxBxv*x+cyBxv*y+czBxv*z; 
	 B_y = By0v+cxByv*x+cyByv*y+czByv*z; 
	 B_z = Bz0v+cxBzv*x+cyBzv*y+czBzv*z;



	 
	 return;

}


void QuadEMfield::Ex0(double val)	{Ex0v=val;}
void QuadEMfield::Ey0(double val)	{Ey0v=val;}
void QuadEMfield::Ez0(double val)	{Ez0v=val;}

void QuadEMfield::Bx0(double val)	{Bx0v=val;}
void QuadEMfield::By0(double val)	{By0v=val;}
void QuadEMfield::Bz0(double val)	{Bz0v=val;}

void QuadEMfield::x0(double val)	{x0v=val;}
void QuadEMfield::y0(double val)	{y0v=val;}
void QuadEMfield::z0(double val)	{z0v=val;}


void QuadEMfield::cxEx(double val)	{cxExv=val;}
void QuadEMfield::cyEx(double val)	{cyExv=val;}
void QuadEMfield::czEx(double val)	{czExv=val;}

void QuadEMfield::cxEy(double val)	{cxEyv=val;}
void QuadEMfield::cyEy(double val)	{cyEyv=val;}
void QuadEMfield::czEy(double val)	{czEyv=val;}

void QuadEMfield::cxEz(double val)	{cxEzv=val;}
void QuadEMfield::cyEz(double val)	{cyEzv=val;}
void QuadEMfield::czEz(double val)	{czEzv=val;}


void QuadEMfield::cxBx(double val)	{cxBxv=val;}
void QuadEMfield::cyBx(double val)	{cyBxv=val;}
void QuadEMfield::czBx(double val)	{czBxv=val;}

void QuadEMfield::cxBy(double val)	{cxByv=val;}
void QuadEMfield::cyBy(double val)	{cyByv=val;}
void QuadEMfield::czBy(double val)	{czByv=val;}

void QuadEMfield::cxBz(double val)	{cxBzv=val;}
void QuadEMfield::cyBz(double val)	{cyBzv=val;}
void QuadEMfield::czBz(double val)	{czBzv=val;}


