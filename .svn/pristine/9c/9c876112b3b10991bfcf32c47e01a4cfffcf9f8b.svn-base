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
#include "ConeLFmode.hh"
#include "MathPolynomial.hh"

#include "orbit_mpi.hh"
#include "OrbitConst.hh"
#include <iostream>



using namespace OrbitUtils;

ConeLFmode::ConeLFmode(double P,double r_x,double r_y,double a_x,double a_y,double la)
{
	U=2*sqrt(P*OrbitConst::c*OrbitConst::permeability/OrbitConst::PI);

	orient = new FieldOrientation();  
	Laser_lambda=la;
	k=2*OrbitConst::PI/Laser_lambda;
	rx=r_x;
	ry=r_y;
	ax=a_x;
	ay=a_y;

	

	//default orientation of laser field polarizations
	nEx=1;nEy=0;nEz=0;
	nHx=0;nHy=1;nHz=0;


	
	
	
}

ConeLFmode::~ConeLFmode()
{	
}









void ConeLFmode::setLaserFieldOrientation(double x0, double y0, double z0,
													double kx, double ky, double kz,
													double mx, double my, double mz,
													double n_Ex, double n_Ey, double n_Ez)	
{
	double a=sqrt(n_Ex*n_Ex+n_Ey*n_Ey+n_Ez*n_Ez);
	orient->setCoefficients(x0, y0, z0,kx, ky, kz, mx, my, mz);

	
	nEx=n_Ex/a;	
	nEy=n_Ey/a;	
	nEz=n_Ez/a;	//Vector of polarization of electric field
	
double	n_Hx=nEz*ky-nEy*kz;
double	n_Hy=nEx*kz-nEz*kx;
double	n_Hz=nEy*kx-nEx*ky;
	
double b=sqrt(n_Hx*n_Hx+n_Hy*n_Hy+n_Hz*n_Hz);
	
	nHx=n_Hx/b;
	nHy=n_Hy/b;
	nHz=n_Hz/b;
	

}









double ConeLFmode::HermiteGaussianOmega(double m,double x, double y, double z, double px, double py, double pz, double t){
	
	
	double a=OrbitConst::c/sqrt(m*m+px*px+py*py+pz*pz);
	
	double vx=px*a;
	double vy=py*a;
	double vz=pz*a;


	return 	k*(OrbitConst::c-vz-x*vx*ax/(z*ax+rx)-y*vy*ay/(z*ay+rx));

}












void ConeLFmode::getLaserElectricMagneticField(double x, double y, double z, double t, 
		tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
		tcomplex& H_x, tcomplex& H_y, tcomplex& H_z){
	


			orient->OrientCoordinates(x,y,z);
tcomplex	E=U*getNonOrientedU(x,y,z,t);


tcomplex	H=E/OrbitConst::c;



	E_x=E*nEx;	E_y=E*nEy;	E_z=E*nEz;
	H_x=H*nHx;	H_y=H*nHy;	H_z=H*nHz;
	
}







double ConeLFmode::getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t){

	
	orient->OrientCoordinates(x,y,z);
	orient->OrientVector2(px,py,pz);
	

	return HermiteGaussianOmega(m,x,y,z,px,py,pz,t);

}











tcomplex ConeLFmode::getNonOrientedU(double x, double y, double z, double t){

	
	
	double rxax= 1/(rx+ax*z);
	double ryay= 1/(ry+ay*z);
		

		tcomplex p=tcomplex(0.,k*(z-t*OrbitConst::c+x*x*ax*rxax/2+y*y*ay*ryay/2));

		
		return	sqrt(rxax*ryay)*exp(-pow(x*rxax,2)-pow(y*ryay,2)-p);
	

	
}










