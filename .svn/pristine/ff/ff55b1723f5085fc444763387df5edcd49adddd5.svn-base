//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.cc
//
// AUTHOR
//    T. Gorlov
//
// CREATED
//    06/28/2008
//
// DESCRIPTION
//    This class provides Lorentz transformations for the electromagnetic field
//    from the laboratory frame to the particle rest frame.
//    mass - mass of the particle in GeV
//    px,py,pz - momentum of the particle the lab frame in GeV/c
//    E_x,E_y,E_z - components of the electric field V/m (parameters are replaced in place) 
//    B_x,B_y,B_z - components of the magnetic field [T] (parameters are replaced in place) 
//
//    OrbitConst::c in [m/sec]
//
///////////////////////////////////////////////////////////////////////////
#include "FieldOrientation.hh"
#include <iostream>


//kx = new z
//ky = new z
//kz = new z

//mx = new x
//my = new x
//mz = new x
FieldOrientation::FieldOrientation()	{

	ax = 1;
	bx = 0;
	cx = 0;
	
	ay = 0;
	by = 1;
	cy = 0;
	
	az = 0;
	bz = 0;
	cz = 1;
	
	ax1 = 1;
	bx1 = 0;
	cx1 = 0;
	
	ay1 = 0;
	by1 = 1;
	cy1 = 0;
	
	az1 = 0;
	bz1 = 0;
	cz1 = 1;
	
	x0 = 0;
	y0 = 0;
	z0 = 0;
	

}

FieldOrientation::~FieldOrientation()	{

}

void FieldOrientation::setCoefficients(double x00, double y00,double z00,
		  								double kx, double ky,double kz,
		  								double mx, double my,double mz)	{


x0 = x00;
y0 = y00;
z0 = z00;

 double a=sqrt((kx*my-ky*mx)*(kx*my-ky*mx)+(kz*mx-kx*mz)*(kz*mx-kx*mz)+(ky*mz-kz*my)*(ky*mz-kz*my));
 double b=sqrt(kx*kx+ky*ky+kz*kz);

 ax = (kz*(kz*mx-kx*mz)-ky*(kx*my-ky*mx))/(a*b);
 bx = (kx*(kx*my-ky*mx)-kz*(ky*mz-kz*my))/(a*b);
 cx = (ky*(ky*mz-kz*my)-kx*(kz*mx-kx*mz))/(a*b);
 
 ay = (ky*mz-kz*my)/a;
 by = (kz*mx-kx*mz)/a;
 cy = (kx*my-ky*mx)/a;
 
 az = kx/b;
 bz = ky/b;
 cz = kz/b;
 
 double d=az*by*cx-ay*bz*cx-az*bx*cy+ax*bz*cy+ay*bx*cz-ax*by*cz;
 
 ax1 = (bz*cy-by*cz)/d;
 bx1 = (bx*cz-bz*cx)/d;
 cx1 = (by*cx-bx*cy)/d;
 
 ay1 = (ay*cz-az*cy)/d;
 by1 = (az*cx-ax*cz)/d;
 cy1 = (ax*cy-ay*cx)/d;
 
 az1 = (az*by-ay*bz)/d;
 bz1 = (ax*bz-az*bx)/d;
 cz1 = (ay*bx-ax*by)/d;
 
return;

}




void FieldOrientation::OrientCoordinates(double& x, double& y,double& z)	{
	

	
double xt=x-x0;
double yt=y-y0;
double zt=z-z0;

x =	ax*xt+bx*yt+cx*zt;
y =	ay*xt+by*yt+cy*zt;
z =	az*xt+bz*yt+cz*zt;



return;

}


void FieldOrientation::OrientVector2(double& x, double& y,double& z)	{
	
double xt=x;
double yt=y;
double zt=z;

x =	ax*xt+bx*yt+cx*zt;
y =	ay*xt+by*yt+cy*zt;
z =	az*xt+bz*yt+cz*zt;

return;

}

void FieldOrientation::OrientVector(double& x, double& y,double& z)	{
	
double xt=x;
double yt=y;
double zt=z;

x =	ax1*xt+bx1*yt+cx1*zt;
y =	ay1*xt+by1*yt+cy1*zt;
z =	az1*xt+bz1*yt+cz1*zt;

return;

}





