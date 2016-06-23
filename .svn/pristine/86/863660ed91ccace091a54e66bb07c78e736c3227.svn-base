//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    LorentzTransformationEM.hh
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
///////////////////////////////////////////////////////////////////////////
#ifndef LASERFIELDORIENTATION_HH_
#define LASERFIELDORIENTATION_HH_


#include <cmath>

class  FieldOrientation
{
public:
	
	FieldOrientation();
	~FieldOrientation();
	
void setCoefficients(double x00, double y00,double z00,
			  		double kx, double ky,double kz,
			  		double mx, double my,double mz);

void OrientCoordinates(double& x, double& y,double& z);
void OrientVector(double& x, double& y,double& z);
void OrientVector2(double& x, double& y,double& z);
	

private:
	
double x0;
double y0;
double z0;
	
double ax;
double bx;
double cx;

double ay;
double by;
double cy;

double az;
double bz;
double cz;

double ax1;
double bx1;
double cx1;

double ay1;
double by1;
double cy1;

double az1;
double bz1;
double cz1;
  


};






#endif /*LASERFIELDORIENTATION_HH_*/
