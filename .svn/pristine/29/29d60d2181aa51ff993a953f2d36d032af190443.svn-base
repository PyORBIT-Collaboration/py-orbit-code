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
#include "HermiteGaussianLFmode.hh"
#include "MathPolynomial.hh"

#include "orbit_mpi.hh"
#include "OrbitConst.hh"
#include <iostream>


#define J tcomplex(0.,1.)

using namespace OrbitUtils;

HermiteGaussianLFmode::HermiteGaussianLFmode(tcomplex C,int n,int m,double w_x,double w_y,double f_x,double f_y,double la,double env__peak, double env__sigma)
{
	Cnm=C;
	Unm=Cnm*sqrt(2*OrbitConst::c*OrbitConst::permeability)*sqrt(pow(w_x,2*n+1)*pow(w_y,2*m+1)*MathPolynomial::Factorial(n)*MathPolynomial::Factorial(m)*pow(2,n+m+1)/(MathPolynomial::Factorial(2*n)*MathPolynomial::Factorial(2*m)*OrbitConst::PI));

	orient = new FieldOrientation();  
	Laser_lambda=la;
	env_peak = env__peak;
	env_sigma = env__sigma;
	fx=f_x;
	fy=f_y;
	wx=w_x;
	wy=w_y;
	n_moda=n;
	m_moda=m;
	k=2*OrbitConst::PI/Laser_lambda;
	
	la_pi_wx = Laser_lambda/(OrbitConst::PI*wx);
	la_pi_wy = Laser_lambda/(OrbitConst::PI*wy);

	//default orientation of laser field polarizations
	nEx=1;nEy=0;nEz=0;
	nHx=0;nHy=1;nHz=0;


	
	
	
}

HermiteGaussianLFmode::~HermiteGaussianLFmode()
{	
}









void HermiteGaussianLFmode::setLaserFieldOrientation(double x0, double y0, double z0,
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




void HermiteGaussianLFmode::setLocalParameters(double rx, double ry, double ax, double ay)	{
	
	fx=-ax*rx*rx*rx/(4/(k*k)+ax*ax*rx*rx);
	fy=-ay*ry*ry*ry/(4/(k*k)+ay*ay*ry*ry);
	wx=2*rx/(k*sqrt(4/(k*k)+ax*ax*rx*rx));
	wy=2*ry/(k*sqrt(4/(k*k)+ay*ay*ry*ry));
	
	int n=n_moda;
	int m=m_moda;
	
	Unm=Cnm*sqrt(2*OrbitConst::c*OrbitConst::permeability)*sqrt(pow(wx,2*n+1)*pow(wy,2*m+1)*MathPolynomial::Factorial(n)*MathPolynomial::Factorial(m)*pow(2,n+m+1)/(MathPolynomial::Factorial(2*n)*MathPolynomial::Factorial(2*m)*OrbitConst::PI));


}




double HermiteGaussianLFmode::HermiteGaussianOmega(double m,double x, double y, double z, double px, double py, double pz, double t){
	
	
	
	double a=OrbitConst::c/sqrt(m*m+px*px+py*py+pz*pz);
	
	double vx=px*a;
	double vy=py*a;
	double vz=pz*a;
	


double ax=k*k*pow(wx,4)+4*pow(fx-z,2);
double ay=k*k*pow(wy,4)+4*pow(fy-z,2);

double bx = 4*x*(fx-z);
double by = 4*y*(fy-z);

	return 	k*(OrbitConst::c+vz*(pow(bx/ax,2)+pow(by/ay,2)-1) + (vz*(wx*wx-2*x*x) + bx*vx)/ax + (vz*(wy*wy-2*y*y) + by*vy)/ay);

}












tcomplex HermiteGaussianLFmode::getLaserEMField(double x, double y, double z, double t, 
		tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
		tcomplex& H_x, tcomplex& H_y, tcomplex& H_z){
	
	tcomplex exp_phasa;

			orient->OrientCoordinates(x,y,z);
			
tcomplex	E=Unm*getNonOrientedU(n_moda,m_moda,x,y,z,t);

			exp_phasa = E/abs(E);
			double absp = abs(exp_phasa); 
			if ((absp != absp)||(fabs(absp-1.)>1.0e-3))
			exp_phasa = tcomplex(0.,0.);

			if(env_sigma>0)	{
			double ex = ((z - env_peak)/OrbitConst::c - t)/env_sigma;
			E *= exp(-ex*ex/4);
			}

tcomplex	H=E/OrbitConst::c;



	E_x=E*nEx;	E_y=E*nEy;	E_z=E*nEz;
	H_x=H*nHx;	H_y=H*nHy;	H_z=H*nHz;
		
	return exp_phasa;
	
}


bool HermiteGaussianLFmode::region(double x, double y, double z){
	
	double ax = (z-fx)*la_pi_wx;
	double ay = (z-fy)*la_pi_wy;
	
	
	orient->OrientCoordinates(x,y,z);
	
	return (x*x/(wx*wx+ax*ax) + y*y/(wy*wy+ay*ay) < 25.);
	
}




double HermiteGaussianLFmode::getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t){

	
	orient->OrientCoordinates(x,y,z);
	orient->OrientVector2(px,py,pz);
	

	return HermiteGaussianOmega(m,x,y,z,px,py,pz,t);

}











tcomplex HermiteGaussianLFmode::getNonOrientedU(int n, int m, double x, double y, double z, double t){

	tcomplex U;
	
	if (n==0&&m==0)	{
		
		tcomplex a=tcomplex(0.,k*(z-t*OrbitConst::c));
		tcomplex funx=pow(tcomplex(wx*wx,-(z-fx)*2/k),-1);
		tcomplex funy=pow(tcomplex(wy*wy,-(z-fy)*2/k),-1);
		
		U = sqrt(funx*funy)*exp(-x*x*funx-y*y*funy-a);	
		if (abs(U) != abs(U)) U = tcomplex(0.,0.);
		
		return	U;	
	
	}	else	{
	
	
	
	tcomplex funx=sqrt(tcomplex(wx*wx-2.0*J*(z-fx)/k));
	tcomplex funy=sqrt(tcomplex(wy*wy-2.0*J*(z-fy)/k));
	tcomplex xf=x/funx;
	tcomplex yf=y/funy;
	
	U = pow(funx,-n-1)*pow(funy,-m-1)*exp(-xf*xf-yf*yf-J*k*(z-t*OrbitConst::c))*MathPolynomial::ComplexHermite(n,xf)*MathPolynomial::ComplexHermite(m,yf);
	if (abs(U) != abs(U)) U = tcomplex(0.,0.);
	
	return	U;
	
	}
	

	
}










