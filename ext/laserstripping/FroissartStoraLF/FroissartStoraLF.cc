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
#include "FroissartStoraLF.hh"
#include "orbit_mpi.hh"
#include <iostream>

#define J tcomplex(0.,1.)

using namespace OrbitUtils;

FroissartStoraLF::FroissartStoraLF(double Om,double G,double El)
{
	
	Omega=Om;
	Gamma=G;
	amp_ELas = El;
	

	//default polarization of laser field
	nEx=1;
	nEy=0;
	nEz=0;
	
}



FroissartStoraLF::~FroissartStoraLF()
{
}









void FroissartStoraLF::setLaserFieldPolarization(double n_Ex, double n_Ey, double n_Ez)		{
	
	double a=sqrt(n_Ex*n_Ex+n_Ey*n_Ey+n_Ez*n_Ez);

	
	nEx=n_Ex/a;	
	nEy=n_Ey/a;	
	nEz=n_Ez/a;	//Vector of polarization of electric field
		

}







tcomplex FroissartStoraLF::getNonOrientedU(double t){
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e11;				//Atomic unit of electric field
	t/=ta;									//convertion of natural time to atomic unit time
	

	return amp_ELas*Ea*exp(J*(Omega*t+Gamma*t*t/2));
	

}








tcomplex FroissartStoraLF::getLaserEMField(double x, double y, double z, double t, 
		tcomplex& E_x, tcomplex& E_y, tcomplex& E_z,
		tcomplex& H_x, tcomplex& H_y, tcomplex& H_z){
	
	
tcomplex	E=getNonOrientedU(t);


	E_x=E*nEx;	E_y=E*nEy;	E_z=E*nEz;
	H_x=0;	H_y=0;	H_z=0;
	
	return E/abs(E);
	
}






double FroissartStoraLF::getFrequencyOmega(double m, double x, double y, double z, double px, double py, double pz, double t){

	double ta=2.418884326505e-17;			//atomic unit of time
	t/=ta;									//convertion of natural time to atomic unit time
	
	return (Omega+Gamma*t)/ta;
	
}











