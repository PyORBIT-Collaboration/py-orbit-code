//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION    
//     
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "TwoLevelAtom.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"


#define Re1(i) AmplAttr->attArr(i)[1]		//i-part index, n,m-attr index
#define Im1(i) AmplAttr->attArr(i)[2]
#define Re2(i) AmplAttr->attArr(i)[3]		//i-part index, n,m-attr index
#define Im2(i) AmplAttr->attArr(i)[4]


#define Ampl_1(i) tcomplex(Re1(i),Im1(i))
#define Ampl_2(i) tcomplex(Re2(i),Im2(i))

//#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]

using namespace LaserStripping;
using namespace OrbitUtils;


TwoLevelAtom::TwoLevelAtom(BaseLaserFieldSource*	BaseLaserField, double delta_E, double dipole_tr)
{
	setName("unnamed");
	
	LaserField=BaseLaserField;
	d_Energy=delta_E;
	dip_transition=dipole_tr;
	
	if(LaserField->getPyWrapper() != NULL){
			Py_INCREF(LaserField->getPyWrapper());
	}
}

TwoLevelAtom::~TwoLevelAtom() 
{
	
			if(LaserField->getPyWrapper() == NULL){
				delete LaserField;
			} else {
				Py_XDECREF(LaserField->getPyWrapper());
			}

}

void TwoLevelAtom::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[1] =  Re1(i)*Re1(i)+Im1(i)*Im1(i);
	PopAttr->attArr(i)[2] =  Re2(i)*Re2(i)+Im2(i)*Im2(i);
	PopAttr->attArr(i)[0] = 1 - PopAttr->attArr(i)[1]-PopAttr->attArr(i)[2];
	
}

void TwoLevelAtom::setupEffects(Bunch* bunch){
	
	
	if(bunch->hasParticleAttributes("Amplitudes")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = 5;
		bunch->addParticleAttributes("Amplitudes",part_attr_dict);
		
		for (int i=0; i<bunch->getSize();i++){
			bunch->getParticleAttributes("Amplitudes")->attValue(i,1) = 1;
		}
	}
	
	
	if(bunch->hasParticleAttributes("Populations")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = 3;
		bunch->addParticleAttributes("Populations",part_attr_dict);
	}
	
	

	AmplAttr = bunch->getParticleAttributes("Amplitudes");
	PopAttr = bunch->getParticleAttributes("Populations");
	
	for (int i=0; i<bunch->getSize();i++)
		CalcPopulations(i, bunch);
}



void TwoLevelAtom::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{
	
	x0 = y_in_vct[0];y0 = y_in_vct[1];z0 = y_in_vct[2];
	px0 = y_in_vct[3];py0 = y_in_vct[4];pz0 = y_in_vct[5];
	x = y_out_vct[0];y = y_out_vct[1];z = y_out_vct[2];
	px = y_out_vct[3];py = y_out_vct[4];pz = y_out_vct[5];
	


		if(LaserField->region(x,y,z))	{
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
			//in natural unts (Volt per meter)	in the frame of particle				
		GetParticleFrameFields(i, t,t_step, bunch);

		
		//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
		//in atomic units in frame of particle		
		GetParticleFrameParameters(i,t,t_step,bunch);	
		
		//	This function provides step solution for density matrix using rk4	method
		AmplSolver4step(i,bunch);	
		
		CalcPopulations(i, bunch);
		
		}
	
}


void TwoLevelAtom::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch)	{
	double E_x,E_y,E_z,nx,ny,nz,Eabs;
	
	for (int j=0; j<3;j++)	{
		
		exp_phasa[j] = LaserField->getLaserEMField(x0+j*(x-x0)/2,y0+j*(y-y0)/2,z0+j*(z-z0)/2,
			t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
		
		LorentzTransformationEM::complex_electric_transform(bunch->getMass(),
			px0,py0,pz0,
			Ex_las[j],Ey_las[j],Ez_las[j],
			Bx_las[j],By_las[j],Bz_las[j]);	
		
		if(j==0)		{
			E_x=abs(Ex_las[j]);
			E_y=abs(Ey_las[j]);
			E_z=abs(Ez_las[j]);
			Eabs=sqrt(pow(E_x,2)+pow(E_y,2)+pow(E_z,2));
		}
		
		Ez_las[j] = Eabs*exp_phasa[j];	
	}
//	cout<<t<<"\t"<<Eabs<<"\n";
}


void	TwoLevelAtom::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e+011;				//Atomic unit of electric field
	double m=bunch->getMass();
	
	//This line calculates relyativistic factor-Gamma
	double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
	double coeff=1./(gamma*ta);
	part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)
	t_part=t*coeff;	 
	
	//Convertion form SI units to atomic units
	for(int j=0;j<3;j++)	{		
		Ez_las[j]/=Ea; 
		
	}
	//cout<<gamma*ta*LaserField->getFrequencyOmega(m,x0(i),y0(i),z0(i),px0(i),py0(i),pz0(i),t)<<"\n";	
}

void TwoLevelAtom::AmplSolver4step(int i, Bunch* bunch)	{

	tcomplex z1,z2;
	double dt;
				
	for(int j=0; j<3;j++)	mu_Elas[j]=dip_transition*conj(Ez_las[j])*J/2.;		
			z1=exp(J*t_part*fabs(d_Energy));
			z2=exp(J*part_t_step*fabs(d_Energy)/2.);

			exp_mu_El[1]=z1*mu_Elas[0];			
			exp_mu_El[2]=z1*z2*mu_Elas[1];
			exp_mu_El[3]=exp_mu_El[2];
			exp_mu_El[4]=z1*z2*z2*mu_Elas[2];
			
			for(int j=1; j<5; j++)	{	
				if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

				k_RungeKutt_1[j]=conj(exp_mu_El[j])*(Ampl_2(i)+k_RungeKutt_2[j-1]*dt);	
				k_RungeKutt_2[j]=-exp_mu_El[j]*(Ampl_1(i)+k_RungeKutt_1[j-1]*dt);	
			}
				
		z1=(k_RungeKutt_1[1]+2.*k_RungeKutt_1[2]+2.*k_RungeKutt_1[3]+k_RungeKutt_1[4])/6.;	
		z2=(k_RungeKutt_2[1]+2.*k_RungeKutt_2[2]+2.*k_RungeKutt_2[3]+k_RungeKutt_2[4])/6.;
		
		Re1(i)+=z1.real()*part_t_step;
		Im1(i)+=z1.imag()*part_t_step;
		
		Re2(i)+=z2.real()*part_t_step;
		Im2(i)+=z2.imag()*part_t_step;
}
