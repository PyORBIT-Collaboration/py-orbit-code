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

#include "TwoLevelStrongField.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"
#include "FieldRotation.hh"


#define Re1(i) AmplAttr->attArr(i)[1]		//i-part index, n,m-attr index
#define Im1(i) AmplAttr->attArr(i)[2]
#define Re2(i) AmplAttr->attArr(i)[3]		//i-part index, n,m-attr index
#define Im2(i) AmplAttr->attArr(i)[4]


#define Ampl_1(i) tcomplex(Re1(i),Im1(i))
#define Ampl_2(i) tcomplex(Re2(i),Im2(i))

//#define k_rk(j,n,m) k_RungeKutt[j][(n-1)*levels+m]

using namespace LaserStripping;
using namespace OrbitUtils;


TwoLevelStrongField::TwoLevelStrongField(BaseLaserFieldSource*	BaseLaserField, StarkStrongField* StarkStr)
{
	setName("unnamed");
	
	LaserField = BaseLaserField;
	Stark = StarkStr;

	
	if(LaserField->getPyWrapper() != NULL){
			Py_INCREF(LaserField->getPyWrapper());
	}
	
	if(Stark->getPyWrapper() != NULL){
			Py_INCREF(Stark->getPyWrapper());
	}
	
}

TwoLevelStrongField::~TwoLevelStrongField() 
{
	
			if(LaserField->getPyWrapper() == NULL){
				delete LaserField;
			} else {
				Py_XDECREF(LaserField->getPyWrapper());
			}
			
			if(Stark->getPyWrapper() == NULL){
				delete Stark;
			} else {
				Py_XDECREF(Stark->getPyWrapper());
			}
			
			

}

void TwoLevelStrongField::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[1] =  Re1(i)*Re1(i)+Im1(i)*Im1(i);
	PopAttr->attArr(i)[2] =  Re2(i)*Re2(i)+Im2(i)*Im2(i);
	PopAttr->attArr(i)[0] = 1 - PopAttr->attArr(i)[1]-PopAttr->attArr(i)[2];
	
}

void TwoLevelStrongField::setupEffects(Bunch* bunch){
	
	int_E = new tcomplex[bunch->getSize()];

	
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
	
	
	nx=new double[bunch->getSize()];
	ny=new double[bunch->getSize()];
	nz=new double[bunch->getSize()];
	
	install_field_dir=new bool [bunch->getSize()];
	
	for (int i=0; i<bunch->getSize();i++)
		install_field_dir[i] = true;
	
	
	for (int i=0; i<bunch->getSize();i++)
		CalcPopulations(i, bunch);
}


void TwoLevelStrongField::finalizeEffects(Bunch* bunch)		{
	
	delete	[]	int_E;

	delete [] nx;
	delete [] ny;
	delete [] nz;
	
	delete [] install_field_dir;
	
}


void TwoLevelStrongField::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{
	
	x0 = y_in_vct[0];y0 = y_in_vct[1];z0 = y_in_vct[2];
	px0 = y_in_vct[3];py0 = y_in_vct[4];pz0 = y_in_vct[5];
	x = y_out_vct[0];y = y_out_vct[1];z = y_out_vct[2];
	px = y_out_vct[3];py = y_out_vct[4];pz = y_out_vct[5];
	


		//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	
		//in natural unts (Volt per meter)	in the frame of particle				
		GetParticleFrameFields(i, t,t_step, bunch,fieldSource);

		
		//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
		//in atomic units in frame of particle		
		GetParticleFrameParameters(i,t,t_step,bunch);	
		
		//	This function provides step solution for density matrix using rk4	method
		AmplSolver4step(i,bunch);	
		
		CalcPopulations(i, bunch);
	
}


void TwoLevelStrongField::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch, OrbitUtils::BaseFieldSource* fieldSource)	{
	double E_x,E_y,E_z,Eabs;
	
	
	fieldSource->getElectricMagneticField(x0,y0,z0,t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);	
	LorentzTransformationEM::transform(bunch->getMass(),px0,py0,pz0,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);	
	
	if (install_field_dir[i])	{
	nx[i] = Ex_stat;	
	ny[i] = Ey_stat;	
	nz[i] = Ez_stat;
	install_field_dir[i]=false;
	}	
	
	
	Ez_stat=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);
	Ex_stat=0;	
	Ey_stat=0;
	
	
	
	if(LaserField->region(x,y,z))	{
	
	for (int j=0; j<3;j++)	{
		
		exp_phasa[j] = LaserField->getLaserEMField(x0+j*(x-x0)/2,y0+j*(y-y0)/2,z0+j*(z-z0)/2,
			t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
		
		LorentzTransformationEM::complex_electric_transform(bunch->getMass(),
			px0,py0,pz0,
			Ex_las[j],Ey_las[j],Ez_las[j],
			Bx_las[j],By_las[j],Bz_las[j]);	
		
		
		FieldRotation::RotateElectricFieldsV(nx[i],ny[i],nz[i],Ex_las[j],Ey_las[j],Ez_las[j]);
		
		if(j==0)		{
			E_x=abs(Ex_las[j]);
			E_y=abs(Ey_las[j]);
			E_z=abs(Ez_las[j]);
			
			
			
			if(Stark->m == 0)
			Eabs = E_z;
			if(Stark->m == 1)
			Eabs=sqrt(pow(E_x,2)+pow(E_y,2));
		}
		
		E_las[j] = Eabs*exp_phasa[j];	
	}
	}
	else 	
		for(int j=0;j<3;j++)	{		
		E_las[j] = 0.; 
	}
		
	
//	if(((int)(t/t_step))%1000 == 0)
//	cout<<t<<"\t"<<Eabs/5.14220642e+011<<"\n";
}


void	TwoLevelStrongField::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e+011;				//Atomic unit of electric field
	double m=bunch->getMass();
	
	//This line calculates relyativistic factor-Gamma
	double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
	double coeff=1./(gamma*ta);
	part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)
	
	//Convertion form SI units to atomic units
	for(int j=0;j<3;j++)	{		
		E_las[j]/=Ea; 
	}
	
	Ex_stat/=Ea;	
	Ey_stat/=Ea;		
	Ez_stat/=Ea;
	
//		if(((int)(t/t_step))%100 == 0)
//		cout<<t<<"\t"<<Ez_stat<<"\n";
	
//	double omega = gamma*ta*LaserField->getFrequencyOmega(m,x0,y0,z0,px0,py0,pz0,t);	
	
	

}

void TwoLevelStrongField::AmplSolver4step(int i, Bunch* bunch)	{

	tcomplex z1,z2,z,JdEdt;
	double dt;
	
	Stark->SetE(Ez_stat);

	
	JdEdt = tcomplex(0.,Stark->deltaEn*part_t_step);
				
	
	
	
	for(int j=0; j<3;j++)	mu_Elas[j]=Stark->dipole_transitionn*conj(E_las[j])*J/2.;	
	

			
			z1=exp(int_E[i]);		
			z2=exp(JdEdt/2.);

			exp_mu_El[1]=z1*mu_Elas[0];			
			exp_mu_El[2]=z1*z2*mu_Elas[1];
			exp_mu_El[3]=exp_mu_El[2];
			exp_mu_El[4]=z1*z2*z2*mu_Elas[2];
			
			int_E[i] += JdEdt;

			
			
			for(int j=1; j<5; j++)	{	
				if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

				k_RungeKutt_1[j]=conj(exp_mu_El[j])*(Ampl_2(i)+k_RungeKutt_2[j-1]*dt);	
				k_RungeKutt_2[j]=-exp_mu_El[j]*(Ampl_1(i)+k_RungeKutt_1[j-1]*dt);	
				
				k_RungeKutt_2[j]-=Stark->Gamman*(Ampl_2(i)+k_RungeKutt_2[j-1]*dt)/2.;
			}
				
			
		z1=(k_RungeKutt_1[1]+2.*k_RungeKutt_1[2]+2.*k_RungeKutt_1[3]+k_RungeKutt_1[4])/6.;	
		z2=(k_RungeKutt_2[1]+2.*k_RungeKutt_2[2]+2.*k_RungeKutt_2[3]+k_RungeKutt_2[4])/6.;
		
		Re1(i)+=z1.real()*part_t_step;
		Im1(i)+=z1.imag()*part_t_step;
		
		Re2(i)+=z2.real()*part_t_step;
		Im2(i)+=z2.imag()*part_t_step;
}
