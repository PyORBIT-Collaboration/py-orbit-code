//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION
//    The base class for C++ implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//     
//         
//       
//     
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>



#include "SchrodingerEquation.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"
#include "FieldRotation.hh"




#define Re(i,m) AmplAttr->attArr(i)[m]		//i-part index, n,m-attr index
#define Im(i,m) AmplAttr->attArr(i)[m+levels]


#define a(i,m) tcomplex(Re(i,m),Im(i,m))


using namespace LaserStripping;
using namespace OrbitUtils;


SchrodingerEquation::SchrodingerEquation(BaseLaserFieldSource*	BaseLaserField, Stark* Stark,double par_res)
{
	setName("unnamed");
	


	StarkEffect=Stark;
	LaserField=BaseLaserField;
	Parameter_resonance=par_res;
	levels = StarkEffect->getStates()*(1+StarkEffect->getStates())*(1+2*StarkEffect->getStates())/6;
	

	




	
//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new tcomplex*[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex[5];
exp_mu_El=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	exp_mu_El[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	exp_mu_El[i][j]=new tcomplex[5];
mu_Elas=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	mu_Elas[i]=new tcomplex*[levels+1];		for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	mu_Elas[i][j]=new tcomplex[3];
cond=new bool*[levels+1];	for (int i=0;i<levels+1;i++)	cond[i]=new bool[levels+1];


for (int i=0; i<levels+1;i++)
	k_RungeKutt[i][0]=tcomplex(0.);

if(LaserField->getPyWrapper() != NULL){
		Py_INCREF(LaserField->getPyWrapper());
	}


if(StarkEffect->getPyWrapper() != NULL){
		Py_INCREF(StarkEffect->getPyWrapper());
	}

}














SchrodingerEquation::~SchrodingerEquation()
{

	for (int i=0;i<levels+1;i++) delete	[]	k_RungeKutt[i];		delete [] k_RungeKutt;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] exp_mu_El[i][j]; for (int i=0;i<levels+1;i++)	delete [] exp_mu_El[i];	delete	[]	exp_mu_El;		
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] mu_Elas[i][j]; for (int i=0;i<levels+1;i++)	delete [] mu_Elas[i];	delete	[]	mu_Elas;

	
	for (int i=0;i<levels+1;i++)	delete	[]	cond[i];		delete	[]	cond;


	

	
	if(LaserField->getPyWrapper() == NULL){
		delete LaserField;
	} else { 
		Py_XDECREF(LaserField->getPyWrapper());
	}
		
	
	if(StarkEffect->getPyWrapper() == NULL){
		delete StarkEffect; 
	} else {
		Py_XDECREF(StarkEffect->getPyWrapper());
	}



}


void SchrodingerEquation::CalcPopulations(int i, Bunch* bunch)	{
		
	PopAttr->attArr(i)[0] = 1;
	for (int j=1; j<levels + 1;j++)	{
	PopAttr->attArr(i)[j] =  Re(i,j)*Re(i,j)+Im(i,j)*Im(i,j);
	PopAttr->attArr(i)[0] -= PopAttr->attArr(i)[j];
	}
	
}




void SchrodingerEquation::setupEffects(Bunch* bunch){	
	
	
	
	
	int_E=new tcomplex**[levels+1];	for (int i=0;i<levels+1;i++)	int_E[i]=new tcomplex*[levels+1]; for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	int_E[i][j]=new tcomplex[bunch->getSize()];

	for(int i=0; i<bunch->getSize();i++)
	for(int n=1; n<levels+1;n++)
	for(int m=1; m<levels+1;m++)
	int_E[n][m][i] = tcomplex(0.,0.);
	

	if(bunch->hasParticleAttributes("Amplitudes")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = 2*levels+1;
		bunch->addParticleAttributes("Amplitudes",part_attr_dict);
		
		for (int i=0; i<bunch->getSize();i++)
		bunch->getParticleAttributes("Amplitudes")->attValue(i,1) = 1;
	}

	
		if(bunch->hasParticleAttributes("Populations")==0)	{
			std::map<std::string,double> part_attr_dict;
			part_attr_dict["size"] = levels+1;
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
		




//-0.00226674 0.000171211 -0.326731 0.000842283-0.000358742 1.69601
	
void SchrodingerEquation::finalizeEffects(Bunch* bunch){
	
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] int_E[i][j]; for (int i=0;i<levels+1;i++)	delete [] int_E[i];	delete	[]	int_E;
	
	
	delete [] nx;
	delete [] ny;
	delete [] nz;
	
	delete [] install_field_dir;
	
}







void SchrodingerEquation::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{
	
x0 = y_in_vct[0];y0 = y_in_vct[1];z0 = y_in_vct[2];
px0 = y_in_vct[3];py0 = y_in_vct[4];pz0 = y_in_vct[5];
x = y_out_vct[0];y = y_out_vct[1];z = y_out_vct[2];
px = y_out_vct[3];py = y_out_vct[4];pz = y_out_vct[5];



					

			GetParticleFrameFields(i, t, t_step,bunch,fieldSource);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetParticleFrameParameters(i,t,t_step,bunch);	
			//	This function provides step solution for density matrix using rk4	method

			AmplSolver4step(i,bunch);	
			
			
			CalcPopulations(i, bunch);		
			


}











void SchrodingerEquation::GetParticleFrameFields(int i,double t, double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	


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
							
		LaserField->getLaserEMField(x0+j*(x-x0)/2,y0+j*(y-y0)/2,z0+j*(z-z0)/2,
				t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
		

			
	LorentzTransformationEM::complex_transform(bunch->getMass(),
																px0,py0,pz0,
																		 Ex_las[j],Ey_las[j],Ez_las[j],
																		 Bx_las[j],By_las[j],Bz_las[j]);	


	FieldRotation::RotateElectricFieldsV(nx[i],ny[i],nz[i],Ex_las[j],Ey_las[j],Ez_las[j]);
	
	}
	}
	else 	
		for(int j=0;j<3;j++)	{		
			Ex_las[j] = 0.;
			Ey_las[j] = 0.;
			Ez_las[j] = 0.;
	}
	
			
}








void	SchrodingerEquation::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
double coeff=1./(gamma*ta);

part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)




for(int j=0;j<3;j++)	
{		
	Ex_las[j]/=Ea; 
	Ey_las[j]/=Ea; 
	Ez_las[j]/=Ea; }

Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;



omega_part=gamma*ta*LaserField->getFrequencyOmega(m,x0,y0,z0,px0,py0,pz0,t);		// frequensy of laser in particle frame (in atomic units)

	
}









void SchrodingerEquation::AmplSolver4step(int i, Bunch* bunch)	{
	


	tcomplex z,z1,z2,JdEdt, mu_x,mu_y,mu_z;
	double dt,dE;
		
	
		StarkEffect->SetE(Ez_stat);	

				

			
			
			for(int n=2; n<levels+1;n++)
			for(int m=1; m<n;m++)	
			if (StarkEffect->field_thresh[n]>Ez_stat && StarkEffect->field_thresh[m]>Ez_stat)	{
				
				StarkEffect->getTransition(n,m,mu_x,mu_y,mu_z);
				dE = fabs(StarkEffect->En[n]-StarkEffect->En[m]);
				JdEdt = tcomplex(0.,dE*part_t_step);
				

				mu_Elas[n][m][1]=(mu_x*conj(Ex_las[1])+mu_y*conj(Ey_las[1])+mu_z*conj(Ez_las[1]))*(J/2.);
				cond[n][m]=fabs(dE - omega_part)<Parameter_resonance*abs(2.*mu_Elas[n][m][1]);

				if (cond[n][m])	{
					
					mu_Elas[n][m][0]=(mu_x*conj(Ex_las[0])+mu_y*conj(Ey_las[0])+mu_z*conj(Ez_las[0]))*(J/2.);	
					mu_Elas[n][m][2]=(mu_x*conj(Ex_las[2])+mu_y*conj(Ey_las[2])+mu_z*conj(Ez_las[2]))*(J/2.);
					
					
					
					z1=exp(int_E[n][m][i]);		
					z2=exp(JdEdt/2.);
				
				exp_mu_El[n][m][1]=z1*mu_Elas[n][m][0];			
				exp_mu_El[n][m][2]=z1*z2*mu_Elas[n][m][1];
				exp_mu_El[n][m][3]=exp_mu_El[n][m][2];
				exp_mu_El[n][m][4]=z1*z2*z2*mu_Elas[n][m][2];
					
					
					
				}
				
				int_E[n][m][i] += JdEdt;
												
			}




		for(int j=1; j<5; j++)
		for(int m=1; m<levels+1; m++)
		if (StarkEffect->field_thresh[m]>Ez_stat)	{

			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;

			k_RungeKutt[m][j]=tcomplex(0);
			for(int k=1;k<m;k++)			if (cond[m][k] && StarkEffect->field_thresh[k]>Ez_stat)		 k_RungeKutt[m][j]-=exp_mu_El[m][k][j]*(a(i,k)+k_RungeKutt[k][j-1]*dt);		
			for(int k=m+1;k<levels+1;k++)	if (cond[k][m] && StarkEffect->field_thresh[k]>Ez_stat) 	 k_RungeKutt[m][j]+=conj(exp_mu_El[k][m][j])*(a(i,k)+k_RungeKutt[k][j-1]*dt); 	
			
			
			k_RungeKutt[m][j]-=StarkEffect->Gamman[m]*(a(i,m)+k_RungeKutt[m][j-1]*dt)/2.;
			

		}


		
		
		
		
		
	for(int m=1;m<levels+1;m++)	{
		
		if (StarkEffect->field_thresh[m]>Ez_stat)	{
		z=part_t_step*(k_RungeKutt[m][1]+2.*k_RungeKutt[m][2]+2.*k_RungeKutt[m][3]+k_RungeKutt[m][4])/6.;	
		Re(i,m)+=z.real();
		Im(i,m)+=z.imag();
		}
		else	{
			Re(i,m) = 0;
			Im(i,m) = 0;
		}
			
	}
	


}


