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



#include "ContinuumSS.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"
#include "FieldRotation.hh"




#define Re(i,m) AmplAttr->attArr(i)[m]		//i-part index, n,m-attr index
#define Im(i,m) AmplAttr->attArr(i)[m+levels+1]


#define a(i,m) tcomplex(Re(i,m),Im(i,m))


using namespace LaserStripping;
using namespace OrbitUtils;


ContinuumSS::ContinuumSS(BaseLaserFieldSource*	BaseLaserField,TDMcontinuum* cont)
{
	

	setName("unnamed");
	

	LaserField = BaseLaserField;

	levels = cont->get_levels();		

		
		dE = new double[levels + 1];
		TDM = new double[levels + 1];

		for (int i=1;i<levels+1;i++)	{
		cont->get_dE_TDM(i,dE[i],TDM[i]);
//		TDM[i] = 60.;
		}



		delE = fabs(dE[2] - dE[1]);
		



	
//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new tcomplex*[levels+1];	for (int i=0;i<levels+1;i++)	k_RungeKutt[i]=new tcomplex[5];
exp_mu_El=new tcomplex*[levels+1];		for (int i=0;i<levels+1;i++)	exp_mu_El[i]=new tcomplex[5]; 
	



for (int i=0; i<levels+1;i++)
	k_RungeKutt[i][0] = 0;

if(LaserField->getPyWrapper() != NULL){
		Py_INCREF(LaserField->getPyWrapper());
	}


}














ContinuumSS::~ContinuumSS()
{

	for (int i=0;i<levels+1;i++) delete	[]	k_RungeKutt[i];		delete [] k_RungeKutt;
	for (int i=0;i<levels+1;i++) delete []	exp_mu_El[i];		delete [] exp_mu_El; 	
	delete [] TDM;
	delete [] dE;
	

	
	if(LaserField->getPyWrapper() == NULL){
		delete LaserField;
	} else { 
		Py_XDECREF(LaserField->getPyWrapper());
	}

}


void ContinuumSS::CalcPopulations(int i, Bunch* bunch)	{
		

	for (int m=0; m<levels + 1;m++)	{
	PopAttr->attArr(i)[m] = pow(abs(a(i,m)),2);
	}


	
/*
	double sum = 0;
	for (int m=1; m<levels + 1;m++){
		sum += pow(abs(a(i,m)),2);
	}
		
		cout.precision(30);
	cout<<pow(abs(a(i,0)),2) + sum*delE<<endl;
*/
//	cout<<pow(abs(a(i,0)),2)<<endl;
}




void ContinuumSS::setupEffects(Bunch* bunch)	{	
	
	
	int_E=new double*[levels+1];	for (int i=0;i<levels+1;i++)	int_E[i]=new double[bunch->getSize()];

	for(int i=0; i<bunch->getSize();i++)
	for(int n=1; n<levels+1;n++)
	int_E[n][i] = 0;
	

	if(bunch->hasParticleAttributes("Amplitudes")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = 2*levels+2;
		bunch->addParticleAttributes("Amplitudes",part_attr_dict);
		
		for (int i=0; i<bunch->getSize();i++)
		bunch->getParticleAttributes("Amplitudes")->attValue(i,0) = 1;
	}

	
		if(bunch->hasParticleAttributes("Populations")==0)	{
			std::map<std::string,double> part_attr_dict;
			part_attr_dict["size"] = levels+1;
			bunch->addParticleAttributes("Populations",part_attr_dict);
		}
	
	

	AmplAttr = bunch->getParticleAttributes("Amplitudes");
	PopAttr = bunch->getParticleAttributes("Populations");
	
	

	
	for (int i=0; i<bunch->getSize();i++)
	CalcPopulations(i, bunch);

}
		





	
void ContinuumSS::finalizeEffects(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)
	CalcPopulations(i, bunch);
	
	for (int i=0;i<levels+1;i++) delete []	int_E[i];			delete [] int_E; 
	

	
}







void ContinuumSS::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{
	
x0 = y_in_vct[0];y0 = y_in_vct[1];z0 = y_in_vct[2];
px0 = y_in_vct[3];py0 = y_in_vct[4];pz0 = y_in_vct[5];
x = y_out_vct[0];y = y_out_vct[1];z = y_out_vct[2];
px = y_out_vct[3];py = y_out_vct[4];pz = y_out_vct[5];

//if(int(t/t_step)%100 == 0)
//cout<<int(t/t_step)<<endl;
					
if(LaserField->region(x,y,z))	{
			GetParticleFrameFields(i, t, t_step,bunch,fieldSource);
	
			//	This function gives parameters Ez_stat	Ex_las[1...3]	Ey_las[1...3]	Ez_las[1...3]	t_part	omega_part	part_t_step 
			//in atomic units in frame of particle		
			GetParticleFrameParameters(i,t,t_step,bunch);	
			//	This function provides step solution for density matrix using rk4	method

			AmplSolver4step(i,bunch);	

			
//			CalcPopulations(i, bunch);	
			PopAttr->attArr(i)[0] = pow(abs(a(i,0)),2);


}

}











void ContinuumSS::GetParticleFrameFields(int i,double t, double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	
	double E_x,E_y,E_z,absE;

		
	for (int j=0; j<3;j++)	{
							
		exp_phasa[j] = LaserField->getLaserEMField(x0+j*(x-x0)/2,y0+j*(y-y0)/2,z0+j*(z-z0)/2,
				t+j*t_step/2,Ex_las[j],Ey_las[j],Ez_las[j],Bx_las[j],By_las[j],Bz_las[j]);
		

			
	LorentzTransformationEM::complex_transform(bunch->getMass(),
																px0,py0,pz0,
																		 Ex_las[j],Ey_las[j],Ez_las[j],
																		 Bx_las[j],By_las[j],Bz_las[j]);	
	

	
	if(j==0)	absE = sqrt(pow(abs(Ex_las[j]),2)+pow(abs(Ey_las[j]),2)+pow(abs(Ez_las[j]),2));

		
			
			E_las[j] = absE*exp_phasa[j];	

	
	}

	
	


	
}








void	ContinuumSS::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
double coeff=1./(gamma*ta);

part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)


for(int j=0;j<3;j++)			
	E_las[j]/=Ea; 

//cout.precision(20);
//omega = gamma*ta*LaserField->getFrequencyOmega(m,x0,y0,z0,px0,py0,pz0,t);		// frequensy of laser in particle frame (in atomic units)
//time_part = t*coeff;
//cout<<omega<<endl;

}






void ContinuumSS::AmplSolver4step(int i, Bunch* bunch)	{
	


	tcomplex z,z1,z2,integral;
	double dt,JdEdt;


	
			for(int n=1; n<levels + 1; n++)	{
				

				JdEdt = dE[n]*part_t_step;		

					z1 = tcomplex(cos(int_E[n][i]),-sin(int_E[n][i]));
					z2 = tcomplex(cos(JdEdt/2),-sin(JdEdt/2));
				
				exp_mu_El[n][1] = z1*TDM[n]*E_las[0];			
				exp_mu_El[n][2] = z1*z2*TDM[n]*E_las[1];
				exp_mu_El[n][3] = exp_mu_El[n][2];
				exp_mu_El[n][4] = z1*z2*z2*TDM[n]*E_las[2];


				
				int_E[n][i] += JdEdt;
												
			}





			for(int j=1; j<5; j++)	{
				
				if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;
				
				
				
				k_RungeKutt[0][j] = 0;
				for (int k=1;k<levels + 1;k++)	{
					k_RungeKutt[0][j] += exp_mu_El[k][j]*(a(i,k)+k_RungeKutt[k][j-1]*dt)*(J/2.)*delE;
					k_RungeKutt[k][j] = conj(exp_mu_El[k][j])*(a(i,0)+k_RungeKutt[0][j-1]*dt)*(J/2.);
				}

				

				

				    

//				double int_TDM2 = 0;
//				for (int k=1;k<levels + 1;k++)	
//					int_TDM2 += TDM[k]*TDM[k]*delE;			
//				cout<<sqrt(int_TDM2)<<"   pert = "<<(-0.526748971193415637860082304527  +  558.748971193415637860082304526749*0.007*0.007)*1.414213562<<endl;
				
				
				 

			
			}
/*
			tcomplex int1 = 0;
			tcomplex int2 = 0;
			double int3 = 0;
			double E;
			double Er = dE[levels/2]-0.5;
			
		    for (int k=1;k<levels + 1;k++)		{
		    	E = dE[k]-0.5;
		    	int1 += E*a(i,k)*TDM[k]*exp(-J*int_E[k][i])*delE;
		    	int2 += a(i,k)*TDM[k]*exp(-J*int_E[k][i])*delE;
//		    	int3 += TDM[k]*(E - Er)*(E - Er)*delE;
		    	
//		    	int1 += E*a(i,k)*(1/sqrt((E-Er)*(E-Er) + 8.54e-5*8.54e-5/4))*exp(-J*int_E[k][i])*delE;
//		    	int2 += a(i,k)*(1/sqrt((E-Er)*(E-Er) + 8.54e-5*8.54e-5/4))*exp(-J*int_E[k][i])*delE;
		    }
		    cout.precision(20);
		cout<<(int1/int2).real()<<"		"<<-2*(int1/int2).imag()<<endl;
*/		
		
		
		
		
	for(int m=0;m<levels+1;m++)	{
		
		z=part_t_step*(k_RungeKutt[m][1]+2.*k_RungeKutt[m][2]+2.*k_RungeKutt[m][3]+k_RungeKutt[m][4])/6.;	
		Re(i,m)+=z.real();
		Im(i,m)+=z.imag();
			
	}
	


}


