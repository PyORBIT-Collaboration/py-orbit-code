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
#include <cstdlib>



#include "HminusStripping.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"






#define pop(i,n) PopAttr->attArr(i)[n]		//i-part index, n,m-attr index


using namespace LaserStripping;
using namespace OrbitUtils;








HminusStripping::HminusStripping(int method)
{
	setName("unnamed");

	index = method;
	
	


//allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
k_RungeKutt=new double[5];

	k_RungeKutt[0]=0;

}

HminusStripping::~HminusStripping()
{
	
	delete	[]	k_RungeKutt;

}





void HminusStripping::setupEffects(Bunch* bunch){	
	

                
		prob = new double[bunch->getSize()];
                bunch->setCharge(-1);

		for (int i=0; i<bunch->getSize();i++)	{
                        bunch->recoverParticle(i);
			prob[i] = (double)rand()/(double)RAND_MAX;
		//std::cout<<prob[i]<<"\n";
		}

		

                if(bunch->hasParticleAttributes("Populations")!=0)
                    if(bunch->getParticleAttributes("Populations")->getAttSize() != 1)
                        bunch->removeAllParticleAttributes();


		if(bunch->hasParticleAttributes("Populations")==0)	{
			std::map<std::string,double> part_attr_dict;
			part_attr_dict["size"] = 1;
			bunch->addParticleAttributes("Populations",part_attr_dict);

			for (int i=0; i<bunch->getSize();i++)
			bunch->getParticleAttributes("Populations")->attValue(i,0) = 0.;
		}
	

	PopAttr = bunch->getParticleAttributes("Populations");
	


}
		

	
void HminusStripping::finalizeEffects(Bunch* bunch){
	
	delete	[]	prob;
        

}







void HminusStripping::applyEffectsForEach(Bunch* bunch, int i, 
	                            double* y_in_vct, double* y_out_vct, 
														  double t, double t_step, 														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{
        
        if(bunch->flag(i) == 0)	{
            
            double** partCoordArr = bunch->coordArr();
            
            y_in_vct[0] = partCoordArr[i][0];
            y_in_vct[1] = partCoordArr[i][2];
            y_in_vct[2] = partCoordArr[i][4];
            y_in_vct[3] = partCoordArr[i][1];
            y_in_vct[4] = partCoordArr[i][3];
            y_in_vct[5] = partCoordArr[i][5];            
         
        double mass = bunch->getMass();		        
        double vpt = mass*mass + y_in_vct[3]*y_in_vct[3] + y_in_vct[4]*y_in_vct[4] + y_in_vct[5]*y_in_vct[5];
        
            vpt = 299792458*t_step/sqrt(vpt);
            y_out_vct[0] = y_in_vct[0] + vpt*y_in_vct[3];
            y_out_vct[1] = y_in_vct[1] + vpt*y_in_vct[4];
            y_out_vct[2] = y_in_vct[2] + vpt*y_in_vct[5];
            y_out_vct[3] = y_in_vct[3];
            y_out_vct[4] = y_in_vct[4];
            y_out_vct[5] = y_in_vct[5];
                
          partCoordArr[i][0] = y_out_vct[0];
          partCoordArr[i][2] = y_out_vct[1];
          partCoordArr[i][4] = y_out_vct[2];
          partCoordArr[i][1] = y_out_vct[3];
          partCoordArr[i][3] = y_out_vct[4];
          partCoordArr[i][5] = y_out_vct[5];  
          
       }
    

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
			AmplSolver4step(t_step, i,bunch);
                        
	


//	cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\n";

	
}
















void HminusStripping::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	

	double Ez;
	
		fieldSource->getElectricMagneticField(x0,y0,z0,t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);		
		LorentzTransformationEM::transform(bunch->getMass(),px0,py0,pz0,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);

		
		
	
Ez=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);
	
	Ex_stat=0;	
	Ey_stat=0;		
	Ez_stat=Ez;
	
			
	
}








void	HminusStripping::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{
	
		
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();

//This line calculates relyativistic factor-Gamma
double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
double coeff=1./(gamma*ta);

part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)



Ex_stat/=Ea;	
Ey_stat/=Ea;		
Ez_stat/=Ea;




	
}








double HminusStripping::Gamma(double E_au)	{
	
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	
	double E = E_au*Ea*1e-8;
	double A = 7.96e-14;
	double C = 44.94;
	
	
	return ta/((A/E)*exp(C/E));
}





void HminusStripping::AmplSolver4step(double t_step, int i, Bunch* bunch)	{
	
	


	double dt;
		

		if(bunch->flag(i) == 1)	{
		for(int j=1; j<5; j++)	{	
			
			if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;	
			k_RungeKutt[j]=0;	
			k_RungeKutt[j] += (1 - (pop(i,0)+k_RungeKutt[j-1]*dt))*Gamma(Ez_stat);		
					
		}
		
		pop(i,0) += part_t_step*(k_RungeKutt[1]+2.*k_RungeKutt[2]+2.*k_RungeKutt[3]+k_RungeKutt[4])/6.;	
		}

		
	

		if((index == 1)&&(pop(i,0) >= prob[i]))
			bunch->deleteParticleFast(i);

		
}


