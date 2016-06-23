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
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "DM_noLaserField.hh"
#include "RungeKuttaTracker.hh"
#include "OrbitConst.hh"
#include "LorentzTransformationEM.hh"

#define pop(i,n) PopAttr->attArr(i)[n]		//i-part index, n,m-attr index


using namespace LaserStripping;
using namespace OrbitUtils;

inline int convert3to1level(int n,int n1, int m){
	return 1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
}

DM_noLaserField::DM_noLaserField(Stark* Starkef, int method)
{
	setName("unnamed");
        index = method;

	StarkEffect=Starkef;
	
	int st = StarkEffect->getStates();
	levels=	st*(1+st)*(1+2*st)/6;
	
  //allocating memory for koefficients of 4-th order Runge-Kutta method and other koeeficients of the master equation
	
	k_RungeKutt=new double*[levels+1]; for (int i=0;i<levels+1;i++) k_RungeKutt[i]=new double[5];
	gamma_ij=new double*[levels+1];	for (int i=0;i<levels+1;i++)	gamma_ij[i]=new double[levels+1];
	cond=new bool*[levels+1];	for (int i=0;i<levels+1;i++)	cond[i]=new bool[levels+1];
	
	for(int n=1;n<st+1;n++){
		for(int m=-(n-1);m<(n-1)+1;m++){
			for(int n1=0;n1<n-abs(m)-1+1;n1++){
				
				for(int ns=1;ns<st+1;ns++){
					for(int ms=-(ns-1);ms<(ns-1)+1;ms++){
						for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++){
							
							cond[convert3to1level(n,n1,m)][convert3to1level(ns,n1s,ms)]=(n!=ns);
							
						}
					}
				}
			}
		}
	}
	
	for (int i=0; i<levels+1;i++)
		k_RungeKutt[i][0]=0;
	
	
	if(StarkEffect->getPyWrapper() != NULL){
		Py_INCREF(StarkEffect->getPyWrapper());
	}
	
}

DM_noLaserField::~DM_noLaserField()
{
	for (int i=0;i<levels+1;i++) 	delete [] k_RungeKutt[i]; 	delete	[]	k_RungeKutt;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_ij[i];	delete	[]	gamma_ij;
	for (int i=0;i<levels+1;i++)	delete	[]	cond[i];		delete	[]	cond;
	
	if(StarkEffect->getPyWrapper() == NULL){
		delete StarkEffect;
	} else {
		Py_XDECREF(StarkEffect->getPyWrapper());
	}
}

void DM_noLaserField::CalcPopulations(int i, Bunch* bunch)	{
	PopAttr->attArr(i)[0] = 1;
	for (int j=1; j<levels + 1;j++)		
	PopAttr->attArr(i)[0] -= PopAttr->attArr(i)[j];
        
}




void DM_noLaserField::setupEffects(Bunch* bunch){	
        

        prob = new double[bunch->getSize()];
        bunch->setCharge(+1);
        for (int i=0; i<bunch->getSize();i++)	{
                bunch->deleteParticleFast(i);
                prob[i] = (double)rand()/(double)RAND_MAX;
//		std::cout<<prob[i]<<"\n";
        }
                
        
        if(bunch->hasParticleAttributes("Populations")!=0)
            if(bunch->getParticleAttributes("Populations")->getAttSize() != levels+1)
                bunch->removeAllParticleAttributes();
        
	
	if(bunch->hasParticleAttributes("Populations")==0)	{
		std::map<std::string,double> part_attr_dict;
		part_attr_dict["size"] = levels+1;
		bunch->addParticleAttributes("Populations",part_attr_dict);
		
		for (int i=0; i<bunch->getSize();i++)
			bunch->getParticleAttributes("Populations")->attValue(i,1) = 1;
	}
	
	
	PopAttr = bunch->getParticleAttributes("Populations");
	
	for (int i=0; i<bunch->getSize();i++)
		CalcPopulations(i, bunch);
}
		


void DM_noLaserField::finalizeEffects(Bunch* bunch){
		
        delete	[]	prob;
        

}


void DM_noLaserField::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{

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
		
			CalcPopulations(i, bunch);

	//cout<<scientific<<setprecision(20)<<bunch->x(0)<<"\t"<<bunch->y(0)<<"\t"<<bunch->z(0)<<"\t"<<PopAttr->attArr(i)[0]<<"\t" <<PopAttr->attArr(i)[12] + PopAttr->attArr(i)[13]+PopAttr->attArr(i)[7] + PopAttr->attArr(i)[8] + PopAttr->attArr(i)[1]<<"\n";
}

void DM_noLaserField::GetParticleFrameFields(int i,double t,double t_step,  Bunch* bunch,  BaseFieldSource* fieldSource)	{
	double Ez;
	
	fieldSource->getElectricMagneticField(x0,y0,z0,t,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);		
	LorentzTransformationEM::transform(bunch->getMass(),px0,py0,pz0,Ex_stat,Ey_stat,Ez_stat,Bx_stat,By_stat,Bz_stat);
	
	Ez=sqrt(Ex_stat*Ex_stat+Ey_stat*Ey_stat+Ez_stat*Ez_stat);
	
	Ex_stat=0;	
	Ey_stat=0;		
	Ez_stat=Ez;
}

void	DM_noLaserField::GetParticleFrameParameters(int i, double t,double t_step, Bunch* bunch)	{		
	
	double ta=2.418884326505e-17;			//atomic unit of time
	double Ea=5.14220642e011;				//Atomic unit of electric field
	double m=bunch->getMass();
	
	//This line calculates relyativistic factor-Gamma
	double gamma=sqrt(m*m+px0*px0+py0*py0+pz0*pz0)/m;
	double coeff=1./(gamma*ta);
	
	part_t_step=t_step*coeff;	//time step in frame of particle (in atomic units)
	t_part=t*coeff;	
	
	
	Ex_stat/=Ea;	
	Ey_stat/=Ea;		
	Ez_stat/=Ea;	
}



void DM_noLaserField::AmplSolver4step(double t_step, int i, Bunch* bunch)	{

	double dt,sum;
        
        
        if(bunch->flag(i) == 0)           {
        
	StarkEffect->SetE(Ez_stat);	
	
	for(int n=2; n<levels+1;n++) {
		for(int m=1; m<n;m++) {
			if (cond[n][m] && StarkEffect->field_thresh[n]>Ez_stat && StarkEffect->field_thresh[m]>Ez_stat)	{
				gamma_ij[n][m] = StarkEffect->getRelax(n,m);

			}
		}
	}
	
        
	for(int j=1; j<5; j++){
		for(int m=1; m<levels+1;m++)	{
			if (StarkEffect->field_thresh[m]>Ez_stat)	{	
				
				if (j==4)	dt=part_t_step;	else dt=part_t_step/2.;
	
				k_RungeKutt[m][j]=0;
				
				for(int k=m+1;k<levels+1;k++)	 if (cond[k][m] && StarkEffect->field_thresh[k]>Ez_stat) 	k_RungeKutt[m][j]+=gamma_ij[k][m]*(pop(i,k)+k_RungeKutt[k][j-1]*dt);			 
				sum=0;for(int k=1;k<m;k++)			 if (cond[m][k] && StarkEffect->field_thresh[k]>Ez_stat) 	sum+=gamma_ij[m][k];
				
				sum+=StarkEffect->Gamman[m];
				k_RungeKutt[m][j]-=(pop(i,m)+k_RungeKutt[m][j-1]*dt)*sum;
				
				
			}
		}
	}
	
	for(int m=1;m<levels+1;m++)	{
		
		if (StarkEffect->field_thresh[m]>Ez_stat)	{
			pop(i,m) += part_t_step*(k_RungeKutt[m][1]+2.*k_RungeKutt[m][2]+2.*k_RungeKutt[m][3]+k_RungeKutt[m][4])/6.;	
		}
		else {
			pop(i,m) = 0;
		}
		
	}


	
        }
        

			
		

	if((index == 1)&&(pop(i,0) >= prob[i]))
		bunch->recoverParticle(i);

        
        
}


