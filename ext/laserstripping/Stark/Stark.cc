//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    HydrogenStarkParam.cc
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
#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include "Stark.hh"
#include "LorentzTransformationEM.hh"
#include <cmath>
#include <fstream>
#include <iostream>

#define MAX_LENGTH_ADDRESS 1024

using namespace OrbitUtils;
using namespace std;

inline int convert3to1level(int n,int n1, int m){
	return 1+m*n+n1+(n*n*n-n)/3-m*abs(m-1)/2;
}




Stark::Stark(std::string addressEG,int states)
{		

	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	std::ifstream file;
	double F;
	char G_ch[200];
	double alpha=7.297352570e-3;
	int k,ks,fi,nn1,nn2,mm,nn1s,nn2s,mms,ff;
	char nameEG[MAX_LENGTH_ADDRESS];
	std::string dump,dump2;

	
	st=states;
	const_relax = 4*alpha*alpha*alpha/3;
	levels=st*(1+st)*(1+2*st)/6;


	
  	field_thresh = new double[levels+1];
  	n_data = new int[levels+1];
  	En = new double[levels+1];
  	Gamman = new double[levels+1];	
	

	
	if(rank_MPI == 0) {
		n_data[0] = 1;
		
		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%s000.txt",addressEG.c_str());	
		file.open(nameEG);	file>>F>>F>>F>>delta_F; file.clear();file.close();

		
		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%stransitions.txt",addressEG.c_str());
		file.open(nameEG);	file>>ff>>ff>>ff>>dump>>ff>>ff>>ff; fi=0;	while(dump2!=dump)	{file>>dump2;fi++;} file.clear();file.close(); order_trans = (fi-4)/3;
		
		
		for(int n=1;n<states+1;n++){
			for(int m=-(n-1);m<(n-1)+1;m++){
				for(int n1=0;n1<n-abs(m)-1+1;n1++)	{

					k=convert3to1level(n,n1,m);		
					snprintf(nameEG,MAX_LENGTH_ADDRESS,"%s%i%i%i.txt",addressEG.c_str(),n1,n-n1-abs(m)-1,abs(m));		
					n_data[k]=-1;	file.open(nameEG);	while(!file.eof())	{file>>field_thresh[k]>>F>>G_ch; n_data[k]++;} file.clear();file.close();	
					
				}
			}
		}
		
	}		


	
	ORBIT_MPI_Bcast(&order_trans, 1, MPI_INT,0,MPI_COMM_WORLD);	
	ORBIT_MPI_Bcast(&delta_F, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
	for (int i=0;i<levels+1;i++) ORBIT_MPI_Bcast(&n_data[i],1 ,MPI_INT,0,MPI_COMM_WORLD);
	for (int i=0;i<levels+1;i++) ORBIT_MPI_Bcast(&field_thresh[i],1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);



	//allocating memory for dynamic massive of data that will be read fronm files	
	E_pow_n = new double[order_trans];
	

	
	energy=new double*[levels+1];	
	for (int i=0;i<levels+1;i++)	energy[i]=new double[n_data[i]];
	
	gamma_autoionization=new double*[levels+1];	
	for (int i=0;i<levels+1;i++)	gamma_autoionization[i]=new double[n_data[i]];
	

	
	dipole_transition_x=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_x[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_x[i][j]=new double[order_trans];
	
	dipole_transition_y=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_y[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_y[i][j]=new double[order_trans];
	
	dipole_transition_z=new double**[levels+1];	
	for (int i=0;i<levels+1;i++)	dipole_transition_z[i]=new double*[levels+1]; 
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	dipole_transition_z[i][j]=new double[order_trans];	
	

  

  
	

				
		for(int n=1;n<states+1;n++){
			for(int m=-(n-1);m<(n-1)+1;m++){
				for(int n1=0;n1<n-abs(m)-1+1;n1++){

					k=convert3to1level(n,n1,m);		
				
					if(rank_MPI == 0) {
					snprintf(nameEG,MAX_LENGTH_ADDRESS,"%s%i%i%i.txt",addressEG.c_str(),n1,n-n1-abs(m)-1,abs(m));						
					file.open(nameEG);for (fi=0; fi<n_data[k]; fi++)	{file>>F>>energy[k][fi]>>G_ch; gamma_autoionization[k][fi] = atof(G_ch);}	file.clear();file.close();
					}
						
					
					for(int fi=0;fi<n_data[k];fi++)	{
					ORBIT_MPI_Bcast(&energy[k][fi],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
					ORBIT_MPI_Bcast(&gamma_autoionization[k][fi],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
					}
					
				}
			}
		}
		
	 



	//this double loop reads and broadcast dipole transitions 
	if(rank_MPI == 0) {
	

		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%stransitions.txt",addressEG.c_str());
		file.open(nameEG);
	
	}

		for(int n=1;n<states+1;n++){
			for(int m=-(n-1);m<(n-1)+1;m++){
				for(int n1=0;n1<n-abs(m)-1+1;n1++){
					
					for(int ns=1;ns<states+1;ns++){
						for(int ms=-(ns-1);ms<(ns-1)+1;ms++){
							for(int n1s=0;n1s<ns-abs(ms)-1+1;n1s++){
								
				
								k=convert3to1level(n,n1,m);
								ks=convert3to1level(ns,n1s,ms);
								
								
								
								if(rank_MPI == 0)	{
								while(true)	{
								file>>nn1>>nn2>>mm>>dump>>nn1s>>nn2s>>mms;
						
								if(nn1==n1&&nn2==n-n1-fabs(m)-1&&mm==m&&nn1s==n1s&&nn2s==ns-n1s-fabs(ms)-1&&mms==ms)	{

									
										for (fi=0; fi<order_trans; fi++) file>>dipole_transition_x[k][ks][fi];
										for (fi=0; fi<order_trans; fi++) file>>dipole_transition_y[k][ks][fi];
										for (fi=0; fi<order_trans; fi++) file>>dipole_transition_z[k][ks][fi];
																				
									break;
									
								}
								else for (fi=0; fi<order_trans; fi++) file>>F>>F>>F;
								}
								}
									
								
								for(int fi=0;fi<order_trans;fi++)	{
								ORBIT_MPI_Bcast(&dipole_transition_x[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
								ORBIT_MPI_Bcast(&dipole_transition_y[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
								ORBIT_MPI_Bcast(&dipole_transition_z[k][ks][fi], 1, MPI_DOUBLE,0,MPI_COMM_WORLD);		
								
//								if(rank_MPI == 0&&k==1&&ks==2)cout<<dipole_transition_x[1][2][fi]<<"\n";
//								if(rank_MPI == 1&&k==1&&ks==2)cout<<dipole_transition_x[1][2][fi]<<"\n";
								}
								


								
								}

	
							}
						}
					}
				}
			}
		if(rank_MPI == 0)	{file.clear();file.close();}
	


	
		
}


Stark::~Stark()	{

	
	for (int i=0;i<levels+1;i++)	delete	[]	energy[i];	delete	[]	energy;
	for (int i=0;i<levels+1;i++)	delete	[]	gamma_autoionization[i];	delete	[]	gamma_autoionization;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_x[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_x[i];	delete	[]	dipole_transition_x;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_y[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_y[i];	delete	[]	dipole_transition_y;
	for (int i=0;i<levels+1;i++) for (int j=0;j<levels+1;j++)	delete [] dipole_transition_z[i][j]; for (int i=0;i<levels+1;i++)	delete [] dipole_transition_z[i];	delete	[]	dipole_transition_z;

	delete [] field_thresh;
	delete [] E_pow_n;
	delete [] n_data;
	delete [] En;
	delete [] Gamman;


}





void Stark::getTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z){
	
	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	double mu_xd=0, mu_yd=0, mu_zd=0;
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] = E_pow_n[fi-1]*Ez_stat;	
	
	
	
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
			
	mu_x = tcomplex(mu_xd,0);
	mu_y = tcomplex(0,mu_yd);
	mu_z = tcomplex(mu_zd,0);
	
	
	
}



double Stark::getRelax(int k,int ks){
	
	double mu_xd = 0, mu_yd = 0, mu_zd = 0;
	
	double absdE = fabs(En[k]-En[ks]);
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] = E_pow_n[fi-1]*Ez_stat;
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
		
	return const_relax*absdE*absdE*absdE*(mu_xd*mu_xd+mu_yd*mu_yd+mu_zd*mu_zd);
}



double Stark::getRelaxTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z){
	
	double mu_xd = 0, mu_yd = 0, mu_zd = 0;
	double absdE = fabs(En[k]-En[ks]);
	
	E_pow_n[0] = 1;
	for (int fi=1; fi<order_trans; fi++)  E_pow_n[fi] = E_pow_n[fi-1]*Ez_stat;
	
	
	for (int fi=0; fi<order_trans; fi++) mu_xd += dipole_transition_x[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_yd += dipole_transition_y[k][ks][fi]*E_pow_n[fi];
	for (int fi=0; fi<order_trans; fi++) mu_zd += dipole_transition_z[k][ks][fi]*E_pow_n[fi];
	
	
	mu_x = tcomplex(mu_xd,0);
	mu_y = tcomplex(0,mu_yd);
	mu_z = tcomplex(mu_zd,0);
		
	return const_relax*absdE*absdE*absdE*(mu_xd*mu_xd+mu_yd*mu_yd+mu_zd*mu_zd);
}




void	Stark::SetE(double E){	
	
	Ez_stat=E;	
	
	double c = Ez_stat/delta_F;
	iEz=(int)c;
	cEz=c-iEz;
	
	for(int k=1; k<levels+1;k++)	{
	if(Ez_stat>field_thresh[k])	{En[k] = 0; Gamman[k] = 0;}	
	else	{
	En[k] = energy[k][iEz]+cEz*(energy[k][iEz+1]-energy[k][iEz]);
	
	if(gamma_autoionization[k][iEz]<1e-21)
	Gamman[k] = 0;
	else	
	Gamman[k] = gamma_autoionization[k][iEz]*pow(gamma_autoionization[k][iEz+1]/gamma_autoionization[k][iEz],cEz);
	}
	}


}

int Stark::getStates() {return st;}



double	Stark::getStarkEnergy(
	double mass, int n1, int n2, int m, 
	double Ex,double Ey,double Ez,
	double Bx,double By,double Bz,
	double px,double py,double pz){



LorentzTransformationEM::transform(mass,px,py,pz,Ex,Ey,Ez,Bx,By,Bz);

SetE(sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/5.14220642e011);	

return En[convert3to1level(n1+n2+abs(m)+1,n1,m)];

	}
	
	
	
	

