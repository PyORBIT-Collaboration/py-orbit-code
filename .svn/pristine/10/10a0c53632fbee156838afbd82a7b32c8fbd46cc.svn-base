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

#include "StarkStrongField.hh"
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




StarkStrongField::StarkStrongField(std::string addressEG,int n1,int n2,int mm)
{		

	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	std::ifstream file;
	std::ifstream file1;
	std::ifstream file2;
	double F,E0,E1;
	char G_ch[200];
	m = abs(mm);

	int fi;
	char nameEG[MAX_LENGTH_ADDRESS];
	char nameEG1[MAX_LENGTH_ADDRESS];
	char nameEG2[MAX_LENGTH_ADDRESS];





	


	

	
	if(rank_MPI == 0) {
		
		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%sStarkEG/000.txt",addressEG.c_str());	
		file.open(nameEG);	file>>F>>F>>F>>delta_F; file.clear();file.close();		

		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%sStarkEG/%i%i%i.txt",addressEG.c_str(),n1,n2,abs(m));		
		n_data=-1;	file.open(nameEG);	while(!file.eof())	{file>>field_thresh>>F>>G_ch; n_data++;} file.clear();file.close();	
					

	}		


	ORBIT_MPI_Bcast(&delta_F, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&n_data,1 ,MPI_INT,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&field_thresh,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);




	
	deltaE = new double[n_data];		
	gamma_autoionization = new double[n_data];	
	dipole_transition = new double[n_data];
	

	
				
		if(rank_MPI == 0) {
		snprintf(nameEG,MAX_LENGTH_ADDRESS,"%sStarkEG/000.txt",addressEG.c_str());
		snprintf(nameEG1,MAX_LENGTH_ADDRESS,"%sStarkEG/%i%i%i.txt",addressEG.c_str(),n1,n2,abs(m));	
		snprintf(nameEG2,MAX_LENGTH_ADDRESS,"%sTDM_StrongField/TDM_%i%i%i.txt",addressEG.c_str(),n1,n2,abs(m));
		file.open(nameEG);
		file1.open(nameEG1);
		file2.open(nameEG2);
		for (fi=0; fi<n_data; fi++)		{file>>F>>E0>>G_ch; 
										file1>>F>>E1>>G_ch; gamma_autoionization[fi] = atof(G_ch);
										file2>>F>>dipole_transition[fi];

										
										deltaE[fi] = E1 - E0;}	
										file.clear();file.close();
										file1.clear();file1.close();
										file2.clear();file2.close();
		}
						
					
			for(int fi=0;fi<n_data;fi++)	{
			ORBIT_MPI_Bcast(&deltaE[fi],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
			ORBIT_MPI_Bcast(&gamma_autoionization[fi],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
			ORBIT_MPI_Bcast(&dipole_transition[fi],1, MPI_DOUBLE,0,MPI_COMM_WORLD);
			}

			
	
}


StarkStrongField::~StarkStrongField()	{

	delete [] deltaE;
	delete [] gamma_autoionization;
	delete [] dipole_transition;

}





void	StarkStrongField::SetE(double E){	
	
	double c = E/delta_F;
	int iEz=(int)c;
	double cEz=c-iEz;
	

	if(E>field_thresh)	{deltaEn = 0; Gamman = 0; dipole_transition = 0;}	
	else	{
	deltaEn = deltaE[iEz]+cEz*(deltaE[iEz+1]-deltaE[iEz]);
	dipole_transitionn = dipole_transition[iEz]+cEz*(dipole_transition[iEz+1]-dipole_transition[iEz]);
	
	if(gamma_autoionization[iEz]<1e-100)
	Gamman = 0;
	else	
	Gamman = gamma_autoionization[iEz]*pow(gamma_autoionization[iEz+1]/gamma_autoionization[iEz],cEz);
	}


}

	

double	StarkStrongField::getStarkStrongFielddeltaE(
	double mass,  
	double Ex,double Ey,double Ez,
	double Bx,double By,double Bz,
	double px,double py,double pz){



LorentzTransformationEM::transform(mass,px,py,pz,Ex,Ey,Ez,Bx,By,Bz);

SetE(sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/5.14220642e011);	

return deltaEn;

	}
	
	
		
	

