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

#include "TDMcontinuum.hh"
#include "LorentzTransformationEM.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <errno.h>










#define MAX_LENGTH_ADDRESS 1024
#define MAX_NUMBER_OF_FILES 1024

using namespace OrbitUtils;
using namespace std;





TDMcontinuum::TDMcontinuum(std::string addr_trans)
{		

	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);	
	ifstream file;
	ifstream files[MAX_NUMBER_OF_FILES];
	double FF[MAX_NUMBER_OF_FILES];
	char addr_file[MAX_LENGTH_ADDRESS];
	DIR *dp;
	struct dirent *dirp;
	double E,temp,temp1,temp2;
	string temp_str;
	int k;


	
	

	if(rank_MPI == 0) {

		
	if((dp = opendir(addr_trans.c_str())) == NULL) 		cout << "Error(" << errno << ") opening " << addr_trans << endl;

	
	num_files = 0;
	levels=0;
	while ((dirp = readdir(dp)) != NULL)	
	if(string(dirp->d_name)[0]==string("T")[0])	{
		num_files++;
		
		snprintf(addr_file,MAX_LENGTH_ADDRESS,"%s%s",addr_trans.c_str(),string(dirp->d_name).c_str());
		files[num_files].open(addr_file);
		if(num_files==1)	{	
		while(!files[num_files].eof())	 
		{if (levels==0) files[num_files]>>temp1>>temp;	 else files[num_files]>>temp2>>temp;	levels++;}
		files[num_files].clear();files[num_files].close(); files[num_files].open(addr_file);	}
		
		deltaE = temp2 - temp1;
		
		temp_str = string(dirp->d_name);	temp_str.erase(0,6);	 temp_str.erase(temp_str.size()-4,4);  	FF[num_files] = atof(temp_str.c_str());
		
		if(num_files == 1)	{Fmin = FF[num_files]; Fmax = FF[num_files];}
		else	
		{
			if(FF[num_files] < Fmin)	Fmin = FF[num_files];
			if(FF[num_files] > Fmax)	Fmax = FF[num_files];
		}
		
	}

	
	closedir(dp);
	}
	
	
	
	ORBIT_MPI_Bcast(&levels,1 ,MPI_INT,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&num_files,1 ,MPI_INT,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&Fmin,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&Fmax,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	ORBIT_MPI_Bcast(&deltaE,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	
	
	dF = (Fmax - Fmin)/(num_files - 1);
	
	F = new double[num_files + 1];
	E_F_tdm = new double*[levels+1];	for (int i=0;i<levels+1;i++)	E_F_tdm[i] = new double[num_files + 1];
	dEmin = new double[num_files + 1];
	


	if(rank_MPI == 0) {		
	for (int i = 1; i<num_files + 1; i++)	{
		temp = 2*Fmax;
	for (int j = 1; j<num_files + 1; j++)	{
		
		if(FF[j]<temp)	{
			k = j;
			temp = FF[j];
		}
	
	}

	for (int n=1;n<levels+1;n++)	 {	files[k]>>temp>>E_F_tdm[n][i]; if (n==1)	dEmin[i] = 0.5 + temp; }	
	
	F[i] = FF[k];
	FF[k] = 2*Fmax;

	files[k].clear();	
	files[k].close();	
	
	}

	

	}
	


	for (int i=1;i<num_files+1;i++)	{
		ORBIT_MPI_Bcast(&dEmin[i],1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
		ORBIT_MPI_Bcast(&F[i],1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		for (int j=1;j<levels+1;j++)
		ORBIT_MPI_Bcast(&E_F_tdm[j][i],1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}


}	
			



TDMcontinuum::~TDMcontinuum()	{

	delete []	F;
	delete []	dEmin;
	for (int i=0;i<levels+1;i++) delete	[]	E_F_tdm[i];		delete [] E_F_tdm;
	
}




void	TDMcontinuum::SetE(double E){
	
	if(E>Fmax||E<Fmin)	{
		cout<<"error:	Electric field in the particle rest frame"<<"  E=  "<<E<<" lies outside of data on TDM of continuum spectra"<<" ("<<Fmin<<","<<Fmax<<") ";
		abort();
	}
	else
	Field = E;
	
	
}


int TDMcontinuum::get_levels()	{

	
	return levels;
}
	

double	TDMcontinuum::setField_returndE(double mass,  
	double Ex,double Ey,double Ez,
	double Bx,double By,double Bz,
	double px,double py,double pz){
	
	
	LorentzTransformationEM::transform(mass,px,py,pz,Ex,Ey,Ez,Bx,By,Bz);

	SetE(sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/5.14220642e011);
	
//	cout<<"E=	"<<sqrt(Ex*Ex+Ey*Ey+Ez*Ez)/5.14220642e011<<endl;
	
	double dE_, TDM_, maxTDM, dEmax;
	maxTDM = 0;
	
	for (int i = 1; i<levels + 1; i++)	{
		get_dE_TDM(i, dE_,  TDM_);
		if(TDM_>=maxTDM){
			maxTDM = TDM_;
			dEmax = dE_;
		}
		}
	
	return dEmax;

}

	
void TDMcontinuum::get_dE_TDM(int i,double& dE, double& TDM)	{
	
	double FF = Field - Fmin;
	int k = (int)(FF/dF);
	if(k<=0)	k=1;
	if(k>=num_files - 2)	k=num_files - 3;
	
	double a = ((2+k)*(k*(1+k)*E_F_tdm[i][k]-3*(k-1)*((1+k)*E_F_tdm[i][k+1]-k*E_F_tdm[i][k+2]))-k*(k*k-1)*E_F_tdm[i][k+3])/6;
	double b = (-(2+3*k*(2+k))*E_F_tdm[i][k]+3*(-1+k*(4+3*k))*E_F_tdm[i][k+1]-3*(-2+k*(2+3*k))*E_F_tdm[i][k+2]+(-1+3*k*k)*E_F_tdm[i][k+3])/(6*dF);
	double c = ((1+k)*E_F_tdm[i][k]-(2+3*k)*E_F_tdm[i][k+1]+(1+3*k)*E_F_tdm[i][k+2]-k*E_F_tdm[i][k+3])/(2*dF*dF);
	double d = (-E_F_tdm[i][k]+3*E_F_tdm[i][k+1]-3*E_F_tdm[i][k+2]+E_F_tdm[i][k+3])/(6*dF*dF*dF);
	
	
	double a1 = ((2+k)*(k*(1+k)*dEmin[k]-3*(k-1)*((1+k)*dEmin[k+1]-k*dEmin[k+2]))-k*(k*k-1)*dEmin[k+3])/6;
	double b1 = (-(2+3*k*(2+k))*dEmin[k]+3*(-1+k*(4+3*k))*dEmin[k+1]-3*(-2+k*(2+3*k))*dEmin[k+2]+(-1+3*k*k)*dEmin[k+3])/(6*dF);
	double c1 = ((1+k)*dEmin[k]-(2+3*k)*dEmin[k+1]+(1+3*k)*dEmin[k+2]-k*dEmin[k+3])/(2*dF*dF);
	double d1 = (-dEmin[k]+3*dEmin[k+1]-3*dEmin[k+2]+dEmin[k+3])/(6*dF*dF*dF);
	
	TDM = a + b*FF + c*FF*FF + d*FF*FF*FF;
	dE = a1 + b1*FF + c1*FF*FF + d1*FF*FF*FF	+	(i-1)*deltaE/(levels - 1);
	
	return;
}
	
	
		
	

