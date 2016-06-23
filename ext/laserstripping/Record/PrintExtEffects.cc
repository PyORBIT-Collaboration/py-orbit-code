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
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "PrintExtEffects.hh"
#include "RungeKuttaTracker.hh"

using namespace LaserStripping;
using namespace OrbitUtils;

PrintExtEffects::PrintExtEffects(std::string name,int i,std::string addr)
{
	setRankSetup(1);
	setRankMemorize(1);
	setRankFinalize(1);
	setRankApply(1);
	
	
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	num_print=i;
	addr_print=addr;
	eff_name = name;
}

PrintExtEffects::~PrintExtEffects() 
{
}

void PrintExtEffects::setupEffects(Bunch* bunch){		
	setup_par = true;
}


	
void PrintExtEffects::applyEffectsForEach(Bunch* bunch, int i, 
				                            double* y_in_vct, double* y_out_vct, 
																	  double t, double t_step, 
																	  OrbitUtils::BaseFieldSource* fieldSource,
																		RungeKuttaTracker* tracker)		{

	if (setup_par == true)	{
		Num = tracker->getStepsNumber();
		
		if(Num<num_print)		{cout<<"The number of record points must be equal or less then the number of steps \n"; abort();}
		if(Num%num_print!=0)	{cout<<"The number of steps must be divisible by the number of record points \n"; abort();}
		
		setup_par = false;
		t_in = t;
	}
	
	if(int((t-t_in+t_step)/t_step+0.5)%(Num/num_print) == 0)		{

			snprintf(addr_name,MAX_LENGTH_ADDRESS,"%s%i.dat",addr_print.c_str(),i*size_MPI+rank_MPI);
			ofstream file(addr_name,ios::app);
			file<<t;
			for (int j=0; j<bunch->getParticleAttributes(eff_name)->getAttSize();j++)	
			{
				file<<"\t"<<bunch->getParticleAttributes(eff_name)->attArr(i)[j];
			}
			file<<"\n";
			file.close();

	}
}
															
															
															
															























