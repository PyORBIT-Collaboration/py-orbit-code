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

#include "RecordEvolution.hh"
#include "RungeKuttaTracker.hh"




#define t_step_ev(i)   	Evol->attArr(i)[num_plot+1]
#define x0_ev(i)   		Evol->attArr(i)[num_plot+2]
#define y0_ev(i)   		Evol->attArr(i)[num_plot+3]
#define z0_ev(i)   		Evol->attArr(i)[num_plot+4]

using namespace LaserStripping;
using namespace OrbitUtils;



RecordEvolution::RecordEvolution(std::string effect,int ind_effect, int num)
{
	

	setRankSetup(1);
	setRankMemorize(1);
	setRankFinalize(1);
	setRankApply(1);

	num_plot=num;
	index_effect=ind_effect;
	effect_name=effect;


}

RecordEvolution::~RecordEvolution() 
{
}

void RecordEvolution::setupEffects(Bunch* bunch){
	setup_par = true;
	
	if(bunch->hasParticleAttributes("Evolution")==1)
		bunch->removeParticleAttributes("Evolution");
	
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = num_plot+5;
	bunch->addParticleAttributes("Evolution",part_attr_dict);
	
	Evol = bunch->getParticleAttributes("Evolution");

	RecEff = bunch->getParticleAttributes(effect_name);
			

	for (int i=0; i<bunch->getSize();i++)                        
	Evol->attArr(i)[0] = RecEff->attArr(i)[index_effect];
}


void RecordEvolution::applyEffectsForEach(Bunch* bunch, int i, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker)		{
	



	if (setup_par == true)	{
		
		x0 = y_in_vct[0];y0 = y_in_vct[1];z0 = y_in_vct[2];

		Num = tracker->getStepsNumber();
		
		if(Num<num_plot)	{cout<<"The number of evolution points must be equal or less then the number of steps \n"; abort();}
		if(Num%num_plot!=0)	{cout<<"The number of steps must be divisible by the number of evolution point \n"; abort();}
		
		for (int i=0; i<bunch->getSize();i++)	{
			x0_ev(i) = x0;
			y0_ev(i) = y0;
			z0_ev(i) = z0;
			t_step_ev(i) = Num*t_step/num_plot;
		}
		setup_par = false;
		t_in = t;
	}

	if(int((t-t_in+t_step)/t_step+0.5)%(Num/num_plot) == 0) {
		for (int i=0; i<bunch->getSize();i++)	{
			Evol->attArr(i)[int((t-t_in+t_step)/t_step+0.5)/(Num/num_plot)] = RecEff->attArr(i)[index_effect];
			//	Evol->attArr(i)[int((t-t_in+t_step)/t_step+0.5)/(Num/num_plot)] = bunch->x(i);

		}
	}
	
}
