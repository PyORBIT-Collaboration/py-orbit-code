#ifndef PRINTEXTEFFECTS_HH_
#define PRINTEXTEFFECTS_HH_

#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include <string>

#define MAX_LENGTH_ADDRESS 1024

using namespace TrackerRK4;

namespace LaserStripping{
	
	class  PrintExtEffects: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			PrintExtEffects(std::string eff_name,int i,std::string addr_print);
			
			/** Destructor. */
			~PrintExtEffects();
			
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);


		/** ??????????????????. */
		void applyEffectsForEach(Bunch* bunch, int i, 
                double* y_in_vct, double* y_out_vct, 
										  double t, double t_step, 
										  OrbitUtils::BaseFieldSource* fieldSource,
											RungeKuttaTracker* tracker);	
		
		  private:
				
			  int Num;
			  int num_print;
			  std::string  addr_print;
			  std::string  eff_name;
			  char addr_name[MAX_LENGTH_ADDRESS];
			  
			  int rank_MPI,size_MPI;
			  bool setup_par;
			  double t_in;	
	};
};



#endif /*PRINTEXTEFFECTS_HH_*/


