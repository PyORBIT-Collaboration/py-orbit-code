#ifndef RECORDEVOLUTION_HH_
#define RECORDEVOLUTION_HH_

#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "ParticleAttributes.hh"
#include <string>

using namespace TrackerRK4;

namespace LaserStripping{
	
	class  RecordEvolution: public ExternalEffects
	{
	public:
		
		/** Constructor. */
		RecordEvolution(std::string effect,int ind_effect, int num);
		
		/** Destructor. */
		~RecordEvolution();
		
		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
	

		
		/** ?????????????????. */
		void applyEffectsForEach(Bunch* bunch, int i, 
                double* y_in_vct, double* y_out_vct, 
										  double t, double t_step, 
										  OrbitUtils::BaseFieldSource* fieldSource,
											RungeKuttaTracker* tracker);	
		
		private:
			
			int num_plot;
			int index_effect;
			std::string effect_name;
			double t_per;
			double t_step;
			double t_in;
			int Num;
			bool setup_par;
			
			  double x0;
			  double y0;
			  double z0;
			  double px0;
			  double py0;
			  double pz0;
			  
			  double x;
			  double y;
			  double z;
			  double px;
			  double py;
			  double pz;
			
			
			ParticleAttributes* RecEff;
			ParticleAttributes* Coords;
			ParticleAttributes* Evol;
	};
};




#endif /*RECORDEVOLUTION_HH_*/
