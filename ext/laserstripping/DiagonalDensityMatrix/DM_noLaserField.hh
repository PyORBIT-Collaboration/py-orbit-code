#ifndef DM_NOLASERFIELD_HH_
#define DM_NOLASERFIELD_HH_

#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "Stark.hh"
#include "ParticleAttributes.hh"

using namespace TrackerRK4;

namespace LaserStripping{
	
	class  DM_noLaserField: public ExternalEffects
	{
		public:
		
			/** Constructor */
			DM_noLaserField(Stark* StarkEffect, int method);

			
			/** Destructor. */
			~DM_noLaserField();
			
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
                
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);
		

		/** It applies the external effects to a particle with certain index. */
		void applyEffectsForEach(Bunch* bunch, int index, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker);	
		
		private:
			
			Stark* StarkEffect;
			
			//this array is used on each step of solution of density matrix equation at definite field  
			
			double** k_RungeKutt;
			double** gamma_ij;
			bool** cond;
                        double* prob;
			
			int levels;
                        int index;
			
			ParticleAttributes* PopAttr;
			
			//time and frequensy of laser in frame of particle
			
			double part_t_step;
			double t_part;
			
			//Fields in the particle frame  
		  
			double Ex_stat;
			double Ey_stat;
			double Ez_stat;
			
			double Bx_stat;
			double By_stat;
			double Bz_stat;
			
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
			
			
			/**Solver for Amplitudes**/
			void AmplSolver4step(double t_step, int i,Bunch* bunch);
			
			/**Calculates populations**/
			void CalcPopulations(int i,Bunch* bunch);
			
			/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
			void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
			
			/*this method gives laser and static fields transformed by rotation relatively z axes  */
			void GetParticleFrameFields(int i, double t,double t_step,   Bunch* bunch,  OrbitUtils::BaseFieldSource* fieldSource);

	};
};



#endif /*DM_NOLASERFIELD_HH_*/


