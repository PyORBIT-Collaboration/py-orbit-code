#ifndef CONTINUUM_SS_HH_
#define CONTINUUM_SS_HH_



#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "Stark.hh"
#include "ParticleAttributes.hh"
#include "TDMcontinuum.hh"


using namespace TrackerRK4;

namespace LaserStripping{
	
	class  ContinuumSS: public ExternalEffects
	{
		public:
		
			/** Constructor */
			ContinuumSS(OrbitUtils::BaseLaserFieldSource*	BaseLaserField,TDMcontinuum* cont);

			
			/** Destructor. */
			~ContinuumSS();
			

		
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

		  
			  
			  OrbitUtils::BaseLaserFieldSource*	LaserField;
			  
			  //this array is used on each step of solution of density matrix equation at definite field  
			  tcomplex** exp_mu_El;
			  tcomplex** k_RungeKutt;

			  		  

			  
			  int levels;
			  double delE;
			  double** int_E;
			  double* dE;
			  double* TDM;

			  
			  ParticleAttributes* AmplAttr;
			  ParticleAttributes* PopAttr;
			  
		  

			  			 
			  //time and frequensy of laser in frame of particle
			  double omega_part;
			  double part_t_step;

			  

			  					  
			  //Fields in the particle frame  
			  tcomplex Ex_las[3];
			  tcomplex Ey_las[3];
			  tcomplex Ez_las[3];
			  
			  tcomplex E_las[3];
			  
			  tcomplex Bx_las[3];
			  tcomplex By_las[3];
			  tcomplex Bz_las[3];
			  
			  tcomplex exp_phasa[3];
			  

			  
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
				void AmplSolver4step(int i,Bunch* bunch);
				
				/**Calculates populations**/
				void CalcPopulations(int i,Bunch* bunch);
							
				/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
				void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetParticleFrameFields(int i, double t,double t_step,  Bunch* bunch,  OrbitUtils::BaseFieldSource* fieldSource);
				
				

				
	};
};



#endif /*ContinuumSS_HH_*/


