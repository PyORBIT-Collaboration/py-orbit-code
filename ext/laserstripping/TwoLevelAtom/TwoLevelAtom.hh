#ifndef TWOLEVELATOM_HH_
#define TWOLEVELATOM_HH_

#include "Python.h"

#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "ParticleAttributes.hh"

using namespace TrackerRK4;

namespace LaserStripping{
	
	class  TwoLevelAtom: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			TwoLevelAtom(OrbitUtils::BaseLaserFieldSource*	BaseLaserField, double delta_E, double dipole_tr);
			
			
			/** Destructor. */
			~TwoLevelAtom();
			

		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		

		/** It applies the external effects to a particle with certain index. */
		void applyEffectsForEach(Bunch* bunch, int index, 
			                            double* y_in_vct, double* y_out_vct, 
																  double t, double t_step, 
																  OrbitUtils::BaseFieldSource* fieldSource,
																	RungeKuttaTracker* tracker);	
		
		  private:
	  
			  OrbitUtils::BaseLaserFieldSource*	LaserField;
			  //this array is used on each step of solution of density matrix equation at definite field  

			  tcomplex mu_Elas[5];
			  tcomplex exp_mu_El[5];
			  tcomplex k_RungeKutt_1[5];
			  tcomplex k_RungeKutt_2[5];
			  double d_Energy;
			  double dip_transition;
			  
			  
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
			  
			  ParticleAttributes* AmplAttr;
			  ParticleAttributes* PopAttr;
		  

			  			 
			  //time and frequensy of laser in frame of particle
			  double part_t_step;
			  double t_part;

			  					  
			  //Fields in the particle frame  
			  
			  tcomplex Ex_las[3];
			  tcomplex Ey_las[3];
			  tcomplex Ez_las[3];
			  
			  tcomplex Bx_las[3];
			  tcomplex By_las[3];
			  tcomplex Bz_las[3];
			  
			  tcomplex exp_phasa[3];
			  
			  double Ex_stat;
			  double Ey_stat;
			  double Ez_stat;
			  
			  double Bx_stat;
			  double By_stat;
			  double Bz_stat;
			  	
				/**Solver for Amplitudes**/
				void AmplSolver4step(int i,Bunch* bunch);
							
				/**Calculates populations**/
				void CalcPopulations(int i,Bunch* bunch);
				
				/*this all parameters in frame of particle in atomic units necesarry for applying external effect  */
				void GetParticleFrameParameters(int i, double t, double t_step, Bunch* bunch);
				
				/*this method gives laser and static fields transformed by rotation relatively z axes  */
				void GetParticleFrameFields(int i, double t,double t_step,  Bunch* bunch);				
	};
};



#endif /*TWOLEVELATOM_HH_*/


