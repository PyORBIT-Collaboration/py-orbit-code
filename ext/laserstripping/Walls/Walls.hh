#ifndef WALLS_HH_
#define WALLS_HH_



#include "Python.h"
#include "ExternalEffects.hh"
#include "BaseLaserFieldSource.hh"
#include "ParticleAttributes.hh"




using namespace TrackerRK4;

namespace LaserStripping{
	
	class  Walls: public ExternalEffects
	{
		public:
		
			/** Constructor. */
			Walls();
			
			
			/** Destructor. */
			~Walls();
			

		
		/** It initializes effects. */
		void setupEffects(Bunch* bunch);
		
		/*it memorizes initial coordinates and impulses before rk step*/
		void memorizeInitParams(Bunch* bunch);
		
		/** It finalizes effects. */
		void finalizeEffects(Bunch* bunch);

		/** It applies the external effects to a particle with certain index. */
		void applyEffects(Bunch* bunch, 
														  double t, double t_step, 
														  OrbitUtils::BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker);	
		


		  private:
			  
			  ParticleAttributes* Coords;
			  double** xyz; 
		 

			  	

				bool crossSurface(int i,Bunch* bunch);
				
				bool insideSphere(int i,Bunch* bunch,double x0, double y0,double z0,double R);
				
				
				
				bool crossPlane(int i,Bunch* bunch,double x0, double y0,double z0,double vx, double vy,double vz);
				
				bool crossSphere(int i,Bunch* bunch,double x0, double y0,double z0,double R);
				
				bool crossCylinder(int i,Bunch* bunch, double x0, double y0,double z0,double vx, double vy,double vz,double R);
				
				bool crossCone(/*Parameters of cone here*/);

							

				
				
	};
};



#endif /*WALLS_HH_*/


