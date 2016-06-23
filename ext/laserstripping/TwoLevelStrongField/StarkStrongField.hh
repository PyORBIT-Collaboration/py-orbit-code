#ifndef STARKSTRONGFIELD_HH_
#define STARKSTRONGFIELD_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>









		
		class  StarkStrongField: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		StarkStrongField(std::string addressEG,int n1,int n2,int m);
		
		/*destructor*/	
		~StarkStrongField();
			

//			/*this method gives frequendy transitions as a function of field */
		double getStarkStrongFielddeltaE(
				double mass, 
				double Ex,double Ey,double Ez,
				double Bx,double By,double Bz,
				double px,double py,double pz);
		
			/*this method setups electrostatic field */
		void SetE(double E);
		
		/*Allows to obtain maximum general quantum number n*/
		int getStates();
		
		double Gamman;
		double deltaEn;
		double dipole_transitionn;
		int m;
			
	
		
		private:  
		  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
			
			
		double* deltaE;		
		double* gamma_autoionization;	
		double* dipole_transition;

		  
		double delta_F;
		double field_thresh;
		int n_data;		 
		
		};




#endif /*StarkStrongField_HH_*/
