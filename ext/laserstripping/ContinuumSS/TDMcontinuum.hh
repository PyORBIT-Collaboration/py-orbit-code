#ifndef TDMCONTINUUM_HH_
#define TDMCONTINUUM_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>









		
		class  TDMcontinuum: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		TDMcontinuum(std::string addr);
		
		/*destructor*/	
		~TDMcontinuum();

//			/*this method gives frequendy transitions as a function of field */
		double setField_returndE(
				double mass, 
				double Ex,double Ey,double Ez,
				double Bx,double By,double Bz,
				double px,double py,double pz);
		
			/*this method setups electrostatic field */
		void SetE(double E);
		int get_levels();
		
		void get_dE_TDM(int i,double& dE, double& TDM);
		
		/*Allows to obtain maximum general quantum number n*/
		int getStates();
		
		double Gamman;
		double deltaEn;
		double dipole_transitionn;
		double deltaE;
		double Fmin;
		double Fmax;
		double dF;
		double Field;
		int m;
		
		int num_files;
			
	
		
		private:  
		  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
			
			
		int levels;		
		double* dEmin;
		double* gamma_autoionization;	
		double* dipole_transition;
		double* F;
		
		double** E_F_tdm;


		  
		double delta_F;
		double field_thresh;
		int n_data;		 
		
		};




#endif /*TDMcontinuum_HH_*/
