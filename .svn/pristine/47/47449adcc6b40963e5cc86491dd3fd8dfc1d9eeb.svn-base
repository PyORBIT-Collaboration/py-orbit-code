#ifndef STARK_HH_
#define STARK_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>









		
		class  Stark: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		Stark(std::string addressEG,int states);
		
		/*destructor*/	
		~Stark();
			
			/*this method gives koefficient of autoionization as a function of field  */
		double getRelax(int k,int ks);
		double getRelaxTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z);
		void   getTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z);

		
//			/*this method gives frequendy transitions as a function of field */
		double getStarkEnergy(
				double mass, int n1, int n2, int m, 
				double Ex,double Ey,double Ez,
				double Bx,double By,double Bz,
				double px,double py,double pz);
		
			/*this method setups electrostatic field */
		void SetE(double E);
		
		/*Allows to obtain maximum general quantum number n*/
		int getStates();
		
		
		double* field_thresh;	
		double* En;
		double* Gamman;
		int* n_data;
		
		
		private:  
		  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
		  
		 double*** dipole_transition_x;
		 double*** dipole_transition_y;
		 double*** dipole_transition_z;
		 double*** gamma_spontaneous_relax;
		 double** energy;
		 double** gamma_autoionization;
		 
		 double* E_pow_n;

		  
		 double delta_F;
		 double const_relax;
		 double Ez_stat;
		 double cEz;
		 
		 int iEz;
		 int st;
		 int levels;
		 int order_trans;
		
		};




#endif /*STARK_HH_*/
