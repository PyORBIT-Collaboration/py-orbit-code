#ifndef HYDROGENSTARKPARAM_HH_
#define HYDROGENSTARKPARAM_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"










		
		class  HydrogenStarkParam: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		HydrogenStarkParam(char* addressEG,int states);
		
		/*destructor*/	
		~HydrogenStarkParam();
			
			/*this method gives koefficient of autoionization as a function of field  */
		void GetEnergyAutoionization(int k,double& E, double& Gamma);
			
			/*this method gives frequendy transitions as a function of field */
		void GetDipoleTransition(int k,int ks,tcomplex& mu_x,tcomplex& mu_y,tcomplex& mu_z);	
			
			/*this method gives frequendy transitions as a function of field */
		double GetRelax(int k,int ks);
		
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
		
		
		private:  
		  //parameters of dipole transition, energy, lifetime, spontaneous relaxation, 
		  
		 double*** dipole_transition_x;
		 double*** dipole_transition_y;
		 double*** dipole_transition_z;
		 double*** gamma_spontaneous_relax;
		 double** energy;
		 double** gamma_autoionization;

		  
		 double delta_F;
		 int n_data;
		 double Ez_stat;
		 int st;
		 int levels;
		
		};




#endif /*HYDROGENSTARKPARAM_HH_*/
