#ifndef FUNCT_HH_
#define FUNCT_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>
#include <mpc.h>






		
		class  Funct: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		Funct(int n11,int n22, int mm, int point11);
		
		/*destructor*/	
		~Funct();
		
		long int calcPrecisionForM(int& err_exp, int exp_minG, int max_err_exp, std::string field,std::string c_energy, std::string c_Z1);
		long int calcPrecisionForN(std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);


			


		std::string getM(std::string c_F,std::string c_energy, std::string c_Z1);
		
		long int getN(  std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);
		long int get_a( std::string&  str_a, std::string& str_der_a,std::string c_F,std::string c_energy, std::string c_Z2);		
		long int get_b( std::string&  str_b, std::string& str_der_b,std::string c_F,std::string c_energy, std::string c_Z2);
		
		std::string getB(std::string c_F,std::string c_energy, std::string c_Z2);
		



		


		private:  
			

			
			int n1;
			int n2;
			int m;
			int n;
			
			long int point1;
			long int point2;
			
			long int pointN1;
			long int pointN2;
			long int out_lenN;
			
			long int precisionN;
			long int prec_bitN;

			long int nsum1;
			long int nsumab;
			long int nsumN;
			
			long int prec_bit;
			long int precision;
			long int out_len;

			

			

			

			mpc_t *C;
			mpc_t *a;
			mpc_t *b;
			

		

			


		
		};




#endif /*FUNCT_HH_*/
