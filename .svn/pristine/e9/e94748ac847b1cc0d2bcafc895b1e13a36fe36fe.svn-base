#ifndef FUNCT_HH_
#define FUNCT_HH_


//#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>
//#include <mpc.h>
#include <gmpfrxx.h>



		
		class  FuncSS: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		FuncSS(int n11,int n22, int mm, int point11);
		
		/*destructor*/	
		~FuncSS();
		
		long int calcPrecisionForM(int max_err_exp, std::string c_field, std::string c_energy, std::string c_Z1);
		long int calcPrecisionForN(std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);


			


		std::string getM(std::string c_F,std::string c_energy, std::string c_Z1);
		mpfr_class integrateM2(mpfr_class F,mpfr_class E, mpfr_class Z1);
		
		long int getN(  std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);
		long int get_ab(std::string&  a, std::string& der_a,std::string&  b, std::string& der_b,std::string c_F,std::string c_energy, std::string c_Z2);

		
		std::string getB(std::string c_F,std::string c_energy, std::string c_Z2);
		std::string C(std::string c_F,std::string c_energy, std::string c_Z2);
		std::string M(std::string c_mu,std::string c_F, std::string c_energy, std::string c_Z1);
		std::string N(std::string nu,std::string c_F, std::string c_energy, std::string c_Z2);
		



		


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



			

			

//			mpc_t *C;
//			mpc_t *a;
//			mpc_t *b;
			

		

			


		
		};




#endif /*FUNCSS_HH_*/
