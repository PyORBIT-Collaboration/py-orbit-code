#ifndef WAVEFUNCTION_HH_
#define WAVEFUNCTION_HH_


#include "tcomplex.hh"
#include "CppPyWrapper.hh"
#include <string>
#include <mpc.h>
#include "MathPolynomial.hh"







		
		class  WaveFunction: public OrbitUtils::CppPyWrapper
		{

			
		public:
		/*constructor*/	
		WaveFunction(int n11,int n22, int mm,  int point11,  int pointN,  int precEZ, std::string c_energy, std::string c_Z1,std::string c_F);
		
		/*destructor*/	
		~WaveFunction();
		

		long int calcPrecisionForN(std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);

		long int getnsumM(std::string c_F, std::string c_energy, std::string c_Z1);
		long int getnsumN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z1);
		
		void getArrayM(std::string c_F,std::string c_energy, std::string c_Z1);
		void getArrayN(std::string c_F,std::string c_energy, std::string c_Z1);
		
		void getArray_a(std::string str_A,std::string c_F,std::string c_energy, std::string c_Z1);
		void getArray_b(std::string str_B,std::string c_F,std::string c_energy, std::string c_Z1);
		
		void fillArrayM();
		void fillArrayN();
		void fillArray_ab(std::string c_F,std::string c_energy);

		void getM(mpfr_t mu,mpc_t M);
		void getN(mpfr_t mu,mpc_t N);
		void get_ab(mpfr_t mu,std::string c_F,std::string c_energy, mpc_t ab);
		
		double ndivN_b(std::string c_energy,std::string c_F);
		double ndivN_h(std::string c_energy,std::string c_F);

		
		std::string getFastM(std::string  mu);
		std::string getFastNbelow(std::string  mu);
		std::string getFast_ab(std::string  mu);
		std::string getFastN(std::string  mu);
		
		long int get_a( std::string&  str_a, std::string& str_der_a,std::string c_F,std::string c_energy, std::string c_Z1);		
		long int get_b( std::string&  str_b, std::string& str_der_b,std::string c_F,std::string c_energy, std::string c_Z1);
		
		void getAB(std::string str_a,std::string str_der_a,std::string str_b,std::string str_der_b,std::string str_N,std::string str_der_N,std::string& str_A,std::string& str_B);
			
		long int getN(  std::string&  str_N, std::string& str_der_N,std::string c_F,std::string c_energy, std::string c_Z2);
		
		int getMode();
		
		void getArrayb();
		double getMN(double mu,double nu);
		double f(long int k, long int _m, double x);

		
		
		
		  


		


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
			
			long int err_exp;
			
			long int prec_bit;
			long int precision;
			long int out_len;
			
			long int ndivM;
			long int ndivN;
			long int ndivab;
			
			int mode;
			
			double Fd;
			double A;
			

			

			

			mpc_t *c;
			mpc_t *d;
			
			mpc_t *Mi;
			mpc_t *Ma;
			mpc_t *Mb;
			mpc_t *Mc;
			mpc_t *Md;
			
			mpc_t *Ni;
			mpc_t *Na;
			mpc_t *Nb;
			mpc_t *Nc;
			mpc_t *Nd;
			
			mpc_t *abi;
			mpc_t *ab_a;
			mpc_t *ab_b;
			mpc_t *ab_c;
			mpc_t *ab_d;
			
			mpc_t *a;
			mpc_t *b;
			mpfr_t h;
			
			long int*** bp;
			

			

		

			


		
		};




#endif /*WAVEWaveFunction_HH_*/
