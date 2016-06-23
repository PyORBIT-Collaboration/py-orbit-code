//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    HydrogenStarkParam.cc
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implements 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include "Functions.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>




using namespace OrbitUtils;
using namespace std;




Functions::Functions(int n11,int n22, int mm, int point11)	{		
	

	n1 = n11;
	n2 = n22;
	m = abs(mm);
	n = n1 + n2 + m + 1; 
	
	point1 = point11;
	point2 = point1*point1;
	
	out_lenN = 100;
	
	pointN1 = 55;
	pointN2 = pointN1*pointN1;
	
	
	precisionN = 1000;
	nsumab = 100000;
	prec_bitN = (int)(precisionN*3.3219280948873626);

	
	nsum1 = 10;
	precision = 10;
	prec_bit = (int)(precision*3.3219280948873626);


}





Functions::~Functions()	{
	

}




long int Functions::calcPrecisionForM(int& err_exp, int exp_minG, int max_err_exp, std::string c_field,std::string c_energy, std::string c_Z1)	{


	std::string line = "empty line", Mstr;
	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;
	int exponent;
	
	mpfr_t temp_mpfr;
	mpc_t temp_mpc;

			
		
	
	
	
	
	for (int j=0;j<2;j++) 	{

		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		

		
		if(j==0) out_len = 10;
		if(j==1) 
			
		{
			if ((exponent - max_err_exp)>=(20 - exp_minG))
				{err_exp = max_err_exp;}
			else
				{err_exp = exponent + exp_minG - 20;}
			out_len = exponent - err_exp;
//			cout<<"out_len = "<<out_len<<"\n";
//			cout<<"err_exp = "<<err_exp<<"\n";
		}


		Mstr = getM(c_field,c_energy, c_Z1);

		
		if (line != Mstr)	{
			line.clear();
			line = Mstr;		
//			cout<<out_len<<" "<<nsum1<<"  "<<precision<<"  "<<line<<"\n";
			count = 0;
			
		} else {
			
			_temp = sum;
			sum = prec;
			prec = _temp;
			count++;
		}
			
		if (count == 2)
		break;	


	
	}

		mpfr_init2(temp_mpfr, prec_bit);
		mpc_init2(temp_mpc, prec_bit);
	
		mpc_set_str(temp_mpc, &Mstr[0], 10, MPC_RNDNN);
		mpc_abs(temp_mpfr,temp_mpc,GMP_RNDN);
		mpfr_log10(temp_mpfr, temp_mpfr,GMP_RNDN);
		
		exponent = (int)mpfr_get_d(temp_mpfr,GMP_RNDN);
		if (exponent<5) exponent = 5;
//		if (j==1) {cout<<"exponent= "<<exponent<<"\n";cout<<"nsum1= "<<nsum1<<"\n";}
		
		mpfr_clear(temp_mpfr);
		mpc_clear(temp_mpc);

		
	}	


	
//	out_len = precision;
	
	return prec_bit + 10;
}










std::string Functions::getM(std::string c_F, std::string c_energy, std::string c_Z1){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z1;
	mpc_t sum_M;
	mpc_t M;
	mpc_t temp;
	mpc_t temp1;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	
	
	std::string str;
	
	mpfr_init2(F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z1, prec_bit);
	mpc_init2(sum_M, prec_bit);
	mpc_init2(M, prec_bit);
	mpc_init2(temp, prec_bit);
	mpc_init2(temp1, prec_bit);
	mpc_init2(Ci_3, prec_bit);
	mpc_init2(Ci_2, prec_bit);
	mpc_init2(Ci_1, prec_bit);
	mpc_init2(Ci, prec_bit);
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);


	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1, 1, MPC_RNDNN);
	


	mpc_set_ui(sum_M,1,MPC_RNDNN);	
	for (int i=1;i<nsum1+1;i++)	{
				
		mpc_mul_fr(Ci,Ci_3,F,MPC_RNDNN); 
		mpc_mul_ui(Ci,Ci,point2,MPC_RNDNN); 
		mpc_mul_ui(Ci,Ci,point2,MPC_RNDNN);
		mpc_mul(temp,Ci_1,Z1,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_mul(temp,Ci_2,E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2*point2,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_mul_ui(Ci,Ci,point2,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,4,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i + m,MPC_RNDNN);

		mpc_add(sum_M,sum_M,Ci,MPC_RNDNN);

			
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

		
	}

	
	mpc_set_ui(temp,point1, MPC_RNDNN);
	mpc_set_ui(temp1,m, MPC_RNDNN);
	mpc_pow(temp, temp,temp1, MPC_RNDNN);
	mpc_mul(M,sum_M,temp,MPC_RNDNN);
	

		

	char* gg = mpc_get_str(10, out_len, M, MPC_RNDNN);
	str = gg;
	delete gg;

	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(sum_M);
	mpc_clear(M);
	mpc_clear(temp);
	mpc_clear(temp1);
	mpc_clear(Ci_3);
	mpc_clear(Ci_2);
	mpc_clear(Ci_1);
	mpc_clear(Ci);
	
	
	return str;

}




long int Functions::calcPrecisionForN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	long int crit;
	


	
	while(true)	{
		crit =  - out_lenN + precisionN - getN(str_N,  str_derN, c_F,  c_energy,  c_Z2);
		if (crit>=500)	{pointN1++; pointN2 = pointN1*pointN1;}
		if (crit<=300)	{pointN1--; pointN2 = pointN1*pointN1;}
		if ((crit>300)&&(crit<500))	break;
	}
//	cout<<"pointN1 = "<<pointN1<<"\n";
//	cout<<"crit = "<<crit<<"\n";
//	return (int)(out_lenN*3.3219280948873626);
	return pointN1;
}



long int Functions::getN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z2)	{


	std::string temp_str_N = "empty";
	std::string temp_str_der_N = "empty";
	
	
	mpfr_t F;
	mpfr_t Ci_max;
	mpfr_t absCi;
	mpfr_t abs_sum_N;
	mpc_t E;
	mpc_t Z2;
	mpc_t N;
	mpc_t der_N;
	mpc_t sum_der_N;
	mpc_t sum_N;
	mpc_t temp;
	mpc_t power;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	



	
	mpfr_init2(F, prec_bitN);
	mpfr_init2(abs_sum_N,prec_bitN);
	mpfr_init2(Ci_max, prec_bitN);
	mpfr_init2(absCi, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(N, prec_bitN);
	mpc_init2(der_N, prec_bitN);
	mpc_init2(sum_der_N, prec_bitN);
	mpc_init2(sum_N, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpc_init2(power, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);

	mpfr_neg(F,F,GMP_RNDN);	
	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1, 1, MPC_RNDNN);

		


	mpfr_set_ui(Ci_max,1,GMP_RNDN);
	mpc_set_ui(sum_N,1,MPC_RNDNN);	
	mpc_set_ui(sum_der_N,2*m + 1,MPC_RNDNN);
	mpc_div_ui(sum_der_N,sum_der_N,2,MPC_RNDNN);
		for (long int i=1;i<1000000000;i++)	{
			
			mpc_mul_fr(Ci,Ci_3,F,MPC_RNDNN); 
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN); 
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN);
			mpc_mul(temp,Ci_1,Z2,MPC_RNDNN);
			mpc_sub(Ci,Ci,temp,MPC_RNDNN);
			mpc_mul(temp,Ci_2,E,MPC_RNDNN); 
			mpc_mul_ui(temp,temp,2*pointN2,MPC_RNDNN);
			mpc_sub(Ci,Ci,temp,MPC_RNDNN);
			mpc_mul_ui(Ci,Ci,pointN2,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,4,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,i,MPC_RNDNN);
			mpc_div_ui(Ci,Ci,i + m,MPC_RNDNN);
			
			mpc_abs(absCi,Ci,GMP_RNDN);
			if(mpfr_greater_p(absCi,Ci_max) != 0)
				mpfr_set(Ci_max,absCi,GMP_RNDN);
			
			mpc_set_ui(temp,4,MPC_RNDNN);
			mpc_mul_ui(temp,temp,i,MPC_RNDNN);
			mpc_add_ui(temp,temp,2*m + 1,MPC_RNDNN);
			mpc_div_ui(temp,temp,2,MPC_RNDNN);
			mpc_mul(temp,temp,Ci,MPC_RNDNN);
			
			
			mpc_add(sum_der_N,sum_der_N,temp,MPC_RNDNN);
			mpc_add(sum_N,sum_N,Ci,MPC_RNDNN);
			
		
			
			if (i%100==0)	{
//				std::cout<<"pointN1 = "<<pointN1<<" i= "<<i<<" "<<"sum_der_N = "<<mpc_get_str(10, out_lenN, sum_der_N, MPC_RNDNN)<<"  sum = "<<mpc_get_str(10, out_lenN, sum_N, MPC_RNDNN)<<"\n";
			
			char* gg = mpc_get_str(10, out_lenN, sum_N, MPC_RNDNN); str_N = gg; delete gg;
			char* hh = mpc_get_str(10, out_lenN, sum_der_N, MPC_RNDNN); str_der_N = hh; delete hh;
			
			if((str_N != temp_str_N)&&(temp_str_der_N != str_der_N))
			{
				temp_str_N = str_N;
				temp_str_der_N = temp_str_der_N;
			}
			else
				break;
			}
			
			
			
			
			mpc_set(Ci_3,Ci_2,MPC_RNDNN);
			mpc_set(Ci_2,Ci_1,MPC_RNDNN);
			mpc_set(Ci_1,Ci,MPC_RNDNN);
			
		}
		
		mpc_set_ui(power,1 + 2*m, MPC_RNDNN);
		mpc_div_ui(power,power,2, MPC_RNDNN);
		mpc_set_ui(temp,pointN1, MPC_RNDNN);
		mpc_pow(temp, temp, power, MPC_RNDNN);
		mpc_mul(N, sum_N, temp, MPC_RNDNN);
		
		mpc_mul(der_N,sum_der_N,temp, MPC_RNDNN);
		mpc_div_ui(der_N,der_N,pointN1, MPC_RNDNN);
		

		
		char* gg = mpc_get_str(10, out_lenN, N, MPC_RNDNN); str_N = gg; delete gg;
		char* hh = mpc_get_str(10, out_lenN, der_N, MPC_RNDNN); str_der_N = hh; delete hh;
		
		
		mpfr_log10(Ci_max,Ci_max,GMP_RNDN);
		mpc_abs(abs_sum_N,sum_N,GMP_RNDN);
		mpfr_log10(abs_sum_N,abs_sum_N,GMP_RNDN);
		long int poww = (long int)mpfr_get_d(Ci_max,GMP_RNDN) - (long int)mpfr_get_d(abs_sum_N,GMP_RNDN);
		
		mpfr_clear(F);
		mpfr_clear(Ci_max);
		mpfr_clear(abs_sum_N);
		mpfr_clear(absCi);
		mpc_clear(E);
		mpc_clear(Z2);
		mpc_clear(N);
		mpc_clear(der_N);
		mpc_clear(sum_der_N);
		mpc_clear(sum_N);
		mpc_clear(temp);
		mpc_clear(power);
		mpc_clear(Ci_3);
		mpc_clear(Ci_2);
		mpc_clear(Ci_1);
		mpc_clear(Ci);
	

	
	return poww;
}





long int Functions::get_a(std::string& str_a, std::string& str_der_a,std::string c_F, std::string c_energy, std::string c_Z2){
	
	
	std::string temp_str_a = "empty";
	std::string temp_str_der_a = "empty";
	long int i;
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z2;
	mpc_t a;
	mpc_t der_a;
	mpc_t temp;
	mpfr_t temp1;
	mpc_t temp_exp;
	mpc_t sum_a;
	mpc_t sum_der_a;
	
	mpc_t Ci;
	mpc_t Ci_1;
	mpc_t Ci_2;
	mpc_t Ci_3;
	
	mpc_t k1;
	mpc_t k2;
	mpc_t k3;

	
	mpc_init2(k1, prec_bitN);
	mpc_init2(k2, prec_bitN);
	mpc_init2(k3, prec_bitN);
	mpfr_init2(F, prec_bitN);
	mpfr_init2(sqrt_F, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(a, prec_bitN);
	mpc_init2(der_a, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpfr_init2(temp1, prec_bitN);
	mpc_init2(temp_exp, prec_bitN);
	mpc_init2(sum_a, prec_bitN);
	mpc_init2(sum_der_a, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);


	
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	mpc_set_fr(temp,F,MPC_RNDNN);

	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4*pointN2,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);	
	mpc_mul_ui(k2,k2,pointN1,MPC_RNDNN);
	mpc_mul_ui(k2,k2,pointN2,MPC_RNDNN);
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);
	mpc_neg(k2,k2,MPC_RNDNN);
	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8*pointN1,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	

	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1,1, MPC_RNDNN);
	mpc_div_ui(Ci_1, Ci_1, pointN1, MPC_RNDNN);

	

	mpc_set(sum_a,Ci_1,MPC_RNDNN);	
	mpc_set(sum_der_a,Ci_1,MPC_RNDNN);
	for (i=2;i<nsumab;i++)	{
		
		mpc_set_si(Ci,2*i-5, MPC_RNDNN);
		mpc_mul_si(Ci,Ci,5-2*i, MPC_RNDNN);
		mpc_add_ui(Ci,Ci,4*m*m, MPC_RNDNN);
		mpc_mul(Ci,Ci,Ci_3,MPC_RNDNN);
		mpc_mul_si(temp,k3,i - 2,MPC_RNDNN);
		mpc_mul(temp,temp,Ci_2,MPC_RNDNN);
		mpc_add(Ci,Ci, temp,MPC_RNDNN);
		mpc_mul(temp, Ci_1,k1,MPC_RNDNN);
		mpc_add(Ci,Ci,temp,MPC_RNDNN);
		mpc_div(Ci,Ci,k2,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i - 1,MPC_RNDNN);
		
		mpc_mul_ui(temp,Ci,i,MPC_RNDNN);
		
		mpc_add(sum_der_a,sum_der_a,temp,MPC_RNDNN);
		mpc_add(sum_a,sum_a,Ci,MPC_RNDNN);
		
		
		
		if (i%10==0)	{
//			std::cout<<"pointN1 = "<<pointN1<<" i= "<<i<<" "<<"sum_der_a = "<<mpc_get_str(10, out_lenN, sum_der_a, MPC_RNDNN)<<"  sum_a = "<<mpc_get_str(10, out_lenN, sum_a, MPC_RNDNN)<<"\n";
//			std::cout<<i<<"\n";
		char* gg  = mpc_get_str(10, out_lenN, sum_a, MPC_RNDNN); str_a = gg; delete gg;
		char* hh = mpc_get_str(10, out_lenN, sum_der_a, MPC_RNDNN); str_der_a = hh; delete hh;
		
		if((str_a != temp_str_a)&&(temp_str_der_a != str_der_a))
		{
			temp_str_a = str_a;
			temp_str_der_a = temp_str_der_a;
		}
		else
			break;
		}
			
		
		

		
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

	}
	

	
	mpfr_mul_ui(temp1,sqrt_F,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,pointN2,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_ui(temp,E,pointN1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,1,MPC_RNDNN);
	mpc_exp(temp_exp,temp,MPC_RNDNN);
	

	mpc_mul(a,sum_a,temp_exp,MPC_RNDNN);
	
	
		
	mpfr_mul_ui(temp1,sqrt_F,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,pointN2,GMP_RNDN);
	mpc_mul_ui(temp,E,pointN1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,1,MPC_RNDNN);

	mpc_mul(der_a,sum_a,temp,MPC_RNDNN);
	mpc_sub(der_a,der_a,sum_der_a,MPC_RNDNN);
	mpc_mul(der_a,der_a,temp_exp,MPC_RNDNN);
	mpc_div_ui(der_a,der_a,pointN1,MPC_RNDNN);

	
	char* gg = mpc_get_str(10, out_lenN, a, MPC_RNDNN); str_a = gg; delete gg;
	char* hh = mpc_get_str(10, out_lenN, der_a, MPC_RNDNN); str_der_a = hh; delete hh;
	

	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(a);
	mpc_clear(der_a);
	mpc_clear(temp);
	mpc_clear(temp_exp);
	mpfr_clear(temp1);
	mpc_clear(sum_a);
	mpc_clear(sum_der_a);
	mpc_clear(Ci);
	mpc_clear(Ci_1);
	mpc_clear(Ci_2);
	mpc_clear(Ci_3);

	return i;
	
}




long int Functions::get_b(std::string& str_b, std::string& str_der_b,std::string c_F, std::string c_energy, std::string c_Z2){
	
	
	std::string temp_str_b = "empty";
	std::string temp_str_der_b = "empty";
	long int i;
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z2;
	mpc_t b;
	mpc_t der_b;
	mpc_t temp;
	mpfr_t temp1;
	mpc_t temp_exp;
	mpc_t sum_b;
	mpc_t sum_der_b;
	
	mpc_t Ci;
	mpc_t Ci_1;
	mpc_t Ci_2;
	mpc_t Ci_3;
	
	mpc_t k1;
	mpc_t k2;
	mpc_t k3;

	
	mpc_init2(k1, prec_bitN);
	mpc_init2(k2, prec_bitN);
	mpc_init2(k3, prec_bitN);
	mpfr_init2(F, prec_bitN);
	mpfr_init2(sqrt_F, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(b, prec_bitN);
	mpc_init2(der_b, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpfr_init2(temp1, prec_bitN);
	mpc_init2(temp_exp, prec_bitN);
	mpc_init2(sum_b, prec_bitN);
	mpc_init2(sum_der_b, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z2, &c_Z2[0], 10, MPC_RNDNN);


	
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	mpc_set_fr(temp,F,MPC_RNDNN);

	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4*pointN2,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);	
	mpc_mul_ui(k2,k2,pointN1,MPC_RNDNN);
	mpc_mul_ui(k2,k2,pointN2,MPC_RNDNN);
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);
	
//	mpc_neg(k2,k2,MPC_RNDNN);

	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8*pointN1,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	

	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1,1, MPC_RNDNN);
	mpc_div_ui(Ci_1, Ci_1, pointN1, MPC_RNDNN);

	

	mpc_set(sum_b,Ci_1,MPC_RNDNN);	
	mpc_set(sum_der_b,Ci_1,MPC_RNDNN);
	for (i=2;i<nsumab;i++)	{
		
		mpc_set_si(Ci,2*i-5, MPC_RNDNN); 
		mpc_mul_si(Ci,Ci,5-2*i, MPC_RNDNN);
		mpc_add_ui(Ci,Ci,4*m*m, MPC_RNDNN);
		mpc_mul(Ci,Ci,Ci_3,MPC_RNDNN);
		mpc_mul_si(temp,k3,2 - i,MPC_RNDNN); 
		mpc_mul(temp,temp,Ci_2,MPC_RNDNN);
		mpc_add(Ci,Ci, temp,MPC_RNDNN);
		mpc_mul(temp, Ci_1,k1,MPC_RNDNN);
		mpc_add(Ci,Ci,temp,MPC_RNDNN);
		mpc_div(Ci,Ci,k2,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i - 1,MPC_RNDNN);
		
		mpc_mul_ui(temp,Ci,i,MPC_RNDNN);
		
		mpc_add(sum_der_b,sum_der_b,temp,MPC_RNDNN);
		mpc_add(sum_b,sum_b,Ci,MPC_RNDNN);
		
//		std::cout<<"C"<<i<<" = "<<mpc_get_str(10, out_lenN, Ci, MPC_RNDNN)<<"\n";
	
		
		if (i%10==0)	{
//			std::cout<<"pointN1 = "<<pointN1<<" i= "<<i<<" "<<"sum_der_b = "<<mpc_get_str(10, out_lenN, sum_der_b, MPC_RNDNN)<<"  sum_b = "<<mpc_get_str(10, out_lenN, sum_b, MPC_RNDNN)<<"\n";
		
		char* gg = mpc_get_str(10, out_lenN, sum_b, MPC_RNDNN); str_b = gg; delete gg;
		char* hh = mpc_get_str(10, out_lenN, sum_der_b, MPC_RNDNN); str_der_b = hh; delete hh;
		
		if((str_b != temp_str_b)&&(temp_str_der_b != str_der_b))
		{
			temp_str_b = str_b;
			temp_str_der_b = temp_str_der_b;
		}
		else
			break;
		}
			
		
		

		
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

	}
	

	
	mpfr_mul_ui(temp1,sqrt_F,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,pointN2,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_ui(temp,E,pointN1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,-1,MPC_RNDNN);
	mpc_exp(temp_exp,temp,MPC_RNDNN);
	

	mpc_mul(b,sum_b,temp_exp,MPC_RNDNN);
	
	
		
	mpfr_mul_ui(temp1,sqrt_F,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp1,temp1,pointN2,GMP_RNDN);
	mpc_mul_ui(temp,E,pointN1,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,-1,MPC_RNDNN);

	mpc_mul(der_b,sum_b,temp,MPC_RNDNN);
	mpc_sub(der_b,der_b,sum_der_b,MPC_RNDNN);
	mpc_mul(der_b,der_b,temp_exp,MPC_RNDNN);
	mpc_div_ui(der_b,der_b,pointN1,MPC_RNDNN);

	
	char* gg = mpc_get_str(10, out_lenN, b, MPC_RNDNN);str_b = gg; delete gg;
	char* hh = mpc_get_str(10, out_lenN, der_b, MPC_RNDNN);str_der_b = hh; delete hh;
	

	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(b);
	mpc_clear(der_b);
	mpc_clear(temp);
	mpc_clear(temp_exp);
	mpfr_clear(temp1);
	mpc_clear(sum_b);
	mpc_clear(sum_der_b);
	mpc_clear(Ci);
	mpc_clear(Ci_1);
	mpc_clear(Ci_2);
	mpc_clear(Ci_3);
	
	return i;
	
}








std::string Functions::getB(std::string c_F,std::string c_energy, std::string c_Z2){
	
	std::string str_B;
	std::string str_N;
	std::string str_der_N;
	std::string str_a;
	std::string str_der_a;
	std::string str_b;
	std::string str_der_b;

	long int na, nb;
		
		mpc_t a;
		mpc_t der_a;
		mpc_t b;
		mpc_t der_b;
		mpc_t N;
		mpc_t der_N;
		mpc_t B;
		mpc_t temp1;
		mpc_t temp2;
	
		mpc_init2(B, prec_bitN);
		mpc_init2(a, prec_bitN);
		mpc_init2(der_a, prec_bitN);
		mpc_init2(b, prec_bitN);
		mpc_init2(der_b, prec_bitN);
		mpc_init2(N, prec_bitN);
		mpc_init2(der_N, prec_bitN);
		mpc_init2(temp1, prec_bitN);
		mpc_init2(temp2, prec_bitN);


	na = get_a(str_a,  str_der_a, c_F,  c_energy, c_Z2);	
	nb = get_b(str_b,  str_der_b, c_F,  c_energy, c_Z2);	
	getN(str_N,str_der_N,c_F,c_energy, c_Z2);	
	
	
	mpc_set_str(a, &str_a[0], 10, MPC_RNDNN);
	mpc_set_str(der_a, &str_der_a[0], 10, MPC_RNDNN);
	mpc_set_str(b, &str_b[0], 10, MPC_RNDNN);
	mpc_set_str(der_b, &str_der_b[0], 10, MPC_RNDNN);
	mpc_set_str(N, &str_N[0], 10, MPC_RNDNN);
	mpc_set_str(der_N, &str_der_N[0], 10, MPC_RNDNN);
	mpc_mul(temp1,N,der_a,MPC_RNDNN);
	mpc_mul(temp2,a,der_N,MPC_RNDNN);
	mpc_sub(B,temp1,temp2,MPC_RNDNN);
	mpc_mul(temp1,b,der_a,MPC_RNDNN);
	mpc_mul(temp2,a,der_b,MPC_RNDNN);
	mpc_sub(temp1,temp1,temp2,MPC_RNDNN);
	mpc_div(B,B,temp1,MPC_RNDNN);

	if((na==nsumab)||(nb==nsumab))
	str_B="(0 0)";
	else{
	char* gg = mpc_get_str(10, out_lenN, B, MPC_RNDNN); str_B = gg; delete gg;
	}
	
	mpc_clear(a);
	mpc_clear(der_a);
	mpc_clear(b);
	mpc_clear(der_b);
	mpc_clear(N);
	mpc_clear(der_N);
	mpc_clear(B);
	mpc_clear(temp1);
	mpc_clear(temp2);
	

	
	return str_B;
}



