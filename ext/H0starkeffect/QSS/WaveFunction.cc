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
#include "WaveFunction.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>



#define MAX_DIV_ARRAY 100000
#define MIN_DIV_ARRAY 10000




using namespace OrbitUtils;
using namespace std;




WaveFunction::WaveFunction(int n11,int n22, int mm,  int point11,  int pointN,  int precEZ, std::string c_energy, std::string c_Z1, std::string c_F)	{		
	




	
	std::string str_A;
	std::string str_B;
	std::string str_N;
	std::string str_der_N;
	std::string str_a;
	std::string str_der_a;
	std::string str_b;
	std::string str_der_b;
	
	int nsuma, nsumb;
	
	mpfr_t sqrtF;
	mpfr_t F;
	

	mode = 3;
	
	n1 = n11;
	n2 = n22;
	m = abs(mm);
	n = n1 + n2 + m + 1; 
	
	
	getArrayb();
	A = sqrt(2*MathPolynomial::FactorialDouble(n1+m)*MathPolynomial::FactorialDouble(n2+m)/(MathPolynomial::FactorialDouble(n1)*MathPolynomial::FactorialDouble(n2)))/(n*n);



	
	precisionN = 1000;
	prec_bitN = (long int)(precisionN*3.3219280948873626);
	
	precision = precEZ;
	prec_bit = (long int)(precision*3.3219280948873626);
	
	out_lenN = 100;
	out_len = 100;	
	
	mpfr_init2(h, prec_bitN);
	mpfr_init2(F, prec_bitN);
	mpfr_init2(sqrtF, prec_bitN);
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpfr_sqrt(sqrtF,F,GMP_RNDN);
	Fd = mpfr_get_d(F,GMP_RNDN);
	
	
	point1 = point11;
	point2 = point1*point1;
			
	pointN1 = pointN;
	pointN2 = pointN1*pointN1;
	
		
	
	
	if (mpfr_cmp_ui(F,0) == 0)
		{mode = 0;}

		else	{
		mpfr_ui_div(h,2,sqrtF,GMP_RNDN);
		mpfr_add_ui(h,h,point1*point1,GMP_RNDN);
		mpfr_sqrt(h,h,GMP_RNDN);
//		mpfr_set_ui(h,50,GMP_RNDN);
		
		if(mpfr_cmp_ui(h,pointN) != 1)	{
			mode = 1;
		if(mpfr_cmp_ui(h,pointN) == -1)	{
			pointN1 = mpfr_get_ui(h,GMP_RNDZ)+1;
			pointN2 = pointN1*pointN1;
		}
		}
		else {mode = 2; mpfr_set_ui(h,mpfr_get_ui(h,GMP_RNDZ)+1,GMP_RNDN);}
		
		}
		

//	mode = 0;
	
	if(mode != 0)	{
		double _ndivN = ndivN_b(c_energy, c_F);

		
		if (_ndivN<MAX_DIV_ARRAY)		
			ndivN = (int)_ndivN;	
		else 
			mode = 0;

		
		if (ndivN<MIN_DIV_ARRAY)
			ndivN = MIN_DIV_ARRAY;

		ndivM = MIN_DIV_ARRAY;
	}
	
	
	

	if(mode == 2){
		double _ndivab = ndivN_h(c_energy, c_F);


		if (_ndivab<MAX_DIV_ARRAY)
			ndivab = (int)_ndivab;
		else
			mode = 0;
		
		
		if (ndivab<MIN_DIV_ARRAY)
			ndivab = MIN_DIV_ARRAY;
	}
	

	
	
	if(mode == 2){
		
		nsuma = get_a(str_a, str_der_a,c_F,c_energy, c_Z1);
		nsumb = get_b(str_b, str_der_b,c_F,c_energy, c_Z1);
		
		if(nsuma>nsumb)	 nsumab = nsuma+100;
		if(nsuma<nsumb)	 nsumab = nsumb+100;
		if(nsuma==nsumb) nsumab = nsuma+100;
		
		if((nsumab-100)>(100000-10))	
			mode = 0;		
	}
	
	
	
	
// mode = 0,1 or 2 is completely defined	
/*	
	cout<<mpfr_get_d(h,GMP_RNDN)<<"\n";
	cout<<"nsuma = "<<nsuma<<"\n";
	cout<<"nsumb = "<<nsumb<<"\n";
	cout<<"ndivab = "<<ndivab<<"\n";
	cout<<"ndivN = "<<ndivN<<"\n";
	cout<<"ndivM = "<<ndivM<<"\n";
	cout<<"pointN1 = "<<pointN1<<"\n";
	cout<<"mode = "<<mode<<"\n";
*/	
	
	
	
	if(mode != 0)	{
		
	
	nsum1 = getnsumM(c_F,c_energy,c_Z1) + 10;
//	cout<<"nsum1 = "<<nsum1<<"\n";
	
	
	c = new mpc_t[nsum1 + 1];
			for (int i=0;i<nsum1+1;i++)	
				mpc_init2(c[i], prec_bit);

	
	Mi = new mpc_t[ndivM + 1];
		for (int i=0;i<ndivM+1;i++)	
			mpc_init2(Mi[i], prec_bit);
		
	Ma = new mpc_t[ndivM + 1];
	Mb = new mpc_t[ndivM + 1];
	Mc = new mpc_t[ndivM + 1];
	Md = new mpc_t[ndivM + 1];
	
		for (int i=0;i<ndivM+1;i++)	{	
				mpc_init2(Ma[i], prec_bit);
				mpc_init2(Mb[i], prec_bit);
				mpc_init2(Mc[i], prec_bit);
				mpc_init2(Md[i], prec_bit);
			}
			
	
		
	getArrayM(c_F,c_energy,c_Z1);
	fillArrayM();
	
	for (int i=0;i<nsum1+1;i++)	
		mpc_clear(c[i]);
	delete [] c;	
	
	for (int i=0;i<ndivM+1;i++)	
		mpc_clear(Mi[i]);
	delete [] Mi;
	
	

	
	
	
	
	
	nsumN = getnsumN(str_N,str_der_N,c_F,c_energy,c_Z1) + 10;
//	cout<<"Nsum = "<<nsumN<<"\n";
	
	c = new mpc_t[nsumN + 1];
			for (int i=0;i<nsumN+1;i++)	
				mpc_init2(c[i], prec_bitN);
	

			
			Ni = new mpc_t[ndivN + 1];
				for (int i=0;i<ndivN+1;i++)	
					mpc_init2(Ni[i], prec_bitN);
				
			Na = new mpc_t[ndivN + 1];
			Nb = new mpc_t[ndivN + 1];
			Nc = new mpc_t[ndivN + 1];
			Nd = new mpc_t[ndivN + 1];
			
				for (int i=0;i<ndivN+1;i++)	{	
						mpc_init2(Na[i], prec_bitN);
						mpc_init2(Nb[i], prec_bitN);
						mpc_init2(Nc[i], prec_bitN);
						mpc_init2(Nd[i], prec_bitN);
					}
		
	getArrayN(c_F,c_energy,c_Z1);
	fillArrayN();
			

	
	for (int i=0;i<nsumN+1;i++)	
		mpc_clear(c[i]);
	delete [] c;	
	
	for (int i=0;i<ndivN+1;i++)	
		mpc_clear(Ni[i]);
	delete [] Ni;
	


if (mode == 2)	{	
	
	abi = new mpc_t[ndivab + 1];
		for (int i=0;i<ndivab+1;i++)	
			mpc_init2(abi[i], prec_bitN);
		
	ab_a = new mpc_t[ndivab + 1];
	ab_b = new mpc_t[ndivab + 1];
	ab_c = new mpc_t[ndivab + 1];
	ab_d = new mpc_t[ndivab + 1];
	
		for (int i=0;i<ndivab+1;i++)	{	
				mpc_init2(ab_a[i], prec_bitN);
				mpc_init2(ab_b[i], prec_bitN);
				mpc_init2(ab_c[i], prec_bitN);
				mpc_init2(ab_d[i], prec_bitN);
			}
			
	
	
	

	
	getAB(str_a, str_der_a, str_b, str_der_b, str_N, str_der_N, str_A, str_B);
	

	
	a = new mpc_t[nsumab + 1];
			for (int i=0;i<nsumab+1;i++)	
				mpc_init2(a[i], prec_bitN);
			
	b = new mpc_t[nsumab + 1];
			for (int i=0;i<nsumab+1;i++)	
				mpc_init2(b[i], prec_bitN);
			
			
	getArray_a(str_A,c_F,c_energy,c_Z1);
	getArray_b(str_B,c_F,c_energy,c_Z1);
	fillArray_ab(c_F, c_energy);
			
	
	
			
	for (int i=0;i<nsumab+1;i++)	
		mpc_clear(a[i]);
	delete [] a;
	
	for (int i=0;i<nsumab+1;i++)	
		mpc_clear(b[i]);
	delete [] b;	
	
	for (int i=0;i<ndivab+1;i++)	
		mpc_clear(abi[i]);
	delete [] abi;
	
}
	}
	
	
	mpfr_clear(F);
	mpfr_clear(sqrtF);
	
}





WaveFunction::~WaveFunction()	{
	
	mpfr_clear(h);
	for (int i=0;i<15;i++) for (int j=0;j<15;j++)	delete [] bp[i][j]; for (int i=0;i<15;i++)	delete [] bp[i];	delete	[]	bp;
	
	if (mode != 0 )	{
	
	for (int i=0;i<ndivM+1;i++)	{	
		mpc_clear(Ma[i]);
		mpc_clear(Mb[i]);
		mpc_clear(Mc[i]);
		mpc_clear(Md[i]);
	}

	delete [] Ma;
	delete [] Mb;
	delete [] Mc;
	delete [] Md;
	
	
	
	for (int i=0;i<ndivN+1;i++)	{	
		mpc_clear(Na[i]);
		mpc_clear(Nb[i]);
		mpc_clear(Nc[i]);
		mpc_clear(Nd[i]);
	}

	delete [] Na;
	delete [] Nb;
	delete [] Nc;
	delete [] Nd;
	
if (mode == 2)		{
	for (int i=0;i<ndivab+1;i++)	{	
		mpc_clear(ab_a[i]);
		mpc_clear(ab_b[i]);
		mpc_clear(ab_c[i]);
		mpc_clear(ab_d[i]);
	}

	delete [] ab_a;
	delete [] ab_b;
	delete [] ab_c;
	delete [] ab_d;
}

	}
	
}






double WaveFunction::ndivN_b(std::string c_energy, std::string c_F){
	
	mpfr_t F;
	mpfr_t sqrtF;
	mpfr_t temp;
	mpfr_t temp2;

	mpc_t E;
	
	mpfr_init2(F, prec_bitN);
	mpfr_init2(temp, prec_bitN);
	mpfr_init2(temp2, prec_bitN);
	mpfr_init2(sqrtF, prec_bitN);
	mpc_init2(E, prec_bitN);
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	
	mpfr_sqrt(sqrtF,F,GMP_RNDN);
	mpfr_mul_ui(temp,sqrtF,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp,temp,pointN1,GMP_RNDN);
	
	mpc_real(temp2,E,GMP_RNDN);
	mpfr_div(temp2,temp2,sqrtF,GMP_RNDN);
	mpfr_add(temp,temp,temp2,GMP_RNDN);
	mpfr_mul_ui(temp,temp,2*pointN1,GMP_RNDN);
	double number = mpfr_get_d(temp,GMP_RNDN);
	
	
	
	mpfr_clear(F);
	mpfr_clear(sqrtF);
	mpfr_clear(temp);
	mpfr_clear(temp2);
	mpc_clear(E);
	
	return number;
}






double WaveFunction::ndivN_h(std::string c_energy, std::string c_F){
	

	
	mpfr_t F;
	mpfr_t sqrtF;
	mpfr_t temp;
	mpfr_t temp2;

	mpc_t E;
	
	mpfr_init2(F, prec_bitN);
	mpfr_init2(temp, prec_bitN);
	mpfr_init2(temp2, prec_bitN);
	mpfr_init2(sqrtF, prec_bitN);
	mpc_init2(E, prec_bitN);
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpfr_sqrt(sqrtF,F,GMP_RNDN);


	
	mpfr_sqrt(sqrtF,F,GMP_RNDN);
	mpfr_mul(temp,sqrtF,h,GMP_RNDN);
	mpfr_mul(temp,temp,h,GMP_RNDN);
	
	mpc_real(temp2,E,GMP_RNDN);
	mpfr_div(temp2,temp2,sqrtF,GMP_RNDN);
	mpfr_add(temp,temp,temp2,GMP_RNDN);
	mpfr_sub_ui(temp2,h,pointN1,GMP_RNDN);
	mpfr_mul(temp,temp,temp2,GMP_RNDN);
	mpfr_mul_ui(temp,temp,2,GMP_RNDN);
	double number = mpfr_get_d(temp,GMP_RNDN);

	
	
	mpfr_clear(F);
	mpfr_clear(sqrtF);
	mpfr_clear(temp);
	mpfr_clear(temp2);
	mpc_clear(E);
	
	return number;
}





long int WaveFunction::getnsumM(std::string c_F, std::string c_energy, std::string c_Z1){
	
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
	int k;
	
	

	std::string temp_str_M = "empty";
	std::string str_M = "empty";
	
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
	for (long int i=1;i<1000000000;i++)	{
				
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
		
		
		
		if (i%10==0)	{				
		char* gg = mpc_get_str(10, out_len, sum_M, MPC_RNDNN); str_M = gg; delete gg;
		k = i;
		if(str_M != temp_str_M)
		{
			temp_str_M = str_M;

		}
		else
			
			break;
	
		}
		

			
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

		
	}

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
	
	
	return k+10;

}


void WaveFunction::getArrayM(std::string c_F, std::string c_energy, std::string c_Z1){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z1;
	mpc_t temp;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	
	
	
	mpfr_init2(F, prec_bit);
	mpc_init2(E, prec_bit);
	mpc_init2(Z1, prec_bit);
	mpc_init2(temp, prec_bit);
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

	mpc_set_ui(c[0],1,MPC_RNDNN);
	

	
	for (long int i=1;i<nsum1+1;i++)	{
		
		mpc_mul_fr(Ci,Ci_3,F,MPC_RNDNN); 
		mpc_mul(temp,Ci_1,Z1,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_mul(temp,Ci_2,E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,4,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i + m,MPC_RNDNN);

		mpc_set(c[i],Ci,MPC_RNDNN);

			
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);
		
		

		
	}


	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(temp);
	mpc_clear(Ci_3);
	mpc_clear(Ci_2);
	mpc_clear(Ci_1);
	mpc_clear(Ci);
	
	
	return;

}


void WaveFunction::getM(mpfr_t mu, mpc_t M){
	


	mpfr_t temp1;
	mpfr_t mu2;
	mpfr_t mu2i;
	mpc_t sum_M;
	mpc_t temp;
	

	std::string str;
	std::string temp_str_M = "empty";
	std::string str_M = "empty";
	


	mpfr_init2(temp1, prec_bit);
	mpfr_init2(mu2, prec_bitN);
	mpfr_init2(mu2i, prec_bitN);
	mpc_init2(sum_M, prec_bit);
	mpc_init2(temp, prec_bit);



	mpfr_mul(mu2,mu,mu,GMP_RNDN);
	mpfr_set_ui(mu2i,1,GMP_RNDN);




	mpc_set_ui(sum_M,1,MPC_RNDNN);	
	for (long int i=1;i<nsum1+1;i++)	{

		mpfr_mul(mu2i,mu2i,mu2,GMP_RNDN);
		mpc_mul_fr(temp,c[i],mu2i,MPC_RNDNN);
		
		if (i%10==0)	{				
		char* gg = mpc_get_str(10, out_len, sum_M, MPC_RNDNN); str_M = gg; delete gg;
		if(str_M != temp_str_M)
		{
			temp_str_M = str_M;

		}
		else
			
			break;
	
		}
		

		mpc_add(sum_M,sum_M,temp,MPC_RNDNN);
	}


	mpfr_pow_ui(temp1, mu,m, GMP_RNDN);
	mpc_mul_fr(M,sum_M,temp1,GMP_RNDN);


	mpc_clear(sum_M);
	mpc_clear(temp);
	mpfr_clear(temp1);
	mpfr_clear(mu2);
	mpfr_clear(mu2i);


	
	return;

}


void WaveFunction::fillArrayM(){
	
	mpfr_t mu;
	mpfr_t dmu;
	mpfr_t k1;
	mpfr_t k2;
	mpfr_t k3;
	mpfr_t k4;
	mpfr_t kt;
	mpc_t temp;
	
	mpfr_init2(mu, prec_bit);
	mpfr_init2(dmu, prec_bit);
	mpfr_init2(k1, prec_bit);
	mpfr_init2(k2, prec_bit);
	mpfr_init2(k3, prec_bit);
	mpfr_init2(k4, prec_bit);
	mpfr_init2(kt, prec_bit);
	mpc_init2(temp, prec_bit);
	
	
	mpfr_set_ui(dmu,point1,GMP_RNDN);
	mpfr_div_ui(dmu,dmu,ndivM,GMP_RNDN);
	
	for (long int i=0;i<ndivM+1;i++)	{	
		mpfr_mul_ui(mu,dmu,i,GMP_RNDN);
		getM(mu,Mi[i]);
//		cout<<"M  "<<i<<"\n";
	}
	


	
	for (long int i=1;i<ndivM-1;i++)	{
		
		mpc_sub(Ma[i],Mi[i+2],Mi[i-1],MPC_RNDNN);
		mpc_sub(temp,Mi[i],Mi[i+1],MPC_RNDNN);
		mpc_mul_ui(temp,temp,3,MPC_RNDNN);
		mpc_add(Ma[i],Ma[i],temp,MPC_RNDNN);
		mpc_div_fr(Ma[i],Ma[i],dmu,MPC_RNDNN);
		mpc_div_fr(Ma[i],Ma[i],dmu,MPC_RNDNN);
		mpc_div_fr(Ma[i],Ma[i],dmu,MPC_RNDNN);
		mpc_div_ui(Ma[i],Ma[i],6,MPC_RNDNN);
	
		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_add_ui(k1,k1,2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_add_ui(k2,k2,1,GMP_RNDN);

		
		mpc_mul_ui(Mb[i],Mi[i-1],1+i,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i],k1,MPC_RNDNN);
		mpc_sub(Mb[i],Mb[i],temp,MPC_RNDNN);	
		mpc_mul_fr(temp,Mi[i+1],k2,MPC_RNDNN);
		mpc_add(Mb[i],Mb[i],temp,MPC_RNDNN);	
		mpc_mul_ui(temp,Mi[i+2],i,MPC_RNDNN);
		mpc_sub(Mb[i],Mb[i],temp,MPC_RNDNN);
		mpc_div_fr(Mb[i],Mb[i],dmu,MPC_RNDNN);
		mpc_div_fr(Mb[i],Mb[i],dmu,MPC_RNDNN);
		mpc_div_ui(Mb[i],Mb[i],2,MPC_RNDNN);
		
		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_sub_ui(k1,k1,1,GMP_RNDN);
		
		mpfr_set_ui(k2,9,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k2,k2,kt,GMP_RNDN);
		mpfr_sub_ui(k2,k2,6,GMP_RNDN);
		
		mpfr_set_ui(k3,9,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_set_ui(kt,12,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k3,k3,kt,GMP_RNDN);
		mpfr_sub_ui(k3,k3,3,GMP_RNDN);
		
		mpfr_set_ui(k4,3,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k4,k4,kt,GMP_RNDN);
		mpfr_add_ui(k4,k4,2,GMP_RNDN);

		
		mpc_mul_fr(Mc[i],Mi[i+2],k1,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i+1],k2,MPC_RNDNN);
		mpc_sub(Mc[i],Mc[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i],k3,MPC_RNDNN);
		mpc_add(Mc[i],Mc[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i-1],k4,MPC_RNDNN);
		mpc_sub(Mc[i],Mc[i],temp,MPC_RNDNN);
		mpc_div_fr(Mc[i],Mc[i],dmu,MPC_RNDNN);
		mpc_div_ui(Mc[i],Mc[i],6,MPC_RNDNN);
		
		
		mpfr_set_ui(k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+1,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+2,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+1,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i-1,GMP_RNDN);
		
		mpfr_set_ui(k3,3,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i+2,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i-1,GMP_RNDN);
		
		mpfr_set_ui(k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i+1,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i-1,GMP_RNDN);

		
		
		mpc_mul_fr(Md[i],Mi[i-1],k1,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i],k2,MPC_RNDNN);
		mpc_sub(Md[i],Md[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i+1],k3,MPC_RNDNN);
		mpc_add(Md[i],Md[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Mi[i+2],k4,MPC_RNDNN);
		mpc_sub(Md[i],Md[i],temp,MPC_RNDNN);
		mpc_div_ui(Md[i],Md[i],6,MPC_RNDNN);
	
	}


	
	mpc_set(Ma[0],Ma[1],MPC_RNDNN);
	mpc_set(Mb[0],Mb[1],MPC_RNDNN);
	mpc_set(Mc[0],Mc[1],MPC_RNDNN);
	mpc_set(Md[0],Md[1],MPC_RNDNN);
	
	mpc_set(Ma[ndivM-1],Ma[ndivM-2],MPC_RNDNN);
	mpc_set(Mb[ndivM-1],Mb[ndivM-2],MPC_RNDNN);
	mpc_set(Mc[ndivM-1],Mc[ndivM-2],MPC_RNDNN);
	mpc_set(Md[ndivM-1],Md[ndivM-2],MPC_RNDNN);

	mpc_set(Ma[ndivM],Ma[ndivM-1],MPC_RNDNN);
	mpc_set(Mb[ndivM],Mb[ndivM-1],MPC_RNDNN);
	mpc_set(Mc[ndivM],Mc[ndivM-1],MPC_RNDNN);
	mpc_set(Md[ndivM],Md[ndivM-1],MPC_RNDNN);
	

	mpfr_clear(mu);
	mpfr_clear(k1);
	mpfr_clear(k2);
	mpfr_clear(k3);
	mpfr_clear(k4);
	mpfr_clear(kt);
	mpfr_clear(dmu);
	mpc_clear(temp);
	
	return;
}




std::string WaveFunction::getFastM(std::string mu0){
	
	std::string str = "empty";
	
	mpfr_t temp;
	mpfr_t mu;
	mpc_t M;
	mpc_t ctemp;
	
	long int index = 0;
	
	mpfr_init2(temp, prec_bit);
	mpfr_init2(mu, prec_bit);
	mpc_init2(M, prec_bit);
	mpc_init2(ctemp, prec_bit);
	
	mpfr_set_str(mu, &mu0[0], 10, GMP_RNDN);
	if(mpfr_cmp_ui(mu,0) == -1)
		mpfr_set_ui(mu,0,GMP_RNDN);
	
	if(mpfr_cmp_ui(mu,point1) == 1)	{
		mpc_set_ui(M,0,GMP_RNDN);
	char* gg = mpc_get_str(10, out_lenN, M, MPC_RNDNN); str = gg; delete gg;
	}	else	{
		
	mpfr_div_ui(temp,mu,point1,GMP_RNDN);
	mpfr_mul_ui(temp,temp,ndivM,GMP_RNDN);
	index = mpfr_get_si(temp,GMP_RNDZ);



	mpc_mul_fr(M,Ma[index],mu,MPC_RNDNN);	
	mpc_mul_fr(M,M,mu,MPC_RNDNN);			
	mpc_mul_fr(M,M,mu,MPC_RNDNN);			
	mpc_mul_fr(ctemp,Mb[index],mu,MPC_RNDNN);
	mpc_mul_fr(ctemp,ctemp,mu,MPC_RNDNN);	
	mpc_add(M,M,ctemp,MPC_RNDNN);			
	mpc_mul_fr(ctemp,Mc[index],mu,MPC_RNDNN);
	mpc_add(M,M,ctemp,MPC_RNDNN);			
	mpc_add(M,M,Md[index],MPC_RNDNN);		

	char* gg = mpc_get_str(10, out_lenN, M, MPC_RNDNN); str = gg; delete gg;
	}
	
	mpfr_clear(temp);
	mpfr_clear(mu);
	mpc_clear(M);
	mpc_clear(ctemp);
	
	return str;
}





long int WaveFunction::getnsumN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z1){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z1;
	mpc_t Z2;
	mpc_t sum_N;
	mpc_t sum_der_N;
	mpc_t der_N;
	mpc_t N;
	mpc_t temp;
	mpc_t temp1;
	mpc_t power;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	int k;
	
	

	std::string temp_str_N = "empty";
	std::string temp_str_der_N = "empty";
	
	mpfr_init2(F, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z1, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(sum_N, prec_bitN);
	mpc_init2(sum_der_N, prec_bitN);
	mpc_init2(N, prec_bitN);
	mpc_init2(der_N, prec_bitN);
	mpc_init2(power, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpc_init2(temp1, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);
	
	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	mpfr_neg(F,F,GMP_RNDN);


	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1, 1, MPC_RNDNN);
	

	

	mpc_set_ui(sum_N,1,MPC_RNDNN);	
	mpc_set_ui(sum_der_N,2*m + 1,MPC_RNDNN);
	mpc_div_ui(sum_der_N,sum_der_N,2,MPC_RNDNN);

	mpc_set_ui(sum_N,1,MPC_RNDNN);	
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
		else	{ 
			k=i;
			break;
		}
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

	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(Z2);
	mpc_clear(N);
	mpc_clear(der_N);
	mpc_clear(sum_N);
	mpc_clear(sum_der_N);
	mpc_clear(power);
	mpc_clear(temp);
	mpc_clear(temp1);
	mpc_clear(Ci_3);
	mpc_clear(Ci_2);
	mpc_clear(Ci_1);
	mpc_clear(Ci);
	
	
	return k+10;

}





void WaveFunction::getArrayN(std::string c_F, std::string c_energy, std::string c_Z1){
	
	mpfr_t F;
	mpc_t E;
	mpc_t Z1;
	mpc_t Z2;
	mpc_t temp;
	mpc_t Ci_3;
	mpc_t Ci_2;
	mpc_t Ci_1;
	mpc_t Ci;
	
	
	
	mpfr_init2(F, prec_bitN);
	mpc_init2(E, prec_bitN);
	mpc_init2(Z1, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);
	
	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	mpfr_neg(F,F,GMP_RNDN);


	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set_ui(Ci_1, 1, MPC_RNDNN);

	mpc_set_ui(c[0],1,MPC_RNDNN);
	

	
	for (long int i=1;i<nsumN+1;i++)	{
		
		mpc_mul_fr(Ci,Ci_3,F,MPC_RNDNN); 
		mpc_mul(temp,Ci_1,Z2,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_mul(temp,Ci_2,E,MPC_RNDNN); 
		mpc_mul_ui(temp,temp,2,MPC_RNDNN);
		mpc_sub(Ci,Ci,temp,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,4,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i,MPC_RNDNN);
		mpc_div_ui(Ci,Ci,i + m,MPC_RNDNN);

		mpc_set(c[i],Ci,MPC_RNDNN);

			
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);
		
		

		
	}


	mpfr_clear(F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(Z2);
	mpc_clear(temp);
	mpc_clear(Ci_3);
	mpc_clear(Ci_2);
	mpc_clear(Ci_1);
	mpc_clear(Ci);
	
	
	return;

}



void WaveFunction::getN(mpfr_t mu, mpc_t N){
	


	mpfr_t temp1;
	mpfr_t mu2;
	mpfr_t mu2i;
	mpc_t sum_N;

	mpc_t temp;

	std::string str;
	std::string temp_str_N = "empty";
	std::string str_N = "empty";
	


	mpfr_init2(temp1, prec_bitN);
	mpfr_init2(mu2, prec_bitN);
	mpfr_init2(mu2i, prec_bitN);
	mpc_init2(sum_N, prec_bitN);
	mpc_init2(temp, prec_bitN);



	mpfr_mul(mu2,mu,mu,GMP_RNDN);
	mpfr_set_ui(mu2i,1,GMP_RNDN);



	mpc_set_ui(sum_N,1,MPC_RNDNN);	
	for (long int i=1;i<nsumN+1;i++)	{
		mpfr_mul(mu2i,mu2i,mu2,GMP_RNDN);
		mpc_mul_fr(temp,c[i],mu2i,MPC_RNDNN);
		
		if (i%100==0)	{				
		char* gg = mpc_get_str(10, out_lenN, sum_N, MPC_RNDNN); str_N = gg; delete gg;
		if(str_N != temp_str_N)
		{
			temp_str_N = str_N;

		}
		else
			
			break;
	
		}
		

		mpc_add(sum_N,sum_N,temp,MPC_RNDNN);
	}


	mpfr_pow_ui(temp1, mu,m, GMP_RNDN);
	mpc_mul_fr(N,sum_N,temp1,GMP_RNDN);


	mpc_clear(sum_N);
	mpc_clear(temp);
	mpfr_clear(temp1);
	mpfr_clear(mu2);
	mpfr_clear(mu2i);


	
	return;

}


void WaveFunction::fillArrayN(){
	
	mpfr_t mu;
	mpfr_t dmu;
	mpfr_t k1;
	mpfr_t k2;
	mpfr_t k3;
	mpfr_t k4;
	mpfr_t kt;
	mpc_t temp;
	
	mpfr_init2(mu, prec_bitN);
	mpfr_init2(dmu, prec_bitN);
	mpfr_init2(k1, prec_bitN);
	mpfr_init2(k2, prec_bitN);
	mpfr_init2(k3, prec_bitN);
	mpfr_init2(k4, prec_bitN);
	mpfr_init2(kt, prec_bitN);
	mpc_init2(temp, prec_bitN);
	
	
	mpfr_set_ui(dmu,pointN1,GMP_RNDN);
	mpfr_div_ui(dmu,dmu,ndivN,GMP_RNDN);
	
	for (long int i=0;i<ndivN+1;i++)	{	
		mpfr_mul_ui(mu,dmu,i,GMP_RNDN);
		getN(mu,Ni[i]);
//		cout<<"N  "<<i<<"\n";
	}
	



	for (long int i=1;i<ndivN-1;i++)	{


		mpc_sub(Na[i],Ni[i+2],Ni[i-1],MPC_RNDNN);
		mpc_sub(temp,Ni[i],Ni[i+1],MPC_RNDNN);
		mpc_mul_ui(temp,temp,3,MPC_RNDNN);
		mpc_add(Na[i],Na[i],temp,MPC_RNDNN);
		mpc_div_fr(Na[i],Na[i],dmu,MPC_RNDNN);
		mpc_div_fr(Na[i],Na[i],dmu,MPC_RNDNN);
		mpc_div_fr(Na[i],Na[i],dmu,MPC_RNDNN);
		mpc_div_ui(Na[i],Na[i],6,MPC_RNDNN);

		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_add_ui(k1,k1,2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_add_ui(k2,k2,1,GMP_RNDN);

		
		mpc_mul_ui(Nb[i],Ni[i-1],1+i,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i],k1,MPC_RNDNN);
		mpc_sub(Nb[i],Nb[i],temp,MPC_RNDNN);	
		mpc_mul_fr(temp,Ni[i+1],k2,MPC_RNDNN);
		mpc_add(Nb[i],Nb[i],temp,MPC_RNDNN);	
		mpc_mul_ui(temp,Ni[i+2],i,MPC_RNDNN);
		mpc_sub(Nb[i],Nb[i],temp,MPC_RNDNN);
		mpc_div_fr(Nb[i],Nb[i],dmu,MPC_RNDNN);
		mpc_div_fr(Nb[i],Nb[i],dmu,MPC_RNDNN);
		mpc_div_ui(Nb[i],Nb[i],2,MPC_RNDNN);
		
		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_sub_ui(k1,k1,1,GMP_RNDN);
		
		mpfr_set_ui(k2,9,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k2,k2,kt,GMP_RNDN);
		mpfr_sub_ui(k2,k2,6,GMP_RNDN);
		
		mpfr_set_ui(k3,9,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_set_ui(kt,12,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k3,k3,kt,GMP_RNDN);
		mpfr_sub_ui(k3,k3,3,GMP_RNDN);
		
		mpfr_set_ui(k4,3,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k4,k4,kt,GMP_RNDN);
		mpfr_add_ui(k4,k4,2,GMP_RNDN);

		
		mpc_mul_fr(Nc[i],Ni[i+2],k1,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i+1],k2,MPC_RNDNN);
		mpc_sub(Nc[i],Nc[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i],k3,MPC_RNDNN);
		mpc_add(Nc[i],Nc[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i-1],k4,MPC_RNDNN);
		mpc_sub(Nc[i],Nc[i],temp,MPC_RNDNN);
		mpc_div_fr(Nc[i],Nc[i],dmu,MPC_RNDNN);
		mpc_div_ui(Nc[i],Nc[i],6,MPC_RNDNN);
		
		
		mpfr_set_ui(k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+1,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+2,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+1,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i-1,GMP_RNDN);
		
		mpfr_set_ui(k3,3,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i+2,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i-1,GMP_RNDN);
		
		mpfr_set_ui(k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i+1,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i-1,GMP_RNDN);

		
		
		mpc_mul_fr(Nd[i],Ni[i-1],k1,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i],k2,MPC_RNDNN);
		mpc_sub(Nd[i],Nd[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i+1],k3,MPC_RNDNN);
		mpc_add(Nd[i],Nd[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,Ni[i+2],k4,MPC_RNDNN);
		mpc_sub(Nd[i],Nd[i],temp,MPC_RNDNN);
		mpc_div_ui(Nd[i],Nd[i],6,MPC_RNDNN);
	
	}


	
	mpc_set(Na[0],Na[1],MPC_RNDNN);
	mpc_set(Nb[0],Nb[1],MPC_RNDNN);
	mpc_set(Nc[0],Nc[1],MPC_RNDNN);
	mpc_set(Nd[0],Nd[1],MPC_RNDNN);
	
	mpc_set(Na[ndivN-1],Na[ndivN-2],MPC_RNDNN);
	mpc_set(Nb[ndivN-1],Nb[ndivN-2],MPC_RNDNN);
	mpc_set(Nc[ndivN-1],Nc[ndivN-2],MPC_RNDNN);
	mpc_set(Nd[ndivN-1],Nd[ndivN-2],MPC_RNDNN);

	mpc_set(Na[ndivN],Na[ndivN-1],MPC_RNDNN);
	mpc_set(Nb[ndivN],Nb[ndivN-1],MPC_RNDNN);
	mpc_set(Nc[ndivN],Nc[ndivN-1],MPC_RNDNN);
	mpc_set(Nd[ndivN],Nd[ndivN-1],MPC_RNDNN);
	

	mpfr_clear(mu);
	mpfr_clear(k1);
	mpfr_clear(k2);
	mpfr_clear(k3);
	mpfr_clear(k4);
	mpfr_clear(kt);
	mpfr_clear(dmu);
	mpc_clear(temp);
	
	return;
}




std::string WaveFunction::getFastNbelow(std::string mu0){
	
	std::string str = "empty";
	
	mpfr_t temp;
	mpfr_t mu;
	mpc_t N;
	mpc_t ctemp;
	
	long int index = 0;
	
	mpfr_init2(temp, prec_bitN);
	mpfr_init2(mu, prec_bitN);
	mpc_init2(N, prec_bitN);
	mpc_init2(ctemp, prec_bitN);
	mpfr_set_str(mu, &mu0[0], 10, GMP_RNDN);
	mpfr_div_ui(temp,mu,pointN1,GMP_RNDN);
	mpfr_mul_ui(temp,temp,ndivN,GMP_RNDN);
	index = mpfr_get_si(temp,GMP_RNDZ);
	if (index<0)	index=0;

	mpc_mul_fr(N,Na[index],mu,MPC_RNDNN);
	mpc_mul_fr(N,N,mu,MPC_RNDNN);
	mpc_mul_fr(N,N,mu,MPC_RNDNN);	
	mpc_mul_fr(ctemp,Nb[index],mu,MPC_RNDNN);
	mpc_mul_fr(ctemp,ctemp,mu,MPC_RNDNN);
	mpc_add(N,N,ctemp,MPC_RNDNN);	
	mpc_mul_fr(ctemp,Nc[index],mu,MPC_RNDNN);
	mpc_add(N,N,ctemp,MPC_RNDNN);
	mpc_add(N,N,Nd[index],MPC_RNDNN);


	
	
	char* gg = mpc_get_str(10, out_lenN, N, MPC_RNDNN); str = gg; delete gg;
	
	mpfr_clear(temp);
	mpfr_clear(mu);
	mpc_clear(N);
	mpc_clear(ctemp);
	
	return str;
}







long int WaveFunction::calcPrecisionForN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	long int crit;
	


	
	while(true)	{
		crit =  out_lenN + precisionN - getN(str_N,  str_derN, c_F,  c_energy,  c_Z2);
		if (crit>=500)	{pointN1++; pointN2 = pointN1*pointN1;}
		if (crit<=300)	{pointN1--; pointN2 = pointN1*pointN1;}
		if ((crit>300)&&(crit<500))	break;
	}

	return (int)(out_lenN*3.3219280948873626);
}



long int WaveFunction::getN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z2)	{


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
		mpfr_log10(abs_sum_N,absCi,GMP_RNDN);
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





long int WaveFunction::get_a(std::string& str_a, std::string& str_der_a,std::string c_F, std::string c_energy, std::string c_Z1){
	
	
	std::string temp_str_a = "empty";
	std::string temp_str_der_a = "empty";
	long int i;
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z1;
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
	mpc_init2(Z1, prec_bitN);
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
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);


	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	
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
	for (i=2;i<100000;i++)	{
		
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
		
//		std::cout<<"C"<<i<<" = "<<mpc_get_str(10, out_lenN, sum_a, MPC_RNDNN)<<"\n";		
		
		
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
	mpc_clear(Z1);
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




long int WaveFunction::get_b(std::string& str_b, std::string& str_der_b,std::string c_F, std::string c_energy, std::string c_Z1){
	
	
	std::string temp_str_b = "empty";
	std::string temp_str_der_b = "empty";
	long int i;
	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z1;
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
	mpc_init2(Z1, prec_bitN);
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
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);

	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	
	
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
	for (i=2;i<100000;i++)	{
		
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
		
//		std::cout<<"C"<<i<<" = "<<mpc_get_str(10, out_lenN, sum_b, MPC_RNDNN)<<"\n";
	
		
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
	mpc_clear(Z1);
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








void WaveFunction::getAB(std::string str_a,std::string str_der_a,std::string str_b,std::string str_der_b,std::string str_N,std::string str_der_N,std::string& str_A,std::string& str_B){
	


	long int na, nb;
		
		mpc_t A;
		mpc_t B;
		mpc_t a;
		mpc_t der_a;
		mpc_t b;
		mpc_t der_b;
		mpc_t N;
		mpc_t der_N;

		mpc_t temp1;
		mpc_t temp2;
	

		mpc_init2(A, prec_bitN);
		mpc_init2(B, prec_bitN);
		mpc_init2(a, prec_bitN);
		mpc_init2(der_a, prec_bitN);
		mpc_init2(b, prec_bitN);
		mpc_init2(der_b, prec_bitN);
		mpc_init2(N, prec_bitN);
		mpc_init2(der_N, prec_bitN);
		mpc_init2(temp1, prec_bitN);
		mpc_init2(temp2, prec_bitN);


	
	
	
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
	
	mpc_mul(temp1,N,der_b,MPC_RNDNN);
	mpc_mul(temp2,b,der_N,MPC_RNDNN);
	mpc_sub(A,temp1,temp2,MPC_RNDNN);
	mpc_mul(temp1,a,der_b,MPC_RNDNN);
	mpc_mul(temp2,b,der_a,MPC_RNDNN);
	mpc_sub(temp1,temp1,temp2,MPC_RNDNN);
	mpc_div(A,A,temp1,MPC_RNDNN);
	


	char* gg = mpc_get_str(10, out_lenN, B, MPC_RNDNN); str_B = gg; delete gg;
	char* hh = mpc_get_str(10, out_lenN, A, MPC_RNDNN); str_A = hh; delete hh;

	
	mpc_clear(A);
	mpc_clear(B);
	mpc_clear(a);
	mpc_clear(der_a);
	mpc_clear(b);
	mpc_clear(der_b);
	mpc_clear(N);
	mpc_clear(der_N);

	mpc_clear(temp1);
	mpc_clear(temp2);
	

	
	return;
}




void WaveFunction::getArray_a(std::string str_A,std::string c_F,std::string c_energy, std::string c_Z1){
	
	
	std::string temp_str_a = "empty";
	std::string temp_str_der_a = "empty";

	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z1;
	mpc_t Z2;
	mpc_t A;


	mpc_t temp;
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
	mpc_init2(Z1, prec_bitN);

	mpc_init2(A, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);
	mpc_set_str(A, &str_A[0], 10, MPC_RNDNN);


	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	mpc_set_fr(temp,F,MPC_RNDNN);

	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);	
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);
	mpc_neg(k2,k2,MPC_RNDNN);
	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	

	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set(Ci_1,A, MPC_RNDNN);

	mpc_set(a[1],A,MPC_RNDNN);

	for (long int i=2;i<nsumab;i++)	{
		
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
				

		mpc_set(a[i],Ci,MPC_RNDNN);

		
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

	}
	

	


	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z2);
	mpc_clear(Z1);

	mpc_clear(A);
	mpc_clear(temp);
	mpc_clear(Ci);
	mpc_clear(Ci_1);
	mpc_clear(Ci_2);
	mpc_clear(Ci_3);

	return;
	
}





void WaveFunction::getArray_b(std::string str_B,std::string c_F,std::string c_energy, std::string c_Z1){
	
	
	std::string temp_str_b = "empty";
	std::string temp_str_der_b = "empty";

	
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;
	mpc_t Z1;
	mpc_t Z2;
	mpc_t B;

	mpc_t temp;

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
	mpc_init2(Z1, prec_bitN);
	mpc_init2(Z2, prec_bitN);
	mpc_init2(B, prec_bitN);

	mpc_init2(temp, prec_bitN);
	mpc_init2(Ci_3, prec_bitN);
	mpc_init2(Ci_2, prec_bitN);
	mpc_init2(Ci_1, prec_bitN);
	mpc_init2(Ci, prec_bitN);
	
	

	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	mpc_set_str(Z1, &c_Z1[0], 10, MPC_RNDNN);
	mpc_set_str(B, &str_B[0], 10, MPC_RNDNN);

	mpc_ui_sub(Z2,4,Z1,MPC_RNDNN);
	
	
	//calculation of sqrt(F)
	mpfr_sqrt(sqrt_F,F,GMP_RNDN);
	mpc_set_fr(temp,F,MPC_RNDNN);

	
	//calculation of constant k1
	mpc_mul(k1,E,E,MPC_RNDNN);
	mpc_div_fr(k1,k1,F,MPC_RNDNN);
	mpc_sub(k1,k1,Z2,MPC_RNDNN);
	mpc_mul_ui(k1,k1,4,MPC_RNDNN);
	
	//calculation of constant k2
	mpc_set_ui_ui(k2,0,8,MPC_RNDNN);	
	mpc_mul_fr(k2,k2,sqrt_F,MPC_RNDNN);
	
//	mpc_neg(k2,k2,MPC_RNDNN);

	
	//calculation of constant k3
	mpc_set_ui_ui(k3,0,8,MPC_RNDNN);
	mpc_mul(k3,k3,E,MPC_RNDNN);
	mpc_div_fr(k3,k3,sqrt_F,MPC_RNDNN);
	

	

	mpc_set_ui(Ci_3, 0, MPC_RNDNN);
	mpc_set_ui(Ci_2, 0, MPC_RNDNN);
	mpc_set(Ci_1,B, MPC_RNDNN);


	mpc_set(b[1],B, MPC_RNDNN);


	for (long int i=2;i<nsumab;i++)	{
		
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
		

		mpc_set(b[i],Ci, MPC_RNDNN);
		
		
		mpc_set(Ci_3,Ci_2,MPC_RNDNN);
		mpc_set(Ci_2,Ci_1,MPC_RNDNN);
		mpc_set(Ci_1,Ci,MPC_RNDNN);

	}
	

	
	mpc_clear(k1);
	mpc_clear(k2);
	mpc_clear(k3);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(Z1);
	mpc_clear(Z2);
	mpc_clear(B);	
	mpc_clear(temp);
	mpc_clear(Ci);
	mpc_clear(Ci_1);
	mpc_clear(Ci_2);
	mpc_clear(Ci_3);
	
	return;
	
}









void WaveFunction::get_ab(mpfr_t nu,std::string c_F,std::string c_energy, mpc_t ab){
	
	
	
	mpfr_t temp1;
	mpfr_t nui;
	mpfr_t sqrt_nu;
	mpc_t sum_a;
	mpc_t sum_b;
	mpc_t temp;
	mpfr_t F;
	mpfr_t sqrt_F;
	mpc_t E;

	std::string str;
	std::string temp_str_a = "empty";
	std::string str_a = "empty";
	std::string temp_str_b = "empty";
	std::string str_b = "empty";


	mpfr_init2(temp1, prec_bitN);
	mpfr_init2(nui, prec_bitN);
	mpfr_init2(sqrt_nu, prec_bitN);
	mpc_init2(sum_a, prec_bitN);
	mpc_init2(sum_b, prec_bitN);
	mpc_init2(temp, prec_bitN);
	mpfr_init2(F, prec_bitN);
	mpfr_init2(sqrt_F, prec_bitN);
	mpc_init2(E, prec_bitN);
	
	
	
	
	mpfr_set_str(F, &c_F[0], 10, GMP_RNDN);
	mpc_set_str(E, &c_energy[0], 10, MPC_RNDNN);
	

	mpfr_sqrt(sqrt_F,F,GMP_RNDN);

	
	
	
	
	mpfr_set(nui,nu,GMP_RNDN);
	mpc_div_fr(sum_a,a[1],nu,MPC_RNDNN);	

	for (long int i=2;i<nsumab;i++)	{
		mpfr_mul(nui,nui,nu,GMP_RNDN);
		mpc_div_fr(temp,a[i],nui,MPC_RNDNN);
		
		if (i%10==0)	{				
			char* gg = mpc_get_str(10, out_lenN, sum_a, MPC_RNDNN); str_a = gg; delete gg;
			if(str_a != temp_str_a)
			{
				temp_str_a = str_a;

			}
			else
				
				break;
		
			}
			

			mpc_add(sum_a,sum_a,temp,MPC_RNDNN);

	}
	

	
	mpfr_mul(temp1,sqrt_F,nu,GMP_RNDN);
	mpfr_mul(temp1,temp1,nu,GMP_RNDN);
	mpfr_mul(temp1,temp1,nu,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_fr(temp,E,nu,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,1,MPC_RNDNN);
	mpc_exp(temp,temp,MPC_RNDNN);
	

	mpc_mul(sum_a,sum_a,temp,MPC_RNDNN);
	
	
	
	
	mpfr_set(nui,nu,GMP_RNDN);
	mpc_div_fr(sum_b,b[1],nu,MPC_RNDNN);	

	for (long int i=2;i<nsumab;i++)	{
		mpfr_mul(nui,nui,nu,GMP_RNDN);
		mpc_div_fr(temp,b[i],nui,MPC_RNDNN);
		
		if (i%10==0)	{				
			char* gg = mpc_get_str(10, out_lenN, sum_b, MPC_RNDNN); str_b = gg; delete gg;
			if(str_b != temp_str_b)
			{
				temp_str_b = str_b;

			}
			else
				
				break;
		
			}
			

			mpc_add(sum_b,sum_b,temp,MPC_RNDNN);

	}
	

	
	mpfr_mul(temp1,sqrt_F,nu,GMP_RNDN);
	mpfr_mul(temp1,temp1,nu,GMP_RNDN);
	mpfr_mul(temp1,temp1,nu,GMP_RNDN);
	mpfr_div_ui(temp1,temp1,3,GMP_RNDN);
	mpc_mul_fr(temp,E,nu,MPC_RNDNN);
	mpc_div_fr(temp,temp,sqrt_F,MPC_RNDNN);
	mpc_add_fr(temp, temp,temp1,MPC_RNDNN);
	mpc_mul_i(temp,temp,-1,MPC_RNDNN);
	mpc_exp(temp,temp,MPC_RNDNN);
	

	mpc_mul(sum_b,sum_b,temp,MPC_RNDNN);
	
	
	

	mpc_add(ab,sum_a,sum_b,MPC_RNDNN);
	mpfr_sqrt(sqrt_nu,nu,GMP_RNDN);
	mpc_div_fr(ab,ab,sqrt_nu,MPC_RNDNN);
	
	
	mpc_clear(temp);
	mpfr_clear(F);
	mpfr_clear(sqrt_F);
	mpc_clear(E);
	mpc_clear(sum_b);
	mpc_clear(sum_a);
	mpfr_clear(nui);
	mpfr_clear(sqrt_nu);
	mpfr_clear(temp1);
	
	return;
	
}


void WaveFunction::fillArray_ab(std::string c_F,std::string c_energy){
	
	mpfr_t nu;
	mpfr_t dnu;
	mpfr_t k1;
	mpfr_t k2;
	mpfr_t k3;
	mpfr_t k4;
	mpfr_t kt;
	mpc_t temp;
	
	mpfr_init2(nu, prec_bit);
	mpfr_init2(dnu, prec_bit);
	mpfr_init2(k1, prec_bit);
	mpfr_init2(k2, prec_bit);
	mpfr_init2(k3, prec_bit);
	mpfr_init2(k4, prec_bit);
	mpfr_init2(kt, prec_bit);
	mpc_init2(temp, prec_bit);
	
	
	mpfr_sub_ui(dnu,h,pointN1,GMP_RNDN);
	mpfr_div_ui(dnu,dnu,ndivab,GMP_RNDN);
	
	for (long int i=0;i<ndivab+1;i++)	{	
		mpfr_mul_ui(nu,dnu,i,GMP_RNDN);
		mpfr_add_ui(nu,nu,pointN1,GMP_RNDN);
		get_ab(nu,c_F,c_energy,abi[i]);
//		cout<<"ab "<<i<<"\n";

	}
	


	
	for (long int i=1;i<ndivab-1;i++)	{
		
		mpc_sub(ab_a[i],abi[i+2],abi[i-1],MPC_RNDNN);
		mpc_sub(temp,abi[i],abi[i+1],MPC_RNDNN);
		mpc_mul_ui(temp,temp,3,MPC_RNDNN);
		mpc_add(ab_a[i],ab_a[i],temp,MPC_RNDNN);
		mpc_div_fr(ab_a[i],ab_a[i],dnu,MPC_RNDNN);
		mpc_div_fr(ab_a[i],ab_a[i],dnu,MPC_RNDNN);
		mpc_div_fr(ab_a[i],ab_a[i],dnu,MPC_RNDNN);
		mpc_div_ui(ab_a[i],ab_a[i],6,MPC_RNDNN);
	
		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_add_ui(k1,k1,2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_add_ui(k2,k2,1,GMP_RNDN);

		
		mpc_mul_ui(ab_b[i],abi[i-1],1+i,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i],k1,MPC_RNDNN);
		mpc_sub(ab_b[i],ab_b[i],temp,MPC_RNDNN);	
		mpc_mul_fr(temp,abi[i+1],k2,MPC_RNDNN);
		mpc_add(ab_b[i],ab_b[i],temp,MPC_RNDNN);	
		mpc_mul_ui(temp,abi[i+2],i,MPC_RNDNN);
		mpc_sub(ab_b[i],ab_b[i],temp,MPC_RNDNN);
		mpc_div_fr(ab_b[i],ab_b[i],dnu,MPC_RNDNN);
		mpc_div_fr(ab_b[i],ab_b[i],dnu,MPC_RNDNN);
		mpc_div_ui(ab_b[i],ab_b[i],2,MPC_RNDNN);
		
		mpfr_set_ui(k1,3,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i,GMP_RNDN);
		mpfr_sub_ui(k1,k1,1,GMP_RNDN);
		
		mpfr_set_ui(k2,9,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k2,k2,kt,GMP_RNDN);
		mpfr_sub_ui(k2,k2,6,GMP_RNDN);
		
		mpfr_set_ui(k3,9,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_set_ui(kt,12,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k3,k3,kt,GMP_RNDN);
		mpfr_sub_ui(k3,k3,3,GMP_RNDN);
		
		mpfr_set_ui(k4,3,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i,GMP_RNDN);
		mpfr_set_ui(kt,6,GMP_RNDN);
		mpfr_mul_ui(kt,kt,i,GMP_RNDN);
		mpfr_add(k4,k4,kt,GMP_RNDN);
		mpfr_add_ui(k4,k4,2,GMP_RNDN);

		
		mpc_mul_fr(ab_c[i],abi[i+2],k1,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i+1],k2,MPC_RNDNN);
		mpc_sub(ab_c[i],ab_c[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i],k3,MPC_RNDNN);
		mpc_add(ab_c[i],ab_c[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i-1],k4,MPC_RNDNN);
		mpc_sub(ab_c[i],ab_c[i],temp,MPC_RNDNN);
		mpc_div_fr(ab_c[i],ab_c[i],dnu,MPC_RNDNN);
		mpc_div_ui(ab_c[i],ab_c[i],6,MPC_RNDNN);
		
		
		mpfr_set_ui(k1,i,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+1,GMP_RNDN);
		mpfr_mul_ui(k1,k1,i+2,GMP_RNDN);
		
		mpfr_set_ui(k2,3,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+2,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i+1,GMP_RNDN);
		mpfr_mul_ui(k2,k2,i-1,GMP_RNDN);
		
		mpfr_set_ui(k3,3,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i+2,GMP_RNDN);
		mpfr_mul_ui(k3,k3,i-1,GMP_RNDN);
		
		mpfr_set_ui(k4,i,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i+1,GMP_RNDN);
		mpfr_mul_ui(k4,k4,i-1,GMP_RNDN);

		
		
		mpc_mul_fr(ab_d[i],abi[i-1],k1,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i],k2,MPC_RNDNN);
		mpc_sub(ab_d[i],ab_d[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i+1],k3,MPC_RNDNN);
		mpc_add(ab_d[i],ab_d[i],temp,MPC_RNDNN);
		mpc_mul_fr(temp,abi[i+2],k4,MPC_RNDNN);
		mpc_sub(ab_d[i],ab_d[i],temp,MPC_RNDNN);
		mpc_div_ui(ab_d[i],ab_d[i],6,MPC_RNDNN);
	
	}


	
	mpc_set(ab_a[0],ab_a[1],MPC_RNDNN);
	mpc_set(ab_b[0],ab_b[1],MPC_RNDNN);
	mpc_set(ab_c[0],ab_c[1],MPC_RNDNN);
	mpc_set(ab_d[0],ab_d[1],MPC_RNDNN);
	
	mpc_set(ab_a[ndivab-1],ab_a[ndivab-2],MPC_RNDNN);
	mpc_set(ab_b[ndivab-1],ab_b[ndivab-2],MPC_RNDNN);
	mpc_set(ab_c[ndivab-1],ab_c[ndivab-2],MPC_RNDNN);
	mpc_set(ab_d[ndivab-1],ab_d[ndivab-2],MPC_RNDNN);

	mpc_set(ab_a[ndivab],ab_a[ndivab-1],MPC_RNDNN);
	mpc_set(ab_b[ndivab],ab_b[ndivab-1],MPC_RNDNN);
	mpc_set(ab_c[ndivab],ab_c[ndivab-1],MPC_RNDNN);
	mpc_set(ab_d[ndivab],ab_d[ndivab-1],MPC_RNDNN);
	

	mpfr_clear(nu);
	mpfr_clear(k1);
	mpfr_clear(k2);
	mpfr_clear(k3);
	mpfr_clear(k4);
	mpfr_clear(kt);
	mpfr_clear(dnu);
	mpc_clear(temp);
	
	return;
}



std::string WaveFunction::getFast_ab(std::string nu0){
	
	std::string str = "empty";
	
	mpfr_t temp;
	mpfr_t temp1;
	mpfr_t nu;
	mpc_t ab;
	mpc_t ctemp;
	
	long int index = 0;
	
	mpfr_init2(temp, prec_bitN);
	mpfr_init2(temp1, prec_bitN);
	mpfr_init2(nu, prec_bitN);
	mpc_init2(ab, prec_bitN);
	mpc_init2(ctemp, prec_bitN);
	
	mpfr_set_str(nu, &nu0[0], 10, GMP_RNDN);
	mpfr_sub_ui(temp,h,pointN1,GMP_RNDN);
	mpfr_sub_ui(temp1,nu,pointN1,GMP_RNDN);
	mpfr_div(temp,temp1,temp,GMP_RNDN);
	mpfr_mul_ui(temp,temp,ndivab,GMP_RNDN);
	index = mpfr_get_si(temp,GMP_RNDZ);
	if (index<0)	index=0;



	mpfr_sub_ui(nu,nu,pointN1,GMP_RNDN);
	mpc_mul_fr(ab,ab_a[index],nu,MPC_RNDNN);
	mpc_mul_fr(ab,ab,nu,MPC_RNDNN);
	mpc_mul_fr(ab,ab,nu,MPC_RNDNN);	
	mpc_mul_fr(ctemp,ab_b[index],nu,MPC_RNDNN);
	mpc_mul_fr(ctemp,ctemp,nu,MPC_RNDNN);
	mpc_add(ab,ab,ctemp,MPC_RNDNN);	
	mpc_mul_fr(ctemp,ab_c[index],nu,MPC_RNDNN);
	mpc_add(ab,ab,ctemp,MPC_RNDNN);
	mpc_add(ab,ab,ab_d[index],MPC_RNDNN);



	
	
	char* gg = mpc_get_str(10, out_lenN, ab, MPC_RNDNN); str = gg; delete gg;
	
	mpfr_clear(temp);
	mpfr_clear(temp1);
	mpfr_clear(nu);
	mpc_clear(ab);
	mpc_clear(ctemp);
	
	return str;
}


std::string WaveFunction::getFastN(std::string nu0){
	
	mpfr_t nu;
	mpfr_init2(nu, prec_bitN);
	mpfr_set_str(nu, &nu0[0], 10, GMP_RNDN);
	
	if(mpfr_cmp_ui(nu,0) == -1)	{
		mpfr_set_ui(nu,0,GMP_RNDN);
		nu0 = "0";
	}
		


	
if(mode==1||mode==2){
	
	if (mode == 1)	{
		mpfr_clear(nu);
		return getFastNbelow(nu0);
	}
	if (mode == 2)
	{
		if(mpfr_cmp_ui(nu,pointN1) != 1)	{
			mpfr_clear(nu);
			return getFastNbelow(nu0);
		}
		
		if(mpfr_cmp_ui(nu,pointN1) == 1)	{
			mpfr_clear(nu);
			return getFast_ab(nu0);
		}
	}
}
else 
	{
	mpfr_clear(nu);
	return "(0 0)";
}

	
}


int WaveFunction::getMode()	{
	
	return mode;
}














///////////////////////////////////////////////////////////////////Below is the functions wave functions for perturbation theory of the second order

void WaveFunction::getArrayb()	{

	bp = new long int**[15];	for (int i=0;i<15;i++)	bp[i]=new long int*[15]; for (int i=0;i<15;i++) for (int j=0;j<15;j++)	bp[i][j]=new long int[15];
	
	
	for (long int i=0;i<15;i++)
		for (long int j=0;j<15;j++)
			for (long int k=0;k<15;k++)
				bp[k][i][j] = 0;
    
    int dn = n1 - n2;
    
    bp[0][0][0] = 1;

    bp[1][-2  +2][0  +2] = n1*(-1+n1);
    bp[1][-1  +2][0  +2] = -2*(-2-dn+2*n)*n1;
    bp[1][0   +2][0  +2] = 6*dn;
    bp[1][1   +2][0  +2] = 2*(2-dn+2*n)*(1+m+n1);
    bp[1][2   +2][0  +2] = -(1+m+n1)*(2+m+n1);
    bp[1][0   +2][-2 +2] = -(-1+n2)*n2;
    bp[1][0   +2][-1 +2] = 2*(-2+dn+2*n)*n2;
    bp[1][0   +2][1  +2] = -2*(2+dn+2*n)*(1+m+n2);
    bp[1][0   +2][2  +2] = (1+m+n2)*(2+m+n2);
    
    bp[2][-2  +4][-2  +4] = -(-1+n1)*n1*(-1+n2)*n2;
    bp[2][-1  +4][-1  +4] = -4*(-2-dn+2*n)*(-2+dn+2*n)*n1*n2;
    bp[2][0  +4][0  +4] = (-1289+74*dn*dn-17*dn*dn*dn*dn+650*m*m+18*dn*dn*m*m-m*m*m*m-1494*n*n+42*dn*dn*n*n+66*m*m*n*n-65*n*n*n*n)/8;
    bp[2][1  +4][1  +4] = -4*(2-dn+2*n)*(2+dn+2*n)*(1+m+n1)*(1+m+n2);
    bp[2][2  +4][2  +4] = -(1+m+n1)*(2+m+n1)*(1+m+n2)*(2+m+n2);  
    bp[2][-4  +4][0  +4] = (-3+n1)*(-2+n1)*(-1+n1)*n1/2;
    bp[2][-3  +4][0  +4] = -2*(-14-3*dn+6*n)*(-2+n1)*(-1+n1)*n1/3;
    bp[2][-2  +4][-1  +4] = 2*(-2+dn+2*n)*(-1+n1)*n1*n2;
    bp[2][-2  +4][0  +4] = 2*(-8-dn+2*n)*(-3-dn+2*n)*(-1+n1)*n1;
    bp[2][-2  +4][1  +4] = -2*(2+dn+2*n)*(-1+n1)*n1*(1+m+n2)   ;     
    bp[2][-2  +4][2  +4] = (-1+n1)*n1*(1+m+n2)*(2+m+n2);
    bp[2][-1  +4][0  +4] = (78+75*dn+36*dn*dn-dn*dn*dn+26*m*m+dn*m*m-234*n-66*dn*n-2*m*m*n+58*n*n+3*dn*n*n+2*n*n*n)*n1/2;
    bp[2][-1  +4][1  +4] = 4*(-2-dn+2*n)*(2+dn+2*n)*n1*(1+m+n2);
    bp[2][-1  +4][2  +4] = -2*(-2-dn+2*n)*n1*(1+m+n2)*(2+m+n2);
    bp[2][0  +4][1  +4] = (-78-75*dn-36*dn*dn+dn*dn*dn-26*m*m-dn*m*m-234*n-66*dn*n-2*m*m*n-58*n*n-3*dn*n*n+2*n*n*n)*(1+m+n2)/2   ;     
    bp[2][0  +4][2  +4] = 2*(3+dn+2*n)*(8+dn+2*n)*(1+m+n2)*(2+m+n2);
    bp[2][0  +4][3  +4] = -2*(14+3*dn+6*n)*(1+m+n2)*(2+m+n2)*(3+m+n2)/3;
    bp[2][0  +4][4  +4] = (1+m+n2)*(2+m+n2)*(3+m+n2)*(4+m+n2)/2;
    bp[2][1  +4][2  +4] = 2*(2-dn+2*n)*(1+m+n1)*(1+m+n2)*(2+m+n2);        
    bp[2][0  +4][-4  +4] = (-3+n2)*(-2+n2)*(-1+n2)*n2/2;
    bp[2][0  +4][-3  +4] = -2*(-14+3*dn+6*n)*(-2+n2)*(-1+n2)*n2/3;
    bp[2][-1  +4][-2  +4] = 2*(-2-dn+2*n)*(-1+n2)*n1*n2;
    bp[2][0  +4][-2  +4] = 2*(-8+dn+2*n)*(-3+dn+2*n)*(-1+n2)*n2;
    bp[2][1  +4][-2  +4] = -2*(2-dn+2*n)*(-1+n2)*n2*(1+m+n1);
    bp[2][2  +4][-2  +4] =  (-1+n2)*n2*(1+m+n1)*(2+m+n1);
    bp[2][0  +4][-1  +4] = (78-75*dn+36*dn*dn+dn*dn*dn+26*m*m-dn*m*m-234*n+66*dn*n-2*m*m*n+58*n*n-3*dn*n*n+2*n*n*n)*n2/2;
    bp[2][1  +4][-1  +4] = 4*(-2+dn+2*n)*(2-dn+2*n)*n2*(1+m+n1);
    bp[2][2  +4][-1  +4] = -2*(-2+dn+2*n)*n2*(1+m+n1)*(2+m+n1);
    bp[2][1  +4][0  +4] = (-78+75*dn-36*dn*dn-dn*dn*dn-26*m*m+dn*m*m-234*n+66*dn*n-2*m*m*n-58*n*n+3*dn*n*n+2*n*n*n)*(1+m+n1)/2;
    bp[2][2  +4][0  +4] = 2*(3-dn+2*n)*(8-dn+2*n)*(1+m+n1)*(2+m+n1);
    bp[2][3  +4][0  +4] = -2*(14-3*dn+6*n)*(1+m+n1)*(2+m+n1)*(3+m+n1)/3;
    bp[2][4  +4][0  +4] = (1+m+n1)*(2+m+n1)*(3+m+n1)*(4+m+n1)/2;
    bp[2][2  +4][1  +4] = 2*(2+dn+2*n)*(1+m+n2)*(1+m+n1)*(2+m+n1);


	
return;	
	
}


double WaveFunction::f(long int k, long int _m, double x){
	
	double val = 0;
	int co = 1;
	for (int p=0;p<k+1;p++)	{
		val += co*MathPolynomial::FactorialDouble(k)*pow(x,p+0.5*_m)/(MathPolynomial::FactorialDouble(k-p)*MathPolynomial::FactorialDouble(_m+p)*MathPolynomial::FactorialDouble(p));
		co = -co;
	}
		
	return val*exp(-x*0.5);
}


double WaveFunction::getMN(double mu,double nu){

	double sum1,sum2,sum3;
	double PI = 3.14159265358979323846264338;
	
		sum1 = 0;
		for ( int N=0;N<2+1;N++)	{
			sum2 = 0;
			for ( int i=-2*N;i<2*N+1;i++){
				sum3 = 0;
				for ( int j=-2*N;j<2*N+1;j++)	{
					
					sum3 += bp[N][i+2*N][j+2*N]*f(n2+j,m,nu*nu/n);
					
				}
				sum2 += sum3*f(n1+i,m,mu*mu/n);
			}
			sum1 += sum2*pow(Fd,N)*pow(n*0.5,3*N);
		}
		
	return A*sum1/sqrt(2*PI);
}




