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
#include "FuncSS.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>





using namespace OrbitUtils;
using namespace std;




FuncSS::FuncSS(int n11,int n22, int mm, int point11)	{		
		
	n1 = n11;
	n2 = n22;
	m = abs(mm);
	n = n1 + n2 + m + 1; 
	
	point1 = point11;
	point2 = point1*point1;
	
	out_len = 10;
	out_lenN = 100;
	
	pointN1 = 50;
	pointN2 = pointN1*pointN1;
	
	
	precisionN = 1000;
	nsumab = 100000;
	prec_bitN = (int)(precisionN*3.3219280948873626);

	
	nsum1 = 100;
	precision = 100;
	prec_bit = (int)(precision*3.3219280948873626);


}





FuncSS::~FuncSS()	{
	
}




long int FuncSS::calcPrecisionForM(int max_err_exp, std::string c_field, std::string c_energy, std::string c_Z1)	{


	std::string line = "empty line", Mstr;
	
	int sum = 0;
	int prec = 1;
	int _temp;
	int count = 0;
	int exponent;
				
		
	
	
	
	for (int j=0;j<2;j++) 	{

		
		
	precision = 10;
	nsum1 = 10;


		
	while (true)	{
		
		
		nsum1 += 10*sum;
		precision += 10*prec;
		
		prec_bit = (int)(precision*3.3219280948873626);
		

		
		if(j==0) out_len = 10;
		if(j==1) out_len = exponent - max_err_exp;
	
//			cout<<"out_len = "<<out_len<<"\n";


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
	mpfr_class temp("0",prec_bit);
	mpfr_class M(Mstr,prec_bit);
	
	temp = log10(abs(M));
	exponent = (int)(temp.get_d());
		if (exponent<5) exponent = 5;		

	}	


//	cout<<out_len<<" "<<nsum1<<"  "<<precision<<endl;
//	out_len = precision;
	
	return prec_bit + 10;
	


}










std::string FuncSS::getM(std::string c_F, std::string c_energy, std::string c_Z1){
	
	std::stringstream out;
	out.precision(out_len);
	
	mpfr_class F(c_F,prec_bit);
	mpfr_class E(c_energy,prec_bit);
	mpfr_class Z1(c_Z1,prec_bit);

    mpfr_class Ci_3("0",prec_bit);
    mpfr_class Ci_2("0",prec_bit);
    mpfr_class Ci_1("1",prec_bit);
    mpfr_class Ci("0",prec_bit);
	
	mpfr_class sum_M("1",prec_bit);

	for (int i=1;i<nsum1+1;i++)	{
				
		Ci = (Ci_3*F*point2*point2 - Ci_1*Z1 - Ci_2*E*2*point2)*point2/4/i/(i + m);
		sum_M += Ci;
		
		Ci_3 = Ci_2;
		Ci_2 = Ci_1;
		Ci_1 = Ci;
	}

	mpfr_class p1("0",prec_bit);
	p1 += point1;
	
    mpfr_class power("0.5",prec_bit);
	power += m;
	mpfr_class M = sum_M*pow(p1,power);

	
	out<<M;
	return out.str();

}


mpfr_class FuncSS::integrateM2(mpfr_class F,mpfr_class E, mpfr_class Z1)	{
	


	mpfr_class C[nsum1 + 1]; for(int i=0;i<nsum1 + 1;i++) C[i].set_prec(prec_bit);


    mpfr_class Ci_3("0",prec_bit);
    mpfr_class Ci_2("0",prec_bit);
    mpfr_class Ci_1("1",prec_bit);
    mpfr_class Ci("0",prec_bit);


	C[0] = "1";

	for (int i=1;i<nsum1+1;i++)	{

		Ci = (Ci_3*F*point2*point2 - Ci_1*Z1 - Ci_2*E*2*point2)*point2/4/i/(i + m);
		C[i] = Ci;

		Ci_3 = Ci_2;
		Ci_2 = Ci_1;
		Ci_1 = Ci;
	}


	
	mpfr_class integral("0",prec_bit);

	for (int n=0;n<nsum1+1;n++)	
		for (int k=0;k<nsum1+1;k++)	{
			
			integral += C[k]*C[n]/2/(1+m+k+n);
			
		}
	
	mpfr_class p1("0",prec_bit);
	p1 += point1;
	integral *= pow(p1,2*(1+m));


	return integral;
}


long int FuncSS::calcPrecisionForN(std::string& str_N, std::string& str_derN,std::string c_F, std::string c_energy, std::string c_Z2)	{
	
	long int crit;
	


	
	while(true)	{
		crit =  precisionN - out_lenN - getN(str_N,  str_derN, c_F,  c_energy,  c_Z2);
		if (crit>=500)	{pointN1++; pointN2 = pointN1*pointN1;}
		if (crit<=300)	{pointN1--; pointN2 = pointN1*pointN1;}
		if ((crit>300)&&(crit<500))	break;
//		cout<<"pointN1 = "<<pointN1<<"\n";
//		cout<<"crit = "<<crit<<"\n";

	}
	
	
	

//	return (int)(out_lenN*3.3219280948873626);
	return pointN1;
}



long int FuncSS::getN(std::string& str_N, std::string& str_der_N,std::string c_F, std::string c_energy, std::string c_Z2)	{


	long int i;
	std::string temp_str_N = "empty";
	std::string temp_str_der_N = "empty";
	std::stringstream out;
	out.precision(out_lenN);

    mpfr_class N("0",prec_bitN);
    mpfr_class der_N("0",prec_bitN);
	mpfr_class abs_sum_N("0",prec_bitN);

	mpfr_class F(c_F,prec_bitN);
	mpfr_class E(c_energy,prec_bitN);
	mpfr_class Z2(c_Z2,prec_bitN);


	
    mpfr_class Ci_3("0",prec_bitN);
    mpfr_class Ci_2("0",prec_bitN);
    mpfr_class Ci_1("1",prec_bitN);
    mpfr_class Ci("0",prec_bitN);
    mpfr_class Ci_max("1",prec_bitN);    
    mpfr_class sum_N("1",prec_bitN);
    mpfr_class sum_der_N("0.5",prec_bitN); sum_der_N += m;

	

		for (i=1;i<1000000000;i++)	{
			
			Ci = (-Ci_3*F*pointN2*pointN2 - Ci_1*Z2 - Ci_2*E*2*pointN2)*pointN2/4/i/(i + m);
			
			if(abs(Ci)>Ci_max)
				Ci_max = abs(Ci);
			
			
			sum_der_N += Ci*(4*i+2*m+1)/2;
			sum_N += Ci;
			
			if (i%100==0)	{
				
				
				out<<sum_N;
				str_N = out.str();
				out.str("");
				out<<sum_der_N;
				str_der_N = out.str();
				out.str("");
			
			if((str_N != temp_str_N)&&(temp_str_der_N != str_der_N))
			{
				temp_str_N = str_N;
				temp_str_der_N = temp_str_der_N;
			}
			else
				break;
			}
			
			
			Ci_3 = Ci_2;
			Ci_2 = Ci_1;
			Ci_1 = Ci;
			
			
		}
//		cout<<i<<endl;
		mpfr_class pN1("0",prec_bitN);
		pN1 += pointN1;
		
	    mpfr_class power("0.5",prec_bitN);
		power += m;
		N = sum_N*pow(pN1,power);
		der_N = sum_der_N*pow(pN1,power)/pointN1;
		
		
		out<<N;
		str_N = out.str();
		out.str("");
		out<<der_N;
		str_der_N = out.str();
		out.str("");
		
		
		Ci_max = log10(Ci_max);
		abs_sum_N = log10(abs(sum_N));
		
		return (long int)(Ci_max.get_d()) - (long int)(abs_sum_N.get_d());


}



long int FuncSS::get_ab(std::string& str_a, std::string& str_der_a,std::string& str_b, std::string& str_der_b,std::string c_F, std::string c_energy, std::string c_Z2){
	
	
	std::string temp_str_a = "empty";
	std::string temp_str_der_a = "empty";
	std::string temp_str_b = "empty";
	std::string temp_str_der_b = "empty";

	
	std::stringstream out;
	out.precision(out_lenN);

	long int i;


  mpfr_class a("0",prec_bitN);
  mpfr_class der_a("0",prec_bitN);	
  mpfr_class b("0",prec_bitN);
  mpfr_class der_b("0",prec_bitN);
  mpfr_class k1("0",prec_bitN);
  mpfr_class k2("0",prec_bitN);
  mpfr_class k3("0",prec_bitN);
  mpfr_class sum_a("0",prec_bitN);
  mpfr_class sum_der_a("0",prec_bitN);
  mpfr_class sum_b("0",prec_bitN);
  mpfr_class sum_der_b("0",prec_bitN);
  mpfr_class phasa("0",prec_bitN);
  mpfr_class temp("0",prec_bitN);
  mpfr_class si("0",prec_bitN);
  mpfr_class co("0",prec_bitN);

	mpfr_class F(c_F,prec_bitN);
	mpfr_class E(c_energy,prec_bitN);
	mpfr_class Z2(c_Z2,prec_bitN);
	

	k1 = (E*E/F - Z2)*4*pointN2;
	k2 = sqrt(F)*8*pointN1*pointN2;
	k3 = E*8*pointN1/sqrt(F);

    mpfr_class ai_3("0",prec_bitN);
    mpfr_class ai_2("0",prec_bitN);
    mpfr_class ai_1("1",prec_bitN);
    mpfr_class ai("0",prec_bitN);	
    
    mpfr_class bi_3("0",prec_bitN);
    mpfr_class bi_2("0",prec_bitN);
    mpfr_class bi_1("1",prec_bitN);
    mpfr_class bi("0",prec_bitN);
	
	ai_1 = ai_1/pointN1;
	bi_1 = bi_1/pointN1;

	sum_a = ai_1;
	sum_der_a = ai_1;
	sum_b = bi_1;
	sum_der_b = bi_1;
	

	for (i=2;i<nsumab;i++)	{
		
		ai = -(bi_3*(4*m*m-(5-2*i)*(5-2*i)) + ai_2*k3*(i-2)+bi_1*k1)/k2/(i-1);
		bi = +(ai_3*(4*m*m-(5-2*i)*(5-2*i)) - bi_2*k3*(i-2)+ai_1*k1)/k2/(i-1);
		
		sum_der_a += ai*i;
		sum_a += ai;
		sum_der_b += bi*i;
		sum_b += bi;
		
//		cout.precision(out_lenN);
//		std::cout<<i<<" "<<bi/2<<"   at  "<<ai/2<<endl;
		
		if (i%10==0)	{
//			std::cout<<"pointN1 = "<<pointN1<<" i= "<<i<<" "<<"sum_der_a = "<<mpc_get_str(10, out_lenN, sum_der_a, MPC_RNDNN)<<"  sum_a = "<<mpc_get_str(10, out_lenN, sum_a, MPC_RNDNN)<<"\n";
//			std::cout<<i<<"\n";
			
			
			out<<sum_a;
			str_a = out.str();
			out.str("");
			out<<sum_der_a;
			str_der_a = out.str();
			out.str("");
			
			out<<sum_b;
			str_b = out.str();
			out.str("");
			out<<sum_der_b;
			str_der_b = out.str();
			out.str("");
		
		if((str_a != temp_str_a)&&(temp_str_der_a != str_der_a)&&(str_b != temp_str_b)&&(temp_str_der_b != str_der_b))
		{
			temp_str_a = str_a;
			temp_str_der_a = temp_str_der_a;
			temp_str_b = str_b;
			temp_str_der_b = temp_str_der_b;
		}
		else
			break;
		}
			
		
		ai_3 = ai_2;
		ai_2 = ai_1;
		ai_1 = ai;	
		
		bi_3 = bi_2;
		bi_2 = bi_1;
		bi_1 = bi;	

	}

	phasa = sqrt(F)*pointN1*pointN2/3 + E*pointN1/sqrt(F);
	temp = sqrt(F)*pointN1*pointN2 + E*pointN1/sqrt(F);
	si = sin(phasa);
	co = cos(phasa);
	
	a = sum_a*si + sum_b*co;
	b = sum_b*si - sum_a*co;
	

	der_a = ((temp*sum_a - sum_der_b)*co - (temp*sum_b + sum_der_a)*si)/pointN1;
	der_b = ((temp*sum_b + sum_der_a)*co + (temp*sum_a - sum_der_b)*si)/pointN1;
	
	out<<a;
	str_a = out.str();
	out.str("");
	out<<der_a;
	str_der_a = out.str();
	out.str("");
	
	out<<b;
	str_b = out.str();
	out.str("");
	out<<der_b;
	str_der_b = out.str();
	out.str("");
	

	return i;


}






std::string FuncSS::getB(std::string c_F,std::string c_energy, std::string c_Z2){
	
	std::string str_B;
	std::string str_N;
	std::string str_der_N;
	std::string str_a;
	std::string str_der_a;
	std::string str_b;
	std::string str_der_b;

	
	std::stringstream out;
	out.precision(out_lenN);
	cout.precision(50);


	mpfr_class F(c_F,prec_bitN);
	mpfr_class Z2(c_Z2,prec_bitN);
	mpfr_class E(c_energy,prec_bitN);
	
	mpfr_class A("0",prec_bitN);
	mpfr_class B("0",prec_bitN);
	
	long int nn;

	nn = get_ab( str_a,  str_der_a, str_b,  str_der_b, c_F,  c_energy, c_Z2);	
	
	mpfr_class a(str_a,prec_bitN);
	mpfr_class der_a(str_der_a,prec_bitN);
	mpfr_class b(str_b,prec_bitN);
	mpfr_class der_b(str_der_b,prec_bitN);

	
	getN(str_N,str_der_N,c_F,c_energy, c_Z2);	
	mpfr_class N(str_N,prec_bitN);
	mpfr_class der_N(str_der_N,prec_bitN);
	
	A = -(-N*der_b+b*der_N)/(-b*der_a+a*der_b);
	B = -(-N*der_a+a*der_N)/(b*der_a-a*der_b);

//	std::string c_Z1;
//	out<<(4 - Z2);
//	c_Z1 = out.str();
//	cout<<"int_c++ =  "<<integrateM2(F,E,4 - Z2)<<"	";
//	cout<<point1<<"	";
//	cout<<c_Z1<<endl;
//	cout<<integrateM2(F,E,4 - Z2)<<"	"<<point1<<"	";
//	cout<<getM( c_F,  c_energy, c_Z1)<<"	";
//	out.str("");


	
	if(nn == nsumab)
	str_B="0";
	
	else	{
		
	out<<sqrt((A*A+B*B)*2);	
		
	str_B = out.str();
	out.str("");
	
	}
	
	return str_B;

}

std::string FuncSS::C(std::string c_F,std::string c_energy, std::string c_Z2){
	
	std::string str_C;
	std::string str_N;
	std::string str_der_N;
	std::string str_a;
	std::string str_der_a;
	std::string str_b;
	std::string str_der_b;

	
	std::stringstream out;
	out.precision(out_lenN);



	mpfr_class F(c_F,prec_bitN);
	mpfr_class Z2(c_Z2,prec_bitN);
	mpfr_class E(c_energy,prec_bitN);
	
	mpfr_class A("0",prec_bitN);
	mpfr_class B("0",prec_bitN);
	mpfr_class a1("0",prec_bitN);
	mpfr_class b1("0",prec_bitN);
	
	long int nn = get_ab( str_a,  str_der_a, str_b,  str_der_b, c_F,  c_energy, c_Z2);

	
	mpfr_class a(str_a,prec_bitN);
	mpfr_class der_a(str_der_a,prec_bitN);
	mpfr_class b(str_b,prec_bitN);
	mpfr_class der_b(str_der_b,prec_bitN);
	
	mpfr_class one("1",prec_bit);

	
	getN(str_N,str_der_N,c_F,c_energy, c_Z2);	
	mpfr_class N(str_N,prec_bitN);
	mpfr_class der_N(str_der_N,prec_bitN);
	
	A = -(-N*der_b+b*der_N)/(-b*der_a+a*der_b);	//Bre
	B = -(-N*der_a+a*der_N)/(b*der_a-a*der_b);	//Bim
	
	a1 = A + B;
	b1 = A - B;

	if(nn == nsumab)
	str_C="0";
	
	else	{
		
	out<<one/sqrt(integrateM2(F,E, 4 - Z2)*(a1*a1 + b1*b1)*sqrt(F)*const_pi()*const_pi());	
		
	str_C = out.str();
	out.str("");
	
	}
	
	return str_C;

}




std::string FuncSS::M(std::string c_mu,std::string c_F, std::string c_energy, std::string c_Z1){
	
	std::stringstream out;
	out.precision(out_len);
	mpfr_class mu(c_mu,prec_bit);
	
	if(mu<=point1)	{
	

	
	mpfr_class F(c_F,prec_bit);
	mpfr_class E(c_energy,prec_bit);
	mpfr_class Z1(c_Z1,prec_bit);
	
	mpfr_class mu2("0",prec_bit);
	mpfr_class mu4("0",prec_bit);
	mu2 = mu*mu;
	mu4 = mu2*mu2;

    mpfr_class Ci_3("0",prec_bit);
    mpfr_class Ci_2("0",prec_bit);
    mpfr_class Ci_1("1",prec_bit);
    mpfr_class Ci("0",prec_bit);

	
	mpfr_class sum_M("1",prec_bit);

	for (int i=1;i<nsum1+1;i++)	{
				
		Ci = (Ci_3*F*mu4 - Ci_1*Z1 - Ci_2*E*2*mu2)*mu2/4/i/(i + m);
		sum_M += Ci;
		
		Ci_3 = Ci_2;
		Ci_2 = Ci_1;
		Ci_1 = Ci;
	}


	
    mpfr_class power("0.5",prec_bit);
	power += m;
	mpfr_class M = sum_M*pow(mu,power);

	
	out<<M;
	}
	else
	out<<"0";
		
	return out.str();

}


std::string FuncSS::N(std::string c_nu,std::string c_F, std::string c_energy, std::string c_Z2)	{


	long int i;
	
	std::stringstream out;
	out.precision(out_lenN);

	mpfr_class nu(c_nu,prec_bitN);
	
	if(nu<=pointN1)		{
	
	std::string temp_str_N = "empty";
	std::string str_N = "empty";
	
    mpfr_class N("0",prec_bitN);
	mpfr_class abs_sum_N("0",prec_bitN);

	mpfr_class F(c_F,prec_bitN);
	mpfr_class E(c_energy,prec_bitN);
	mpfr_class Z2(c_Z2,prec_bitN);
	
	mpfr_class nu2("0",prec_bit);
	mpfr_class nu4("0",prec_bit);
	nu2 = nu*nu;
	nu4 = nu2*nu2;


	
    mpfr_class Ci_3("0",prec_bitN);
    mpfr_class Ci_2("0",prec_bitN);
    mpfr_class Ci_1("1",prec_bitN);
    mpfr_class Ci("0",prec_bitN);

    mpfr_class sum_N("1",prec_bitN);


	

		for (i=1;i<1000000000;i++)	{
			
			Ci = (-Ci_3*F*nu4 - Ci_1*Z2 - Ci_2*E*2*nu2)*nu2/4/i/(i + m);
			
			sum_N += Ci;
			
			if (i%100==0)	{
				
				out<<sum_N;
				str_N = out.str();
				out.str("");

			
			if(str_N != temp_str_N)
			{
				temp_str_N = str_N;
			}
			
			else
				break;
			}
			
			
			Ci_3 = Ci_2;
			Ci_2 = Ci_1;
			Ci_1 = Ci;
			
			
		}
		
	    mpfr_class power("0.5",prec_bitN);
		power += m;
		N = sum_N*pow(nu,power);

		
		
		out<<N;
		
	}
	
	else
		out<<"0";
		

		
		return out.str();


}


