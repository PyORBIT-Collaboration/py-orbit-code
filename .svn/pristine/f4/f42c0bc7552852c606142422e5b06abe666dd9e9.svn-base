//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppBaseFieldSource.cc
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

//#include <fstream>
#include "orbit_mpi.hh"
//#include <iostream>
////#include <iomanip>
//#include <cmath>
#include <fstream>
#include <string>

#include "SNSstrippingMagnet.hh"



using namespace OrbitUtils;

SNSstrippingMagnet::SNSstrippingMagnet(std::string field_data)
{      
        pos = 0;
        len = 0;
	std::ifstream file;      
	file.open(field_data.c_str());
        std::string line;        
        while (std::getline(file,line)) { len++;}
        file.close();
        
        z = new double [len];
        Bx = new double[len];

        
        file.open(field_data.c_str());               
        for(int i = 0; i < len; i++){file>>z[i]>>Bx[i];}                
        file.close();             
        

        A = 1;
        l = z[len- 1] - z[0];
        dz = l/(len - 1);           

}




SNSstrippingMagnet::~SNSstrippingMagnet()
{ 

    delete [] z;
    delete [] Bx;

}



double SNSstrippingMagnet::getBx(double zi){
    
	
    if (zi <= z[0] || zi >= z[len - 1])	{return 0; }
    else{
    int i = (int)((zi - z[0])/dz); 
        return Bx[i] + (Bx[i+1] - Bx[i])*(zi - z[i])/dz;
        }
}



void SNSstrippingMagnet::setA(double _A){A = _A;}
void SNSstrippingMagnet::setPosition(double _pos){	pos = _pos;}

void SNSstrippingMagnet::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& B_x, double& B_y, double& B_z)

{       
        
        B_x = getBx(z - pos);
        B_y = 0;
        B_z = 0;
        E_x = 0;
        E_y = 0;
        E_z = 0;



return;

}




 