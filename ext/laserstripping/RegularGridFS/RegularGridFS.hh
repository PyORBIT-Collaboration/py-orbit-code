//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    PyBaseFieldSource.hh
//
// CREATED
//    04/21/2003
//
// DESCRIPTION
//    The base class for Python implementation of a field source. 
//    It should be sub-classed on Python level and implement 
//    getElectricField(x,y,z,t) and getMagneticField (x,y,z,t) methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//    ATTENTION: Using this class in real calculations is not wise! 
//               It is slow, because it delegates the field calculations
//               to the Python level. It is for prototyping and 
//               debugging only.
//
///////////////////////////////////////////////////////////////////////////
#ifndef REGULARGRIDFS_HH_
#define REGULARGRIDFS_HH_

#include "Python.h"
#include "BaseFieldSource.hh"
#include "FieldOrientation.hh"
#include <vector>
#include <string>

namespace OrbitUtils{
	
	class  RegularGridFS: public BaseFieldSource
	{
		public:
		
			RegularGridFS(std::string file_name, std::string type_of_field, double coord_mult, double field_mul);
			~RegularGridFS();
			
			void setFieldOrientation(double x_0, double y_0, double z_0,
										double k_x, double k_y, double k_z,
										double m_x, double m_y, double m_z);
		
			void getElectricMagneticField(double x, double y, double z, double t, double& E_x, double& E_y, double& E_z, double& H_x, double& H_y, double& H_z);

			void getConstructedNonOrientedField(double& H_x,double& H_y,double& H_z, 
											double*** B_x,double***  B_y,double***  B_z, 
											double x,double y,double z,double& t);
			
			double getNonOrientedBn(double*** B_n,double x, double y, double z, double t);
			
		private:
			
			 
			 FieldOrientation*   orient; 
			
			 double*** B_x;
			 double*** B_y;
			 double*** B_z;

			  
			 double delta_x;
			 double delta_y;
			 double delta_z;
			 
			 double x_min;
			 double y_min;
			 double z_min;
			 		 
			 double x_max;
			 double y_max;
			 double z_max;
			 
			 int n_data;
			 
			 int nx;
			 int ny;
			 int nz;
			 
			 int nx_min;
			 int ny_min;
			 int nz_min;
			 
			 double Ez_stat;
			 int st;
			 int levels;
			 
			 bool e_field;
			 double xyz_mult;
			 double field_mult;

			
	};
};




#endif /*REGULARGRIDFS_HH_*/



