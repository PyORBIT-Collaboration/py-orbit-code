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
#include "RegularGridFS.hh"

#include "orbit_mpi.hh"
#include "BufferStore.hh"

#include <cmath>
#include <fstream>
#include <iostream>


using namespace OrbitUtils;

RegularGridFS::RegularGridFS(std::string file_name,std::string type_of_field, double coord_mult, double field_mul)
{
	orient = new FieldOrientation();
	
	if (type_of_field[0]=='e')
		e_field=true;
	if (type_of_field[0]=='m')
		e_field=false;

	
	xyz_mult=coord_mult;
	field_mult = field_mul;
	
	
	int rank_MPI,size_MPI;
	ORBIT_MPI_Comm_size(MPI_COMM_WORLD, &size_MPI);
	ORBIT_MPI_Comm_rank(MPI_COMM_WORLD, &rank_MPI);
	
	std::ifstream file;
	double F,x,y,z;
	int fi;
	

	if(rank_MPI == 0) {
		file.open(file_name.c_str());	file>>x_min>>y_min>>z_min>>F>>F>>F;  file.close();
		
		file.open(file_name.c_str());	while(true)	{file>>x>>y>>z>>F>>F>>F;if (x>x_min+fabs(x_min)*1.0e-14) {delta_x=x-x_min; break;}} file.close();
		file.open(file_name.c_str());	while(true)	{file>>x>>y>>z>>F>>F>>F;if (y>y_min+fabs(y_min)*1.0e-14) {delta_y=y-y_min; break;}} file.close();
		file.open(file_name.c_str());	while(true)	{file>>x>>y>>z>>F>>F>>F;if (z>z_min+fabs(z_min)*1.0e-14) {delta_z=z-z_min; break;}} file.close();		
		
		file.open(file_name.c_str());fi=0;	while(!file.eof())	{file>>x_max>>y_max>>z_max>>F>>F>>F;fi++;} file.close();	n_data=fi-1;
		
		
		x_min *=xyz_mult; y_min *=xyz_mult; z_min *=xyz_mult;
		x_max *=xyz_mult; y_max *=xyz_mult; z_max *=xyz_mult;
		delta_x *=xyz_mult; delta_y *=xyz_mult; delta_z *=xyz_mult;
	}


  ORBIT_MPI_Bcast(&n_data, 1, MPI_INT,0,MPI_COMM_WORLD);	
  
  ORBIT_MPI_Bcast(&x_min, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&y_min, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&z_min, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  ORBIT_MPI_Bcast(&x_max, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&y_max, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&z_max, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  ORBIT_MPI_Bcast(&delta_x, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&delta_y, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  ORBIT_MPI_Bcast(&delta_z, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  nx=int((x_max-x_min)/delta_x+0.5)+1;
  ny=int((y_max-y_min)/delta_y+0.5)+1;
  nz=int((z_max-z_min)/delta_z+0.5)+1;
  

  nx_min = int(x_min/delta_x*(1+0.5/fabs(x_min/delta_x)));
  ny_min = int(y_min/delta_y*(1+0.5/fabs(y_min/delta_y)));
  nz_min = int(z_min/delta_z*(1+0.5/fabs(z_min/delta_z)));
  


  
	B_x=new double**[nx+1];	
	for (int i=0;i<nx+1;i++)	B_x[i]=new double*[ny+1]; 
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	B_x[i][j]=new double[nz+1];
	
	B_y=new double**[nx+1];	
	for (int i=0;i<nx+1;i++)	B_y[i]=new double*[ny+1]; 
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	B_y[i][j]=new double[nz+1];
	
	B_z=new double**[nx+1];	
	for (int i=0;i<nx+1;i++)	B_z[i]=new double*[ny+1]; 
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	B_z[i][j]=new double[nz+1];



	int buff_index = 0;
	double* dump_arr = BufferStore::getBufferStore()->getFreeDoubleArr(buff_index,nz);
 


	//this loop reads magnetic fields
	if(rank_MPI == 0) {
		file.open(file_name.c_str()); 
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
				for(int k=0;k<nz;k++)	{

					file>>x>>y>>z>>B_x[i][j][k]>>B_y[i][j][k]>>B_z[i][j][k];
					B_x[i][j][k] *= field_mult; B_y[i][j][k]*= field_mult; B_z[i][j][k]*= field_mult;
				}
			}
		}
		file.close();
	}
	

	

	for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
					for(int k=0;k<nz;k++)	{
						dump_arr[k]=B_x[i][j][k];
					}
					
	ORBIT_MPI_Bcast(dump_arr, nz, MPI_DOUBLE,0,MPI_COMM_WORLD);
					
					for(int k=0;k<nz;k++)	{
						B_x[i][j][k]=dump_arr[k];
					}
				}
			}
	
	for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
					for(int k=0;k<nz;k++)	{
						dump_arr[k]=B_y[i][j][k];
					}
					
	ORBIT_MPI_Bcast(dump_arr, nz, MPI_DOUBLE,0,MPI_COMM_WORLD);
					
					for(int k=0;k<nz;k++)	{
						B_y[i][j][k]=dump_arr[k];
					}
				}
			}
	
	for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
					for(int k=0;k<nz;k++)	{
						dump_arr[k]=B_z[i][j][k];
					}
					
	ORBIT_MPI_Bcast(dump_arr, nz, MPI_DOUBLE,0,MPI_COMM_WORLD);
					
					for(int k=0;k<nz;k++)	{
						B_z[i][j][k]=dump_arr[k];
					}
				}
			}
	
	
	BufferStore::getBufferStore()->setUnusedDoubleArr(buff_index);	

}




RegularGridFS::~RegularGridFS()
{ 
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	delete [] B_x[i][j]; for (int i=0;i<nx+1;i++)	delete [] B_x[i];	delete	[]	B_x;
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	delete [] B_y[i][j]; for (int i=0;i<nx+1;i++)	delete [] B_y[i];	delete	[]	B_y;
	for (int i=0;i<nx+1;i++) for (int j=0;j<ny+1;j++)	delete [] B_z[i][j]; for (int i=0;i<nx+1;i++)	delete [] B_z[i];	delete	[]	B_z;

}


void RegularGridFS::setFieldOrientation(double x0, double y0, double z0,
													double kx, double ky, double kz,
													double mx, double my, double mz)	
{
	orient->setCoefficients(x0, y0, z0,kx, ky, kz, mx, my, mz);	
	return;	
}




double RegularGridFS::getNonOrientedBn(double*** B_n,double x, double y, double z, double t){
	

	if ((x_min<= x&&x<= x_max)&&(y_min<= y&&y <= y_max)&&(z_min<= z&&z <= z_max))	{
		
	

		double x_t = x/delta_x;
		double y_t = y/delta_y;
		double z_t = z/delta_z;
				
		int nx = int(x_t);
		int ny = int(y_t);
		int nz = int(z_t);
				
/*		
		for (int i=0; i<2; i++)
			for (int j=0; j<2; j++)
				for (int k=0; k<2; k++)
					
		sum += pow(-1,i+j+k)*B_n[nx-nx_min+i][ny-ny_min+j][nz-nz_min+k]*(x_t-(nx+1-i))*(y_t-(ny+1-j))*(z_t-(nz+1-k));
*/		return 

		-B_n[nx-nx_min][ny-ny_min][nz-nz_min]*(x_t-(nx+1))*(y_t-(ny+1))*(z_t-(nz+1))
		+B_n[nx-nx_min][ny-ny_min][nz-nz_min+1]*(x_t-(nx+1))*(y_t-(ny+1))*(z_t-nz)
		+B_n[nx-nx_min][ny-ny_min+1][nz-nz_min]*(x_t-(nx+1))*(y_t-ny)*(z_t-(nz+1))
		-B_n[nx-nx_min][ny-ny_min+1][nz-nz_min+1]*(x_t-(nx+1))*(y_t-ny)*(z_t-nz)
		+B_n[nx-nx_min+1][ny-ny_min][nz-nz_min]*(x_t-nx)*(y_t-(ny+1))*(z_t-(nz+1))
		-B_n[nx-nx_min+1][ny-ny_min][nz-nz_min+1]*(x_t-nx)*(y_t-(ny+1))*(z_t-nz)
		-B_n[nx-nx_min+1][ny-ny_min+1][nz-nz_min]*(x_t-nx)*(y_t-ny)*(z_t-(nz+1))
		+B_n[nx-nx_min+1][ny-ny_min+1][nz-nz_min+1]*(x_t-nx)*(y_t-ny)*(z_t-nz);

	}
	
	else
		return 0;
		
	
	
	

}


void RegularGridFS::getConstructedNonOrientedField(double& H_x,double& H_y,double& H_z, 
								double*** B_x,double***  B_y,double***  B_z, 
								double x,double y,double z,double& t){
	
	if (z>=0)	{
	
	 H_x=getNonOrientedBn(B_x,x,y,z,t);
	 H_y=getNonOrientedBn(B_y,x,y,z,t);
	 H_z=getNonOrientedBn(B_z,x,y,z,t);
	 
	}
	
	if (z<0)	{
	
	 H_x = getNonOrientedBn(B_x,x,y,-z,t);
	 H_y = getNonOrientedBn(B_y,x,y,-z,t);
	 H_z =-getNonOrientedBn(B_z,x,y,-z,t);
	 
	}
	
	
	return;
	
}
		



void RegularGridFS::getElectricMagneticField(double x, double y, double z, double t, 
		double& E_x, double& E_y, double& E_z,
		double& H_x, double& H_y, double& H_z)

{
	
	double F_x, F_y, F_z;

	 orient->OrientCoordinates(x,y,z);
	 getConstructedNonOrientedField(F_x, F_y, F_z, B_x, B_y, B_z, x, y, z, t);	 	
	 orient->OrientVector(F_x,F_y,F_z);
	 

	
	if (e_field)	 { H_x = 0.0; H_y = 0.0; H_z = 0.0; E_x = F_x; E_y = F_y; E_z = F_z;}
	if (!e_field)	 { E_x = 0.0; E_y = 0.0; E_z = 0.0; H_x = F_x; H_y = F_y; H_z = F_z;}

	
//    std::cout<< e_field<<E_x<<"  "<<E_y<<"  "<<E_z<<"  "<<H_x<<"  "<<H_y<<"  "<<H_z<<"\n";	 


return;

}




