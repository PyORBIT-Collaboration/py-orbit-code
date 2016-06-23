#include "FieldRotation.hh"




double	FieldRotation::RotateElectricFieldsV(double vx,double vy,double vz, tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl){
	
	double v = sqrt(vx*vx+vy*vy+vz*vz);
	
	double nx = vx/v;
	double ny = vy/v;
	double nz = vz/v;
	
	double n = 1./(nz + 1.);
	double nxy = -nx*ny*n;


		tcomplex Exll=Exl;
		tcomplex Eyll=Eyl;
		tcomplex Ezll=Ezl;
		
					
		Exl=Exll*(1-nx*nx*n) + Eyll*nxy - Ezll*nx;
		Eyl=Exll*nxy + Eyll*(1-ny*ny*n) - Ezll*ny;
		Ezl=Exll*nx + Eyll*ny + Ezll*nz;

return v;

}



void	FieldRotation::RotateElectricFieldsN(double nx,double ny,double nz, tcomplex& Exl,tcomplex& Eyl,tcomplex& Ezl){
	

	double n = 1/(nz + 1);
	double nxy = -nx*ny*n;


		tcomplex Exll=Exl;
		tcomplex Eyll=Eyl;
		tcomplex Ezll=Ezl;
		
					
		Exl=Exll*(1-nx*nx*n)+ Eyll*nxy-Ezll*nx;
		Eyl=Exll*nxy+ Eyll*(1-ny*ny*n)-Ezll*ny;
		Ezl=Exll*nx+ Eyll*ny+Ezll*nz;


return;

}

