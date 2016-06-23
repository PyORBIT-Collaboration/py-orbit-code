//////////////////////////////// -*- C++ -*- //////////////////////////////
//
// FILE NAME
//    CppExternalEffects.cc
//
// CREATED
//    06/27/2008
//
// DESCRIPTION
//    The base class for C++ implementation of a external effects 
//    during the transport of particles through the external field. 
//    It should be sub-classed on Python level and implement
//    setupEffects(Bunch* bunch)
//    finalizeEffects(Bunch* bunch)
//    applyEffects(Bunch* bunch, int index, 
//	                            double* y_in_vct, double* y_out_vct, 
//														  double t, double t_step, 
//														  OrbitUtils::BaseFieldSource* fieldSource)
//    methods.
//    The results of these methods will be available from the c++ level.
//    This is an example of embedding Python in C++ Orbit level.
//
//     
//         
//       
//     
//
///////////////////////////////////////////////////////////////////////////
#include "orbit_mpi.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>



#include "Walls.hh"
#include "RungeKuttaTracker.hh"







#define x0(i) 	 Coords->attArr(i)[0]
#define y0(i)  	 Coords->attArr(i)[2]
#define z0(i)  	 Coords->attArr(i)[4]

#define x(i) 	 bunch->coordArr()[i][0]
#define y(i)  	 bunch->coordArr()[i][2]
#define z(i)  	 bunch->coordArr()[i][4]

#define px(i) 	 bunch->coordArr()[i][1]
#define py(i)  	 bunch->coordArr()[i][3]
#define pz(i)  	 bunch->coordArr()[i][5]





using namespace LaserStripping;
using namespace OrbitUtils;


Walls::Walls()
{

}




Walls::~Walls() 
{

}
















void Walls::setupEffects(Bunch* bunch){
	

	
	if (bunch->hasParticleAttributes("pq_coords")==0)	{
	std::map<std::string,double> part_attr_dict;
	part_attr_dict["size"] = 6;
	bunch->addParticleAttributes("pq_coords",part_attr_dict);
	}

	Coords = bunch->getParticleAttributes("pq_coords");


}



void Walls::memorizeInitParams(Bunch* bunch){
	
	for (int i=0; i<bunch->getSize();i++)		{
		x0(i)=bunch->coordArr()[i][0];
		y0(i)=bunch->coordArr()[i][2];
		z0(i)=bunch->coordArr()[i][4];

	}
	
	
}
		
	
void Walls::finalizeEffects(Bunch* bunch) {
	
	if(bunch->hasParticleAttributes("pq_coords")==1)
		bunch->removeParticleAttributes("pq_coords");
	
}





/*

bool Walls::crossPlane(int i,Bunch* bunch,double xs, double ys,double zs,double vx, double vy,double vz){
	
	double vr = (x0(i)-x(i))*vx+(y0(i)-y(i))*vy+(z0(i)-z(i))*vz;
	
	
	if(vr==0)
		return false;
	
	else
	{
		double coeff = (vx*(xs-x0(i))+vy*(ys-y0(i))+vz*(zs-z0(i)))/vr;
		double x_cross = x0(i)+(x0(i)-x(i))*coeff;
		double y_cross = y0(i)+(y0(i)-y(i))*coeff;
		double z_cross = z0(i)+(z0(i)-z(i))*coeff;
		
		if ((x0(i)<=x_cross<=x(i))&&(y0(i)<=y_cross<=y(i))&&(z0(i)<=z_cross<=z(i)))
			return true;
		else
			return false;
		
		
	}
	

}

*/



bool Walls::crossSurface(int i,Bunch* bunch)	{

double ax = -0.310729;
double bx = 11.9591;
double ay = -0.218678;
double by = 13.4845;

double gammax = (1+ax*ax)/bx;
double gammay = (1+ay*ay)/by;

double ex = 240.0e-6;
double ey = 240.0e-6;

double x_bump = 0.015;
double y_bump = 0.016;

double xf = -0.09679;			//position of input flange
double yf = 0;					//position of input flange
double zf = 4.430187;			//position of input flange

double x0 = xf+x_bump;			//position of Cl orb bmp end at the emittance aceptance point
double y0 = yf+y_bump;			//position of Cl orb bmp end at the emittance aceptance point
double z0 = zf;					//position of Cl orb bmp end at the emittance aceptance point





double xp = px(i)/pz(i);
double yp = py(i)/pz(i);

double x1 = x(i)-xf;
double y1 = y(i)-yf;

double x = x(i)-x0;
double y = y(i)-y0;



double R = 0.093115;


/*
cout<<x(i)<<" "<<z(i)<<"\n";

if(((z0(i)<=2.936809)&&(2.936809<=z(i))&&(0.118921>=x(i))&&(x(i)>=0.026719)))	cout<<"foil"<<"\n";
	

if(((z0(i)<=z0)&&(z0<=z(i))		//intersection of flange z plane
				
		&&
				
		(x1*x1+y1*y1<=R*R)		//position of ring flange 
					
		&&
		
		(x*x*gammax+2*ax*x*xp+bx*xp*xp < ex)&&(y*y*gammay+2*ay*y*yp+by*yp*yp < ey)))	cout<<"ring"<<"\n";
*/

if(
	
		((z0(i)<=2.936809)&&(2.936809<=z(i))&&(0.118921>=x(i))&&(x(i)>=0.026719))		//position of thick foil
		
		||
		
//		((z0(i)<=4.442887)&&(4.442887<=z(i))&&((x(i)+0.09679)*(x(i)+0.09679)+(y(i)-0)*(y(i)-0)<=0.0103476))	//position of ring flange
		
		((z0(i)<=z0)&&(z0<=z(i))		//intersection of flange z plane
				
		&&
				
		(x1*x1+y1*y1<=R*R)		//position of ring flange 
					
		&&
		
		(x*x*gammax+2*ax*x*xp+bx*xp*xp < ex)&&(y*y*gammay+2*ay*y*yp+by*yp*yp < ey))	//emittance of the ring beam 
		


//		((z0(i)<=2.936809)&&(2.936809<=z(i))&&(0.118921>=x(i))&&(x(i)>=0.026719))
	
		
//		((z0(i)<=0.25)&&(0.25<=z(i)))
		

)
	

		
	
	return true;
else
	return false;

}






void Walls::applyEffects(Bunch* bunch,
														  double t, double t_step, 
														  BaseFieldSource* fieldSource,
															RungeKuttaTracker* tracker)			{



		for (int i=0; i<bunch->getSize();i++)	{

			if(crossSurface(i,bunch))
				bunch->deleteParticleFast(i);
		}	

	
}






