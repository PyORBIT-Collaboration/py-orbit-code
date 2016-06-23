#include "Bunch.hh"
#include "Nll.hh"

#include "ParticleAttributesFactory.hh"
#include "OrbitConst.hh"
#include "StringUtils.hh"
#include "BufferStore.hh"

#include <iomanip>
#include <string>
#include <stdio.h>
#include <math.h>

using namespace OrbitUtils;

//Nll::Nll(): CppPyWrapper(NULL){
Nll::Nll(){}

void Nll::TRACK_EXT (Bunch* bunch){
	
	double L = 1;
	double N = 100;
	double LS = .4;
	double C = 1;	
	
  double dd, u, v, dUu, dUv, dux, duy, dvx, dvy, dt, k, x, y, px, py, pi = 3.141592653589;
  
	//Non-linear Lens type 3
	//U(u,v)=(u*sqrt(u**2-1)+v*sqrt(1-v**2)*(-pi/2+acos(v))/(u**2-v**2);
	//L = length, N = number of slices, LS = Lens Strength, c = distance between poles

	for(int ip = 0, nParts = bunch->getSize(); ip < nParts; ip++){
		x = bunch->x(ip);
		y = bunch->y(ip);
		px = bunch->px(ip);
		py = bunch->py(ip);
		
		u=0.5*sqrt(pow((x-1),2)+pow(y,2))+0.5*sqrt(pow((x+1),2)+pow(y,2));
		v=0.5*sqrt(pow((x+1),2)+pow(y,2))-0.5*sqrt(pow((x-1),2)+pow(y,2));
		if (u==1) {
			dd=0;
		}else dd=pow(u,2)*acosh(u)/sqrt(pow(u,2)-1);
		
		dUu=(u+acosh(u)*sqrt(pow(u,2)-1)+dd)/(pow(u,2)-pow(v,2))-2*u*(u*acosh(u)*sqrt(pow(u,2)-1)+v*(acos(v)-0.5*pi)*sqrt(1-pow(v,2)))/pow((pow(u,2)-pow(v,2)),2);
		dUv=2*v*(u*acosh(u)*sqrt(pow(u,2)-1)+v*(acos(v)-0.5*pi)*sqrt(1-pow(v,2)))/pow((pow(u,2)-pow(v,2)),2)-(v-(acos(v)-0.5*pi)*sqrt(1-pow(v,2))+pow(v,2)*(acos(v)-0.5*pi)/sqrt(1-pow(v,2)))/(pow(u,2)-pow(v,2));
		dux=0.5*(x-1)/sqrt(pow((x-1),2)+pow(y,2))+0.5*(x+1)/sqrt(pow((x+1),2)+pow(y,2));
		duy=0.5*y/sqrt(pow((x-1),2)+pow(y,2))+0.5*y/sqrt(pow((x+1),2)+pow(y,2));
		dvx=0.5*(x+1)/sqrt(pow((x+1),2)+pow(y,2))-0.5*(x-1)/sqrt(pow((x-1),2)+pow(y,2));
		dvy=0.5*y/sqrt(pow((x+1),2)+pow(y,2))-0.5*y/sqrt(pow((x-1),2)+pow(y,2));
		
		dt = L/N;
		k = LS/C;
		
		x = x; //+ dt*(px - dt*(k*(dUu*dux+dUv*dvx)));
		y = y; //+ dt*(py - dt*(k*(dUu*duy+dUv*dvy)));
		px = px - dt*(k*(dUu*dux+dUv*dvx));
		py = py - dt*(k*(dUu*duy+dUv*dvy));
		
		bunch->x(ip) = x;
		bunch->y(ip) = y;
		bunch->px(ip) = px;
		bunch->py(ip) = py;
		
	}

}

