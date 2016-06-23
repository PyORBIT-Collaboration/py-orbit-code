#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_cone_lf_mode.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "ConeLFmode.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_cone_lf_mode{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* ConeLFmode_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int ConeLFmode_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
	  double P=0;
	  double rx=0;
	  double ry=0;
	  double a_x=0;
	  double a_y=0;
	  double lambda=0;
	  
		 if(!PyArg_ParseTuple(	args,"dddddd:",&P,&rx,&ry,&a_x,&a_y,&lambda)){
			 		          error("LaserExternalEfects(P,rx,ry,ax,ay,lambda) - params. are needed");
		 } else	{
			 
		self->cpp_obj = new ConeLFmode(P,rx,ry,a_x,a_y,lambda);
		((ConeLFmode*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }

    return 0;
  }
  
  
  
  
  
  
  static PyObject* ConeLFmode_setLaserFieldOrientation(PyObject *self, PyObject *args){
	  
 	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
   		double x0;
   		double y0;
   		double z0;
   		double kx;
   		double ky;
   		double kz;
   		double mx;
   		double my;
   		double mz;
   		double nEx;
   		double nEy;
   		double nEz;
   		double flag=true;

   		        if(!PyArg_ParseTuple(	args,"dddddddddddd:",&x0, &y0, &z0, &kx, &ky, &kz, &mx, &my, &mz, &nEx, &nEy, &nEz))
   		          {error("LaserExternalEfects - setLaserHalfAngle(comments is not written yet) - param. is needed");flag=false;}
  		        
   		        if(fabs(kx*mx+ky*my+kz*mz)/sqrt(kx*kx+ky*ky+kz*kz)/sqrt(mx*mx+my*my+mz*mz)>1.e-15)
   		          {error("Please be shure that kx*mx+ky*my+kz*mz==0");flag=false;}
   		        if(fabs(kx*nEx+ky*nEy+kz*nEz)/sqrt(kx*kx+ky*ky+kz*kz)/sqrt(nEx*nEx+nEy*nEy+nEz*nEz)>1.e-15)
   		          {error("Please be shure that kx*Ex+ky*Ey+kz*Ez==0");flag=false;}
	

   		        if(flag)	
   		        	LaserField->setLaserFieldOrientation(x0,y0, z0, kx, ky, kz, mx, my, mz, nEx, nEy, nEz);
   		      
   	   
   		    Py_INCREF(Py_None);
   		    return Py_None;	  
    }	
  
		
  
  /*
 	void	LasStripExternalEffects::setLaserHalfAngle(double a);
 	double	LasStripExternalEffects::getLaserHalfAngle();
 	void	LasStripExternalEffects::setLaserPower(double a);
 	double	LasStripExternalEffects::getLaserPower();
 	void	LasStripExternalEffects::setLaser_lambda(double a);
 	double	LasStripExternalEffects::getLaser_lambda();

   
 static PyObject* ConeLFmode_setLaserHalfAngle(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
  		  double half_Angle;	       
  		        if(!PyArg_ParseTuple(	args,"d:",&half_Angle)){
  		          error("LaserExternalEfects - setLaserHalfAngle(half_Angle) - param. is needed");
  		        }
  		        else LaserField->setLaserHalfAngle(half_Angle);
  		      
  	   
  		    Py_INCREF(Py_None);
  		    return Py_None;	  
   }		

 static PyObject* ConeLFmode_setLaserPower(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 		  double LaserPower;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&LaserPower)){
 		          error("LaserExternalEfects - setLaserPower(LaserPower) - param. is needed");
 		        }
 		        else LaserField->setLaserPower(LaserPower);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
 }

 static PyObject* ConeLFmode_setLaser_lambda(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 		  double Laser_lambda;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&Laser_lambda)){
 		          error("LaserExternalEfects - setLaser_lambda(Laser_lambda) - param. is needed");
 		        }
 		        else LaserField->setLaserPower(Laser_lambda);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
 }
 
 
 
 static PyObject* ConeLFmode_setLaserFieldOrientation(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 		  double Laser_lambda;	       
 		        if(!PyArg_ParseTuple(	args,"d:",&Laser_lambda)){
 		          error("LaserExternalEfects - setLaser_lambda(Laser_lambda) - param. is needed");
 		        }
 		        else LaserField->setLaserPower(Laser_lambda);
 		      
 	   
 		    Py_INCREF(Py_None);
 		    return Py_None;	  
 }




 static PyObject* ConeLFmode_getLaserHalfAngle(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaserHalfAngle()); 
 }		
     

 static PyObject* ConeLFmode_getLaserPower(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaserPower()); 
 }		    

 static PyObject* ConeLFmode_getLaser_lambda(PyObject *self, PyObject *args){
	 ConeLFmode* LaserField = (ConeLFmode*)((pyORBIT_Object*) self)->cpp_obj;
 				       
 		  return Py_BuildValue("d",LaserField->getLaser_omega()); 
 }		    

   */

  
  
  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void ConeLFmode_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((ConeLFmode*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef ConeLFmodeClassMethods[] = {
		    { "setLaserFieldOrientation",  ConeLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
/*			{ "setLaserHalfAngle",         ConeLFmode_setLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserPower",         ConeLFmode_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaser_lambda",         ConeLFmode_setLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserHalfAngle",         ConeLFmode_getLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserPower",         ConeLFmode_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaser_lambda",         ConeLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  ConeLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef ConeLFmodeClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_ConeLFmode_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ConeLFmode", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ConeLFmode_del , /*tp_dealloc*/
		0, /*tp_print*/
		0, /*tp_getattr*/
		0, /*tp_setattr*/
		0, /*tp_compare*/
		0, /*tp_repr*/
		0, /*tp_as_number*/
		0, /*tp_as_sequence*/
		0, /*tp_as_mapping*/
		0, /*tp_hash */
		0, /*tp_call*/
		0, /*tp_str*/
		0, /*tp_getattro*/
		0, /*tp_setattro*/
		0, /*tp_as_buffer*/
		Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
		"The ConeLFmode python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ConeLFmodeClassMethods, /* tp_methods */
		ConeLFmodeClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ConeLFmode_init, /* tp_init */
		0, /* tp_alloc */
		ConeLFmode_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initConeLFmode(PyObject* module){
		if (PyType_Ready(&pyORBIT_ConeLFmode_Type) < 0) return;
		Py_INCREF(&pyORBIT_ConeLFmode_Type);
		PyModule_AddObject(module, "ConeLFmode", (PyObject *)&pyORBIT_ConeLFmode_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
