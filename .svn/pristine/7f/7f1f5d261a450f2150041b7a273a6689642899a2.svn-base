#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_froissart_stora_lf.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "FroissartStoraLF.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_froissart_stora_lf{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* FroissartStoraLF_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int FroissartStoraLF_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
	  double Omega;
	  double Gamma;
	  double Elas;
	  
		 if(!PyArg_ParseTuple(	args,"ddd:",&Omega,&Gamma,&Elas)){
			 		          error("FroissartStoraLF(Omega,Gamma,ampl_Elas) - params. are needed");
		 } 
		 else	{
			 
		self->cpp_obj = new FroissartStoraLF(Omega,Gamma,Elas);
		 ((FroissartStoraLF*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }

    return 0;
  }
  
  
  
  
  
  
  static PyObject* FroissartStoraLF_setLaserFieldPolarization(PyObject *self, PyObject *args){
	  FroissartStoraLF* LaserField = (FroissartStoraLF*)((pyORBIT_Object*) self)->cpp_obj;
 	 

   		double nEx;
   		double nEy;
   		double nEz;


   		        if(!PyArg_ParseTuple(	args,"ddd:", &nEx, &nEy, &nEz))
   		          {error("LaserExternalEfects - setLaserFieldPolarization(nEx, nEy, nEz) - params. afe needed");}
  		        
   		        else	
   		        LaserField->setLaserFieldPolarization(nEx, nEy, nEz);
   		      
   	   
   		    Py_INCREF(Py_None);
   		    return Py_None;	  
    }	
  
		
  
  
  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void FroissartStoraLF_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((FroissartStoraLF*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef FroissartStoraLFClassMethods[] = {
		    { "setLaserFieldPolarization",  FroissartStoraLF_setLaserFieldPolarization,         METH_VARARGS,"Sets or returns the name of effects."},
/*			{ "setLaserHalfAngle",         HermiteGaussianLFmode_setLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserPower",         HermiteGaussianLFmode_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaser_lambda",         HermiteGaussianLFmode_setLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserHalfAngle",         HermiteGaussianLFmode_getLaserHalfAngle,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaserPower",         HermiteGaussianLFmode_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaser_lambda",         HermiteGaussianLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  HermiteGaussianLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef FroissartStoraLFClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_FroissartStoraLF_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"FroissartStoraLF", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) FroissartStoraLF_del , /*tp_dealloc*/
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
		"The FroissartStoraLF python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FroissartStoraLFClassMethods, /* tp_methods */
		FroissartStoraLFClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) FroissartStoraLF_init, /* tp_init */
		0, /* tp_alloc */
		FroissartStoraLF_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initFroissartStoraLF(PyObject* module){
		if (PyType_Ready(&pyORBIT_FroissartStoraLF_Type) < 0) return;
		Py_INCREF(&pyORBIT_FroissartStoraLF_Type);
		PyModule_AddObject(module, "FroissartStoraLF", (PyObject *)&pyORBIT_FroissartStoraLF_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
