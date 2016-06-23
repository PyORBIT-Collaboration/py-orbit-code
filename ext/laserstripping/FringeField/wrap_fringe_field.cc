#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_fringe_field.hh"

#include <iostream>
#include <string>

#include "FringeField.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_fringe_field{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* FringeField_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int FringeField_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  

	  

	  double d;
	  double F;
	  int der;

		 if(!PyArg_ParseTuple(	args,"ddi:",&d, &F, &der)){
			 		          error("FringeField(gap , By_max, derivative of module |By(z)| (plus one or minus one)) - params. are needed");
			 			 		        }  
		 else	{

		 self->cpp_obj =  new  FringeField(d, F, der);
		 ((FringeField*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

		
    return 0;
  }
  
  
  
 
  
  static PyObject* FringeField_getField(PyObject *self, PyObject *args){
	  FringeField* cpp_fields = (FringeField*)((pyORBIT_Object*) self)->cpp_obj;
  				       
  		
       int nVars = PyTuple_Size(args);
       
       double x;
       double y;
       double z;
       double t;
       double Ex;
       double Ey;
       double Ez;
       double Bx;
       double By;
       double Bz;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"dddd:",&x,&y,&z,&t))
             error(" getField(x,y,z,t) - parameters are needed");
           else
           cpp_fields->getElectricMagneticField(x,y,z,t,Ex,Ey,Ez,Bx,By,Bz);
           return Py_BuildValue("d",Bx);

  }
  
  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void FringeField_del(pyORBIT_Object* self){
		//std::cerr<<"The FringeField __del__ has been called!"<<std::endl;
		delete ((FringeField*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef FringeFieldClassMethods[] = {
/*		    { "setFieldOrientation",  FringeField_setFieldOrientation,         METH_VARARGS,"Sets field orientation."},*/
			{ "getField",         		FringeField_getField,         		 	METH_VARARGS,"Returns By field."},
/*			{ "setLaserPower",         HermiteGaussianLFmode_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
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
	static PyMemberDef FringeFieldClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_FringeField_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"FringeField", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) FringeField_del , /*tp_dealloc*/
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
		"The FringeField python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FringeFieldClassMethods, /* tp_methods */
		FringeFieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) FringeField_init, /* tp_init */
		0, /* tp_alloc */
		FringeField_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initFringeField(PyObject* module){
		if (PyType_Ready(&pyORBIT_FringeField_Type) < 0) return;
		Py_INCREF(&pyORBIT_FringeField_Type);
		PyModule_AddObject(module, "FringeField", (PyObject *)&pyORBIT_FringeField_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
