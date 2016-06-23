#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_stark_strong_field.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "StarkStrongField.hh"

//using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_stark_strong_field{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* StarkStrongField_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int StarkStrongField_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  int n1,n2,m;
	  const char* addressEG;

		 if(!PyArg_ParseTuple(	args,"siii:",&addressEG,&n1,&n2,&m)){
			 		          error("StarkStrongField(address,n1,n2,m) - params. are needed");
			 			 		        }  
		 else	{
			 
			 std::string addr(addressEG);

		 self->cpp_obj =  new  StarkStrongField(addr,n1,n2,m);
		 ((StarkStrongField*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
    
    
  }
  
  
  
  
  static PyObject* StarkStrongField_deltaE(PyObject *self, PyObject *args){
	  StarkStrongField* cpp_StarkStrongField = (StarkStrongField*)((pyORBIT_Object*) self)->cpp_obj;
  				       
  		
       int nVars = PyTuple_Size(args);
       double val;
       double mass;
       double E_x;
       double E_y;
       double E_z;
       double B_x;
       double B_y;
       double B_z;
       double px;
       double py;
       double pz;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"dddddddddd:",&mass,&E_x,&E_y,&E_z,&B_x,&B_y,&B_z,&px,&py,&pz))
             error(" getStarkStrongFieldEnergy(mass,E_x,E_y,E_z,B_x,B_y,B_z,px,py,pz) - parameters are needed");
           else
           val=cpp_StarkStrongField->getStarkStrongFielddeltaE(mass,E_x,E_y,E_z,B_x,B_y,B_z,px,py,pz);
           return Py_BuildValue("d",val);
  }
  
  
  	
  

  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void StarkStrongField_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((StarkStrongField*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef StarkStrongFieldClassMethods[] = {
		    { "deltaE",  StarkStrongField_deltaE,         		METH_VARARGS,"gets StarkStrongField Energy"},
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
	static PyMemberDef StarkStrongFieldClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_StarkStrongField_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"StarkStrongField", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) StarkStrongField_del , /*tp_dealloc*/
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
		"The StarkStrongField python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		StarkStrongFieldClassMethods, /* tp_methods */
		StarkStrongFieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) StarkStrongField_init, /* tp_init */
		0, /* tp_alloc */
		StarkStrongField_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initStarkStrongField(PyObject* module){
		if (PyType_Ready(&pyORBIT_StarkStrongField_Type) < 0) return;
		Py_INCREF(&pyORBIT_StarkStrongField_Type);
		PyModule_AddObject(module, "StarkStrongField", (PyObject *)&pyORBIT_StarkStrongField_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
