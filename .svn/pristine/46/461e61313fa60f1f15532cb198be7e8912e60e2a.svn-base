#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_stark.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "Stark.hh"

//using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_stark{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* Stark_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int Stark_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  int states;
	  const char* addressEG;

		 if(!PyArg_ParseTuple(	args,"si:",&addressEG,&states)){
			 		          error("Stark(address,states) - params. are needed");
			 			 		        }  
		 else	{
			 
			 std::string addr(addressEG);

		 self->cpp_obj =  new  Stark(addr,states);
		 ((Stark*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
    
    
  }
  
  
  
  
  static PyObject* Stark_getStarkEnergy(PyObject *self, PyObject *args){
	  Stark* cpp_Stark = (Stark*)((pyORBIT_Object*) self)->cpp_obj;
  				       
  		
       int nVars = PyTuple_Size(args);
       double val;
       double mass;
       int n1;
       int n2;
       int m;
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
           if(!PyArg_ParseTuple(	args,"diiiddddddddd:",&mass,&n1,&n2,&m,&E_x,&E_y,&E_z,&B_x,&B_y,&B_z,&px,&py,&pz))
             error(" getStarkEnergy(mass,n1,n2,m,E_x,E-y,E-z,B_x,B_y,B_z,px,py,pz) - parameters are needed");
           else
           val=cpp_Stark->getStarkEnergy(mass,n1,n2,m,E_x,E_y,E_z,B_x,B_y,B_z,px,py,pz);
           return Py_BuildValue("d",val);
  }
  
  
  	
  

  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void Stark_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((Stark*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef StarkClassMethods[] = {
		    { "getStarkEnergy",  Stark_getStarkEnergy,         		METH_VARARGS,"gets Stark Energy"},
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
	static PyMemberDef StarkClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_Stark_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Stark", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Stark_del , /*tp_dealloc*/
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
		"The Stark python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		StarkClassMethods, /* tp_methods */
		StarkClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Stark_init, /* tp_init */
		0, /* tp_alloc */
		Stark_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initStark(PyObject* module){
		if (PyType_Ready(&pyORBIT_Stark_Type) < 0) return;
		Py_INCREF(&pyORBIT_Stark_Type);
		PyModule_AddObject(module, "Stark", (PyObject *)&pyORBIT_Stark_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
