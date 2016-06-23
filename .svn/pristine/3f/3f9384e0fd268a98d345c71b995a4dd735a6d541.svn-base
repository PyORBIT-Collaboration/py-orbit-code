#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_continuum_ss.hh"

#include <iostream>
#include <string>

#include "ContinuumSS.hh"
#include "BaseLaserFieldSource.hh"
#include "TDMcontinuum.hh"


using namespace OrbitUtils;
using namespace LaserStripping;

namespace wrap_continuum_ss{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python ContinuumSS class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppExternalEffects instance
	//It never will be called directly
	static PyObject* ContinuumSS_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  PyExternalEffects class
  //this is implementation of the __init__ method
  static int ContinuumSS_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  



	  PyObject*	pyBaseLaserField=NULL;
	  PyObject*	pycontinuum=NULL;

		 if(!PyArg_ParseTuple(	args,"OO:",&pyBaseLaserField,&pycontinuum)){
			 		          error("ContinuumSS(LaserField,StarkEffect) - params. are needed");
			 			 		        }  
		 else	{
		 BaseLaserFieldSource* lfs = (BaseLaserFieldSource*) ((pyORBIT_Object*) pyBaseLaserField)->cpp_obj;
		 TDMcontinuum* con = (TDMcontinuum*) ((pyORBIT_Object*) pycontinuum)->cpp_obj;
		 self->cpp_obj =  new  ContinuumSS(lfs, con);
		 ((ContinuumSS*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
  }
  
		
  
	// name([name]) - sets or returns the name of the External Effeects class 
  static PyObject* ContinuumSS_name(PyObject *self, PyObject *args){
	  ContinuumSS* cpp_ContinuumSS = (ContinuumSS*) ((pyORBIT_Object*) self)->cpp_obj;
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      error("ContinuumSS - call should be - name([name]).");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_ContinuumSS->setName(name_str);
		}
		return Py_BuildValue("s",cpp_ContinuumSS->getName().c_str());
  }	
  
  
  
  
  
	
  //-----------------------------------------------------
  //destructor for python PyExternalEffects class (__del__ method).
  //-----------------------------------------------------
  static void ContinuumSS_del(pyORBIT_Object* self){
		//std::cerr<<"The DensityMatrix __del__ has been called!"<<std::endl;
		delete ((ContinuumSS*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyExternalEffects wrapper class
	// they will be vailable from python level
  static PyMethodDef ContinuumSSClassMethods[] = {
		{ "name",        			 ContinuumSS_name,        		METH_VARARGS,"Sets or returns the name of effects."},


    {NULL}
  };

	// defenition of the memebers of the python PyExternalEffects wrapper class
	// they will be vailable from python level
	static PyMemberDef ContinuumSSClassMembers [] = {
		{NULL}
	};

	//new python PyExternalEffects wrapper type definition
	static PyTypeObject pyORBIT_ContinuumSS_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ContinuumSS", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ContinuumSS_del , /*tp_dealloc*/
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
		"The ContinuumSS python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ContinuumSSClassMethods, /* tp_methods */
		ContinuumSSClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ContinuumSS_init, /* tp_init */
		0, /* tp_alloc */
		ContinuumSS_new, /* tp_new */
	};	


		
	//--------------------------------------------------
	//Initialization function of the pyPyExternalEffects class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initContinuumSS(PyObject* module){
		if (PyType_Ready(&pyORBIT_ContinuumSS_Type) < 0) return;
		Py_INCREF(&pyORBIT_ContinuumSS_Type);
		PyModule_AddObject(module, "ContinuumSS", (PyObject *)&pyORBIT_ContinuumSS_Type);
				
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_tracker3dfield
}
