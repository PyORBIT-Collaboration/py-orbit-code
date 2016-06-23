#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_DM_noLaserField.hh"

#include <iostream>
#include <string>

#include "DM_noLaserField.hh"
#include "BaseLaserFieldSource.hh"
#include "HydrogenStarkParam.hh"

using namespace OrbitUtils;
using namespace LaserStripping;

namespace wrap_DM_noLaserField{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppExternalEffects class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppExternalEffects instance
	//It never will be called directly
	static PyObject* DM_noLaserField_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  PyExternalEffects class
  //this is implementation of the __init__ method
  static int DM_noLaserField_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  PyObject*	pyStarkEffect=NULL;
		
	  int nVars = PyTuple_Size(args);
          int method = 0;
	  
	  if(!PyArg_ParseTuple(	args,"Oi:",&pyStarkEffect, &method))
		{
			error("DM_noLaserField(StarkEffect) - param.  needed");
		}	
	  
	  else	{
			
			Stark* Starkef = (Stark*) ((pyORBIT_Object*) pyStarkEffect)->cpp_obj;
			self->cpp_obj =  new  DM_noLaserField(Starkef, method);
			((DM_noLaserField*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		}
		
    return 0;
  }
  
	
  
	// name([name]) - sets or returns the name of the External Effeects class 
  static PyObject* DM_noLaserField_name(PyObject *self, PyObject *args){
		DM_noLaserField* cpp_DM_noLaserField = (DM_noLaserField*) ((pyORBIT_Object*) self)->cpp_obj;
    const char* name = NULL;
    if(!PyArg_ParseTuple(	args,"|s:name",&name)){
      error("DM_noLaserField - call should be - name([name]).");
    }
		if(name != NULL){
      std::string name_str(name);
      cpp_DM_noLaserField->setName(name_str);
		}
		return Py_BuildValue("s",cpp_DM_noLaserField->getName().c_str());
  }	
  
  //-----------------------------------------------------
  //destructor for python PyExternalEffects class (__del__ method).
  //-----------------------------------------------------
  static void DM_noLaserField_del(pyORBIT_Object* self){
		//std::cerr<<"The DM_noLaserField __del__ has been called!"<<std::endl;
		delete ((DM_noLaserField*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyExternalEffects wrapper class
	// they will be vailable from python level
  static PyMethodDef DM_noLaserFieldClassMethods[] = {
		{ "name",        			 DM_noLaserField_name,        		METH_VARARGS,"Sets or returns the name of effects."},
    {NULL}
  };

	// defenition of the memebers of the python PyExternalEffects wrapper class
	// they will be vailable from python level
	static PyMemberDef DM_noLaserFieldClassMembers [] = {
		{NULL}
	};

	//new python PyExternalEffects wrapper type definition
	static PyTypeObject pyORBIT_DM_noLaserField_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"DM_noLaserField", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) DM_noLaserField_del , /*tp_dealloc*/
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
		"The DM_noLaserField python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		DM_noLaserFieldClassMethods, /* tp_methods */
		DM_noLaserFieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) DM_noLaserField_init, /* tp_init */
		0, /* tp_alloc */
		DM_noLaserField_new, /* tp_new */
	};	


		
	//--------------------------------------------------
	//Initialization function of the pyPyExternalEffects class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initDM_noLaserField(PyObject* module){
		if (PyType_Ready(&pyORBIT_DM_noLaserField_Type) < 0) return;
		Py_INCREF(&pyORBIT_DM_noLaserField_Type);
		PyModule_AddObject(module, "DM_noLaserField", (PyObject *)&pyORBIT_DM_noLaserField_Type);
				
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_tracker3dfield
}
