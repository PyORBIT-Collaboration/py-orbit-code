#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_snsstrippingmagnet.hh"

#include <iostream>
#include <string>

#include "SNSstrippingMagnet.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_snsstrippingmagnet{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* SNSstrippingMagnet_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
	
	

  static int SNSstrippingMagnet_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  

	  
	  const char* file_name;	  


		 if(!PyArg_ParseTuple(	args,"s:",&file_name)){
		error("SNSstrippingMagnet(file_name ) - params. are needed");
			 			 		        }  
		 else	{
			 
			  std::string name_str_name(file_name);


		 self->cpp_obj =  new  SNSstrippingMagnet(name_str_name);
		 ((SNSstrippingMagnet*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

		
    return 0;
  }
  
  

  static PyObject* SNSstrippingMagnet_setA(PyObject *self, PyObject *args){
	  SNSstrippingMagnet* cpp_field = (SNSstrippingMagnet*)((pyORBIT_Object*) self)->cpp_obj;
	  double A;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&A))
             error(" setA(z) - parameters are needed");
           else
        	   cpp_field->setA(A);
           
  		    Py_INCREF(Py_None);
  		    return Py_None;	 
           
  }  
  

  
  
  static PyObject* SNSstrippingMagnet_setPosition(PyObject *self, PyObject *args){
	  SNSstrippingMagnet* cpp_field = (SNSstrippingMagnet*)((pyORBIT_Object*) self)->cpp_obj;
	  double pos;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&pos))
             error(" setPhase(z) - parameters are needed");
           else
        	   cpp_field->setPosition(pos);
 		    Py_INCREF(Py_None);
 		    return Py_None;	
           
  }   
    
  	// name([name]) - sets or returns the name of the External Effeects class 
  static PyObject* SNSstrippingMagnet_getBx(PyObject *self, PyObject *args){

		SNSstrippingMagnet* cpp_field = (SNSstrippingMagnet*)((pyORBIT_Object*) self)->cpp_obj;
                double z;

    if(!PyArg_ParseTuple(	args,"d",&z)){
      error("getBz(z) must have z value");
    } else
	return Py_BuildValue("d",cpp_field->getBx(z));
  }

  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void SNSstrippingMagnet_del(pyORBIT_Object* self){
		//std::cerr<<"The SNSstrippingMagnet __del__ has been called!"<<std::endl;
		delete ((SNSstrippingMagnet*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef SNSstrippingMagnetClassMethods[] = {
			{ "setA",         		SNSstrippingMagnet_setA,         		 METH_VARARGS,"Returns tuple of EM fields."},
			{ "getBx",         SNSstrippingMagnet_getBx,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setPosition",         SNSstrippingMagnet_setPosition,         METH_VARARGS,"Sets or returns the name of effects."},
/*			{ "getLaserPower",         HermiteGaussianLFmode_getLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "getLaser_lambda",         HermiteGaussianLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  HermiteGaussianLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef SNSstrippingMagnetClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_SNSstrippingMagnet_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SNSstrippingMagnet", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SNSstrippingMagnet_del , /*tp_dealloc*/
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
		"The SNSstrippingMagnet python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SNSstrippingMagnetClassMethods, /* tp_methods */
		SNSstrippingMagnetClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SNSstrippingMagnet_init, /* tp_init */
		0, /* tp_alloc */
		SNSstrippingMagnet_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initSNSstrippingMagnet(PyObject* module){
		if (PyType_Ready(&pyORBIT_SNSstrippingMagnet_Type) < 0) return;
		Py_INCREF(&pyORBIT_SNSstrippingMagnet_Type);
		PyModule_AddObject(module, "SNSstrippingMagnet", (PyObject *)&pyORBIT_SNSstrippingMagnet_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
