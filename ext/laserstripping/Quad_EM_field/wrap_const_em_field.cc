#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_const_em_field.hh"

#include <iostream>
#include <string>

#include "ConstEMfield.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_const_em_field{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* ConstEMfield_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int ConstEMfield_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
//	  self->cpp_obj = new LSFieldSource();
	 
	  double E_x;
	  double E_y;
	  double E_z;
	  double B_x;
	  double B_y;
	  double B_z;
	  
	  int nVars = PyTuple_Size(args);
	     

//	 if (PyArg_ParseTuple(	args,"")) {self->cpp_obj = new LSFieldSource();}
//	 if(PyArg_ParseTuple(args,"dddddd:",&E_x,&E_y,&E_z,&B_x,&B_y,&B_z))	{ self->cpp_obj = new LSFieldSource(E_x,E_y,E_z,B_x,B_y,B_z);}
	  
	  if(nVars==0)	 if (!PyArg_ParseTuple(	args,"")) {} else {self->cpp_obj = new ConstEMfield();}
	  if(nVars==6)	 if (!PyArg_ParseTuple(args,"dddddd:",&E_x,&E_y,&E_z,&B_x,&B_y,&B_z)) 
	  {error("Parameters  (E_x,E_y,E_z,B_x,B_y,B_z) -are needed");} else {self->cpp_obj = new ConstEMfield(E_x,E_y,E_z,B_x,B_y,B_z);}

	     ((ConstEMfield*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		
    return 0;
  }
  
  
  static PyObject* ConstEMfield_Ex(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ex;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ex))
             error(" SetupPrint(Ex - parameter is needed");
           else 	  
           cpp_ConstEMfield->Ex(Ex);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* ConstEMfield_Ey(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ey;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ey))
             error(" SetupPrint(Ey - parameter is needed");
           else 	  
           cpp_ConstEMfield->Ey(Ey);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* ConstEMfield_Ez(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ez;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ez))
             error(" SetupPrint(Ez - parameter is needed");
           else 	  
           cpp_ConstEMfield->Ez(Ez);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* ConstEMfield_Bx(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Bx;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Bx))
             error(" SetupPrint(Bx - parameter is needed");
           else 	  
           cpp_ConstEMfield->Bx(Bx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* ConstEMfield_By(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double By;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&By))
             error(" SetupPrint(By - parameter is needed");
           else 	  
           cpp_ConstEMfield->By(By);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* ConstEMfield_Bz(PyObject *self, PyObject *args){
	  ConstEMfield* cpp_ConstEMfield = (ConstEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Bz;
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Bz))
             error(" SetupPrint(Bz - parameter is needed");
           else 	  
           cpp_ConstEMfield->Bz(Bz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void ConstEMfield_del(pyORBIT_Object* self){
		//std::cerr<<"The ConstEMfield __del__ has been called!"<<std::endl;
		delete ((ConstEMfield*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef ConstEMfieldClassMethods[] = {
		{ "Ex",        			 ConstEMfield_Ex,        		METH_VARARGS,"Sets parameter Ex of EM field."},
		{ "Ey",        			 ConstEMfield_Ey,        		METH_VARARGS,"Sets parameter Ey of EM field."},
		{ "Ez",        			 ConstEMfield_Ez,        		METH_VARARGS,"Sets parameter Ez of EM field."},
		{ "Bx",        			 ConstEMfield_Bx,        		METH_VARARGS,"Sets parameter Bx of EM field."},
		{ "By",        			 ConstEMfield_By,        		METH_VARARGS,"Sets parameter By of EM field."},
		{ "Bz",        			 ConstEMfield_Bz,        		METH_VARARGS,"Sets parameter Bz of EM field."},
		
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef ConstEMfieldClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_ConstEMfield_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"ConstEMfield", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) ConstEMfield_del , /*tp_dealloc*/
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
		"The ConstEMfield python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		ConstEMfieldClassMethods, /* tp_methods */
		ConstEMfieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) ConstEMfield_init, /* tp_init */
		0, /* tp_alloc */
		ConstEMfield_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initConstEMfield(PyObject* module){
		if (PyType_Ready(&pyORBIT_ConstEMfield_Type) < 0) return;
		Py_INCREF(&pyORBIT_ConstEMfield_Type);
		PyModule_AddObject(module, "ConstEMfield", (PyObject *)&pyORBIT_ConstEMfield_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
