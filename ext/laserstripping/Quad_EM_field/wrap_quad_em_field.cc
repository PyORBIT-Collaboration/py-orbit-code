#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_quad_em_field.hh"

#include <iostream>
#include <string>

#include "QuadEMfield.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_quad_em_field{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* QuadEMfield_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int QuadEMfield_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  
	  self->cpp_obj = new QuadEMfield();
	  
	  ((QuadEMfield*) self->cpp_obj)->setPyWrapper((PyObject*) self);
	 
    return 0;
  }
  
  
  static PyObject* QuadEMfield_x0(PyObject *self, PyObject *args){
 	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
        double x0;
            //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
            if(!PyArg_ParseTuple(	args,"d:",&x0))
              error(" SetupPrint(x0 - parameter is needed");
            else 	  
            cpp_QuadEMfield->x0(x0);
          
   		    Py_INCREF(Py_None);
   		    return Py_None;
   }
   
   static PyObject* QuadEMfield_y0(PyObject *self, PyObject *args){
 	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
        double y0;
            //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
            if(!PyArg_ParseTuple(	args,"d:",&y0))
              error(" SetupPrint(y0 - parameter is needed");
            else 	  
            cpp_QuadEMfield->y0(y0);
          
   		    Py_INCREF(Py_None);
   		    return Py_None;
   }
   
   static PyObject* QuadEMfield_z0(PyObject *self, PyObject *args){
 	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
        double z0;
            //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
            if(!PyArg_ParseTuple(	args,"d:",&z0))
              error(" SetupPrint(z0 - parameter is needed");
            else 	  
            cpp_QuadEMfield->z0(z0);
          
   		    Py_INCREF(Py_None);
   		    return Py_None;
   }
   
  
  
  static PyObject* QuadEMfield_Ex0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ex0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ex0))
             error(" SetupPrint(Ex0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->Ex0(Ex0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_Ey0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ey0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ey0))
             error(" SetupPrint(Ey0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->Ey0(Ey0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_Ez0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Ez0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Ez0))
             error(" SetupPrint(Ez0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->Ez0(Ez0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_Bx0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Bx0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Bx0))
             error(" SetupPrint(Bx0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->Bx0(Bx0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_By0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double By0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&By0))
             error(" SetupPrint(By0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->By0(By0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_Bz0(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double Bz0;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&Bz0))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->Bz0(Bz0);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxEx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxEx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxEx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxEx(cxEx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyEx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyEx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyEx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyEx(cyEx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czEx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czEx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czEx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czEx(czEx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxEy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxEy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxEy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxEy(cxEy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyEy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyEy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyEy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyEy(cyEy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czEy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czEy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czEy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czEy(czEy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxEz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxEz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxEz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxEz(cxEz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyEz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyEz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyEz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyEz(cyEz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czEz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czEz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czEz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czEz(czEz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxBx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxBx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxBx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxBx(cxBx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyBx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyBx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyBx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyBx(cyBx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czBx(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czBx;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czBx))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czBx(czBx);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxBy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxBy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxBy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxBy(cxBy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyBy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyBy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyBy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyBy(cyBy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czBy(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czBy;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czBy))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czBy(czBy);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cxBz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cxBz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cxBz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cxBz(cxBz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_cyBz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double cyBz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&cyBz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->cyBz(cyBz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  static PyObject* QuadEMfield_czBz(PyObject *self, PyObject *args){
	  QuadEMfield* cpp_QuadEMfield = (QuadEMfield*)((pyORBIT_Object*) self)->cpp_obj;			       
       double czBz;
           //NO NEW OBJECT CREATED By0 PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"d:",&czBz))
             error(" SetupPrint(Bz0 - parameter is needed");
           else 	  
           cpp_QuadEMfield->czBz(czBz);
         
  		    Py_INCREF(Py_None);
  		    return Py_None;
  }
  
  
  
  
  
  
  
  
  
  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void QuadEMfield_del(pyORBIT_Object* self){
		//std::cerr<<"The LSFieldSource __del__ has been called!"<<std::endl;
		delete ((QuadEMfield*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef QuadEMfieldClassMethods[] = {
			{ "x0",        				 QuadEMfield_x0,        		METH_VARARGS,"Sets parameter x0 of EM field."},
			{ "y0",        			 	 QuadEMfield_y0,        		METH_VARARGS,"Sets parameter y0 of EM field."},
			{ "z0",        			 	 QuadEMfield_z0,        		METH_VARARGS,"Sets parameter z0 of EM field."},
			{ "Ex0",        			 QuadEMfield_Ex0,        		METH_VARARGS,"Sets parameter Ex0 of EM field."},
			{ "Ey0",        			 QuadEMfield_Ey0,        		METH_VARARGS,"Sets parameter Ey0 of EM field."},
			{ "Ez0",        			 QuadEMfield_Ez0,        		METH_VARARGS,"Sets parameter Ez0 of EM field."},
			{ "Bx0",        			 QuadEMfield_Bx0,        		METH_VARARGS,"Sets parameter Bx0 of EM field."},
			{ "By0",        			 QuadEMfield_By0,        		METH_VARARGS,"Sets parameter By0 of EM field."},
			{ "Bz0",        			 QuadEMfield_Bz0,        		METH_VARARGS,"Sets parameter Bz0 of EM field."},
			{ "cxEx",        			 QuadEMfield_cxEx,        		METH_VARARGS,"Sets parameter cxEx of EM field."},
			{ "cyEx",        			 QuadEMfield_cyEx,        		METH_VARARGS,"Sets parameter cyEx of EM field."},
			{ "czEx",        			 QuadEMfield_czEx,        		METH_VARARGS,"Sets parameter czEx of EM field."},
			{ "cxEy",        			 QuadEMfield_cxEy,        		METH_VARARGS,"Sets parameter cxEy of EM field."},
			{ "cyEy",        			 QuadEMfield_cyEy,        		METH_VARARGS,"Sets parameter cyEy of EM field."},
			{ "czEy",        			 QuadEMfield_czEy,        		METH_VARARGS,"Sets parameter czEy of EM field."},
			{ "cxEz",        			 QuadEMfield_cxEz,        		METH_VARARGS,"Sets parameter cxEz of EM field."},
			{ "cyEz",        			 QuadEMfield_cyEz,        		METH_VARARGS,"Sets parameter cyEz of EM field."},
			{ "czEz",        			 QuadEMfield_czEz,        		METH_VARARGS,"Sets parameter czEz of EM field."},
			{ "cxBx",        			 QuadEMfield_cxBx,        		METH_VARARGS,"Sets parameter cxBx of EM field."},
			{ "cyBx",        			 QuadEMfield_cyBx,        		METH_VARARGS,"Sets parameter cyBx of EM field."},
			{ "czBx",        			 QuadEMfield_czBx,        		METH_VARARGS,"Sets parameter czBx of EM field."},
			{ "cxBy",        			 QuadEMfield_cxBy,        		METH_VARARGS,"Sets parameter cxBy of EM field."},
			{ "cyBy",        			 QuadEMfield_cyBy,        		METH_VARARGS,"Sets parameter cyBy of EM field."},
			{ "czBy",        			 QuadEMfield_czBy,        		METH_VARARGS,"Sets parameter czBy of EM field."},
			{ "cxBz",        			 QuadEMfield_cxBz,        		METH_VARARGS,"Sets parameter cxBz of EM field."},
			{ "cyBz",        			 QuadEMfield_cyBz,        		METH_VARARGS,"Sets parameter cyBz of EM field."},
			{ "czBz",        			 QuadEMfield_czBz,        		METH_VARARGS,"Sets parameter czBz of EM field."},
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef QuadEMfieldClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_QuadEMfield_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"QuadEMfield", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) QuadEMfield_del , /*tp_dealloc*/
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
		"The QuadEMfield python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		QuadEMfieldClassMethods, /* tp_methods */
		QuadEMfieldClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) QuadEMfield_init, /* tp_init */
		0, /* tp_alloc */
		QuadEMfield_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initQuadEMfield(PyObject* module){
		if (PyType_Ready(&pyORBIT_QuadEMfield_Type) < 0) return;
		Py_INCREF(&pyORBIT_QuadEMfield_Type);
		PyModule_AddObject(module, "QuadEMfield", (PyObject *)&pyORBIT_QuadEMfield_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
