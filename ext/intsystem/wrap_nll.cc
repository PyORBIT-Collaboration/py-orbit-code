#include "wrap_nll.hh"
#include "wrap_bunch.hh"
#include "wrap_syncpart.hh"

#include "pyORBIT_Object.hh"

#include "Bunch.hh"
#include "Nll.hh"
#include "ParticleAttributesFactory.hh"

namespace wrap_utils_nll{
	
	void error(const char* msg){ ORBIT_MPI_Finalize(msg); }
	
static PyObject* Nll_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

static int Nll_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		//instantiation of a new c++ Nll
		self->cpp_obj = (void*) new Nll();
		((Nll*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		//This is the way to create new class instance from the C-level
		// Template: PyObject* PyObject_CallMethod(	PyObject *o, char *method, char *format, ...)
		//see Python/C API documentation
		PyObject* mod = PyImport_ImportModule("nll");
    return 0;		
  }

static PyObject* nll_TRACK_EXT(PyObject *self, PyObject *args){
	PyObject *pyIn;
	if(!PyArg_ParseTuple(args,"O:TRACK_EXT",&pyIn)){
		error("PyMatrix - TRACK_EXT(Bunch) - Bunch is needed.");
	}			
	PyObject* pyBunchType = wrap_orbit_bunch::getBunchType("Bunch");
	if((!PyObject_IsInstance(pyIn,pyBunchType))){
		error("PyMatrix - TRACK_EXT(Bunch) - input parameter is not Bunch");
	}		
	Bunch* bunch = (Bunch*) ((pyORBIT_Object*) pyIn)->cpp_obj;
	Nll* cpp_nll = (Nll*) ((pyORBIT_Object *) self)->cpp_obj;
  cpp_nll->TRACK_EXT(bunch);
  Py_INCREF(Py_None);
  return Py_None;
}
	
	static void Nll_del(pyORBIT_Object* self){
		Nll* cpp_nll = (Nll*) self->cpp_obj;
		delete cpp_nll;
		self->ob_type->tp_free((PyObject*)self);
	}	

static PyMethodDef NllClassMethods[] = {
  {"TRACK_EXT", nll_TRACK_EXT, METH_VARARGS, "Track particle through non-linear lattice"},
  {NULL,NULL}
};

static PyMemberDef NllClassMembers [] = {
	{NULL}
};
	
	static PyTypeObject pyORBIT_Nll_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"Nll", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) Nll_del , /*tp_dealloc*/
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
		"The Nll python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		NllClassMethods, /* tp_methods */
		NllClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) Nll_init, /* tp_init */
		0, /* tp_alloc */
		Nll_new, /* tp_new */
	};	

static PyMethodDef NllModuleMethods[] = { 
	{NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  /*void initnll(){
		//check that the Nll wrapper is ready
		if (PyType_Ready(&pyORBIT_nll_Type) < 0) return;
		Py_INCREF(&pyORBIT_nll_Type);
    //create new module
    PyObject* module = Py_InitModule("nll",NLLModuleMethods);
		PyModule_AddObject(module, "nll", (PyObject *)&pyORBIT_nll_Type);			
		//add the SyncParticle python class
		//wrap_orbit_syncpart::initsyncpart(module);
  }*/
	
	void initNll(PyObject* module){
		if (PyType_Ready(&pyORBIT_Nll_Type) < 0) return;
		Py_INCREF(&pyORBIT_Nll_Type);
		PyModule_AddObject(module, "Nll", (PyObject *)&pyORBIT_Nll_Type);
	}

  	PyObject* getNLLType(const char* name){
		PyObject* mod = PyImport_ImportModule("nll");
		PyObject* pyType = PyObject_GetAttrString(mod,name);
		Py_DECREF(mod);
		Py_DECREF(pyType);
		return pyType;
		}	
	

#ifdef __cplusplus
}
#endif
}
