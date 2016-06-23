#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_wave_function.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "WaveFunction.hh"

//using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_wave_function{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* WaveFunction_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int WaveFunction_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  int n1,n2,m;
      int point1; 
      int pointN;
      int precEZ;
      const char* c_energy;
      const char* c_Z1;
      const char* c_F;



		 if(!PyArg_ParseTuple(	args,"iiiiiisss:",&n1, &n2, &m, &point1, &pointN,&precEZ, &c_energy, &c_Z1, &c_F)){
			 		          error("WaveFunction(n1, n2, m, point1, pointN,precEZ, c_energy, c_Z1, c_F) - params. are needed");
			 			 		        }  
		 else	{
			 

		 self->cpp_obj =  new  WaveFunction(n1,n2,m,point1, pointN, precEZ, c_energy, c_Z1, c_F);
		 ((WaveFunction*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
    
    
  }
  
  
  
 
  static PyObject* WaveFunction_M(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* mu;
	  std::string  M;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"s:",&mu))
             error(" getM(mu) - parameters are needed");
           else
        	   M = cpp_WaveFunction->getFastM(mu);

           return Py_BuildValue("s",M.c_str());
  }
  
  
  static PyObject* WaveFunction_N(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* mu;
	  std::string N;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"s:",&mu))
             error(" Nb(mu) - parameters are needed");
           else
        	   N = cpp_WaveFunction->getFastN(mu);

           return Py_BuildValue("s",N.c_str());
  }
  
  
  static PyObject* WaveFunction_getMN(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
  				       

	  double psi,mu,nu;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"dd:",&mu,&nu))
             error(" getMN(mu.nu) - parameters are needed");
           else
        	   psi = cpp_WaveFunction->getMN(mu,nu);

           return Py_BuildValue("d",psi);
  }
  
  
  static PyObject* WaveFunction_getMode(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* mu;
	  int mode;

       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,""))
             error(" getMode() -no parameters are needed");
           else
        	   
        	   mode = cpp_WaveFunction->getMode();

           return Py_BuildValue("i",mode);
  }
  
  /* 
  
  static PyObject* WaveFunction_calcPrecisionForN(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z1;
	  const char* c_F;
	  int prec;
	  std::string N = "fff";
	  std::string derN = "hhh";


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z1))
       				error(" getN(c_F, c_energy, c_Z1) - parameters are needed");
       		else	{
       			prec = cpp_WaveFunction->calcPrecisionForN(N, derN, c_F, c_energy, c_Z1);
 
       		}
       		
                
       		return Py_BuildValue("iss",prec, N.c_str(), derN.c_str());

    }	
  
  

  static PyObject* WaveFunction_get_a(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z2;
	  const char* c_F;
	  std::string a;
	  std::string der_a;


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z2))
       				error(" get_a(c_F, c_energy, c_Z2) - parameters are needed");
       		else	
       			cpp_WaveFunction->get_a(a, der_a, c_F, c_energy, c_Z2);
                
       		return Py_BuildValue("ss", a.c_str(), der_a.c_str());

    }	
  
  
  static PyObject* WaveFunction_get_b(PyObject *self, PyObject *args){
	  WaveFunction* cpp_WaveFunction = (WaveFunction*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z2;
	  const char* c_F;
	  std::string b;
	  std::string der_b;


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z2))
       				error(" get_b(c_F, c_energy, c_Z2) - parameters are needed");
       		else	
       			cpp_WaveFunction->get_b(b, der_b, c_F, c_energy, c_Z2);
                
       		return Py_BuildValue("ss", b.c_str(), der_b.c_str());

    }	
  
  

  
  */
  
	
  
  



  	
  

  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void WaveFunction_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((WaveFunction*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef WaveFunctionClassMethods[] = {
		    { "M",  					WaveFunction_M,         								METH_VARARGS,"gets M-function "},
			{ "N",      				WaveFunction_N,				       					 	METH_VARARGS,"Calculates N-function below point N1"},
			{ "getMN",     				WaveFunction_getMN,         							METH_VARARGS,"gets MN-perturbation functions"},
//			{ "get_a",         			WaveFunction_get_a,         							METH_VARARGS,"Gets value and derivative of a - function."},
//			{ "get_b",         			WaveFunction_get_b,         							METH_VARARGS,"Gets value and derivative of b - function."},
			{ "getMode",         		WaveFunction_getMode,         							METH_VARARGS,"mode of wave function ."},
/*			{ "getLaser_lambda",         HermiteGaussianLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  HermiteGaussianLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef WaveFunctionClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_WaveFunction_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"WaveFunction", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) WaveFunction_del , /*tp_dealloc*/
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
		"The WaveFunction python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		WaveFunctionClassMethods, /* tp_methods */
		WaveFunctionClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) WaveFunction_init, /* tp_init */
		0, /* tp_alloc */
		WaveFunction_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initWaveFunction(PyObject* module){
		if (PyType_Ready(&pyORBIT_WaveFunction_Type) < 0) return;
		Py_INCREF(&pyORBIT_WaveFunction_Type);
		PyModule_AddObject(module, "WaveFunction", (PyObject *)&pyORBIT_WaveFunction_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
