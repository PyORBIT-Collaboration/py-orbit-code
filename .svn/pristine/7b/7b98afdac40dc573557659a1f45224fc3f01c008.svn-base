#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_funcSS.hh"

#include <iostream>
#include <string>
#include <cmath>


#include "FuncSS.hh"

//using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_FuncSS{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* FuncSS_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
  static int FuncSS_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){

	  int n1,n2,m;
      int point1; 



      
		 if(!PyArg_ParseTuple(	args,"iiii:",&n1, &n2, &m, &point1)){
			 		          error("FuncSS(n1, n2, m, point1,err_exp, exp_minG) - params. are needed");
			 			 		        }  
		 else	{
			 

		 self->cpp_obj =  new  FuncSS(n1,n2,m, point1);
		 ((FuncSS*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

    return 0;
    
    
  }
  
  
  
  
  static PyObject* FuncSS_getM(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* c_energy;
	  const char* c_Z1;
	  const char* c_F;
	  std::string M;


       
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z1))
             error(" getM(k, step, nsum) - parameters are needed");
           else
        	   M = cpp_FuncSS->getM(c_F, c_energy, c_Z1);

           return Py_BuildValue("s",M.c_str());
  }
  
  
  
  
  static PyObject* FuncSS_M(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* c_energy;
	  const char* c_Z1;
	  const char* c_F;
	  const char* c_mu;
	  std::string M;


       
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"ssss:",&c_mu,&c_F,&c_energy,&c_Z1))
             error(" M(mu, k, step, nsum) - parameters are needed");
           else
        	   M = cpp_FuncSS->M(c_mu, c_F, c_energy, c_Z1);

           return Py_BuildValue("s",M.c_str());
  }
  
  
  static PyObject* FuncSS_N(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
  				       
	  const char* c_energy;
	  const char* c_Z2;
	  const char* c_F;
	  const char* c_nu;
	  std::string N;

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"ssss:",&c_nu,&c_F,&c_energy,&c_Z2))
             error(" M(mu, k, step, nsum) - parameters are needed");
           else
        	   N = cpp_FuncSS->N(c_nu, c_F, c_energy, c_Z2);

           return Py_BuildValue("s",N.c_str());
  }
  
  
  
  static PyObject* FuncSS_calcPrecisionForN(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z1;
	  const char* c_F;
	  int pointN;
	  std::string N = "fff";
	  std::string derN = "hhh";


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z1))
       				error(" getN(c_F, c_energy, c_Z1) - parameters are needed");
       		else	{
       			pointN = cpp_FuncSS->calcPrecisionForN(N, derN, c_F, c_energy, c_Z1);
 
       		}
       		
                
       		return Py_BuildValue("i",pointN);

    }	
  
  
 
  
  static PyObject* FuncSS_getB(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z2;
	  const char* c_F;

	  std::string B;


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z2))
       				error(" get_b(c_F, c_energy, c_Z2) - parameters are needed");
       		else	
       				B = cpp_FuncSS->getB(c_F, c_energy, c_Z2);
                
       		return Py_BuildValue("s", B.c_str());

    }	
  
  
  static PyObject* FuncSS_C(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;
 	 

	  const char* c_energy;
	  const char* c_Z2;
	  const char* c_F;

	  std::string C;


	  	//NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
       		if(!PyArg_ParseTuple(	args,"sss:",&c_F,&c_energy,&c_Z2))
       				error(" get_b(c_F, c_energy, c_Z2) - parameters are needed");
       		else	
       				C = cpp_FuncSS->C(c_F, c_energy, c_Z2);
                
       		return Py_BuildValue("s", C.c_str());

    }	
  
  
  

	
  
  
  
  static PyObject* FuncSS_calcPrecisionForM(PyObject *self, PyObject *args){
	  FuncSS* cpp_FuncSS = (FuncSS*)((pyORBIT_Object*) self)->cpp_obj;			       
	  const char* c_field;
	  const char* c_energy;
	  const char* c_Z1;
	  int val;
	  int max_err_exp;

	  
           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"isss:", &max_err_exp, &c_field,&c_energy,&c_Z1))
             error(" calcPrecisionForM(max_err_exp,c_field,c_energy,c_Z1) - parameter are needed");
           else 	  
        	   val = cpp_FuncSS->calcPrecisionForM(max_err_exp,c_field,c_energy,c_Z1);
         
  		  return Py_BuildValue("i",val);
  		    
  		    
  }
  
  


  	
  

  

  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void FuncSS_del(pyORBIT_Object* self){
		//std::cerr<<"The CppBaseFieldSource __del__ has been called!"<<std::endl;
		delete ((FuncSS*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef FuncSSClassMethods[] = {
		    { "getM",  					FuncSS_getM,         							METH_VARARGS,"gets M-FuncSSion Modulus"},
			{ "calcPrecisionForM",      FuncSS_calcPrecisionForM,				        METH_VARARGS,"Sets the parameter of electric field."},
			{ "calcPrecisionForN",      FuncSS_calcPrecisionForN,         				METH_VARARGS,"gets N-FuncSSion Modulus"},
			{ "getB",         			FuncSS_getB,         							METH_VARARGS,"Gets amplitude of incoming wave ."},
			{ "M",         				FuncSS_M,         								METH_VARARGS,"gets M-FuncSS ."},
			{ "N",         				FuncSS_N,         								METH_VARARGS,"gets N-FuncSS ."},
			{ "C",         				FuncSS_C,         								METH_VARARGS,"gets C-norm. coeff ."},
/*			{ "getLaser_lambda",         HermiteGaussianLFmode_getLaser_lambda,         METH_VARARGS,"Sets or returns the name of effects."},
			{ "setLaserFieldOrientation",  HermiteGaussianLFmode_setLaserFieldOrientation,         METH_VARARGS,"Sets or returns the name of effects."},
*/
    {NULL}
  };

	// defenition of the memebers of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
	static PyMemberDef FuncSSClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_FuncSS_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"FuncSS", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) FuncSS_del , /*tp_dealloc*/
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
		"The FuncSS python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		FuncSSClassMethods, /* tp_methods */
		FuncSSClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) FuncSS_init, /* tp_init */
		0, /* tp_alloc */
		FuncSS_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization FuncSSion of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initFuncSS(PyObject* module){
		if (PyType_Ready(&pyORBIT_FuncSS_Type) < 0) return;
		Py_INCREF(&pyORBIT_FuncSS_Type);
		PyModule_AddObject(module, "FuncSS", (PyObject *)&pyORBIT_FuncSS_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
