#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"
#include "wrap_regular_grid_fs.hh"

#include <iostream>
#include <string>

#include "RegularGridFS.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_regular_grid_fs{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python CppBaseFieldSource class definition
	//---------------------------------------------------------

	//constructor for python class wrapping CppBaseFieldSource instance
	//It never will be called directly
	static PyObject* RegularGridFS_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}

  //initializator for python  CppBaseFieldSource class
  //this is implementation of the __init__ method
	
	

  static int RegularGridFS_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
	  

	  
	  const char* file_name;
	  const char* type_of_field;
	  double field_mult;
	  double coord_mult;
	  


		 if(!PyArg_ParseTuple(	args,"ssdd:",&file_name, &type_of_field, &coord_mult, &field_mult)){
			 		          error("HydrogenStarkParam(file_name , type_of_field, coodr_mult, field_mult) - params. are needed");
			 			 		        }  
		 else	{
			 
			  std::string name_str_name(file_name);
			  std::string name_str_field(type_of_field);

		 self->cpp_obj =  new  RegularGridFS(name_str_name, name_str_field, coord_mult, field_mult);
		 ((RegularGridFS*) self->cpp_obj)->setPyWrapper((PyObject*) self);
		 }
	

		
    return 0;
  }
  
  
  
  
  static PyObject* RegularGridFS_setFieldOrientation(PyObject *self, PyObject *args){
	  RegularGridFS* RegGrid = (RegularGridFS*)((pyORBIT_Object*) self)->cpp_obj;
 	 
   		double x0;
   		double y0;
   		double z0;
   		double kx;
   		double ky;
   		double kz;
   		double mx;
   		double my;
   		double mz;

   		double flag=true;

   		        if(!PyArg_ParseTuple(	args,"ddddddddd:",&x0, &y0, &z0, &kx, &ky, &kz, &mx, &my, &mz))
   		          {error("LaserExternalEfects - setLaserHalfAngle(x0, y0, z0, kx, ky, kz, mx, my, mz) - params. are needed");flag=false;}
  		        
   		        if(fabs(kx*mx+ky*my+kz*mz)/sqrt(kx*kx+ky*ky+kz*kz)/sqrt(mx*mx+my*my+mz*mz)>1.e-15)
   		          {error("Please be shure that kx*mx+ky*my+kz*mz==0");flag=false;}
		        
   		        if(flag)	
   		        	RegGrid->setFieldOrientation(x0,y0, z0, kx, ky, kz, mx, my, mz);
   		      
   	   
   		    Py_INCREF(Py_None);
   		    return Py_None;	  
    }	
  

  
  
  static PyObject* RegularGridFS_getFields(PyObject *self, PyObject *args){
	  RegularGridFS* cpp_fields = (RegularGridFS*)((pyORBIT_Object*) self)->cpp_obj;
  				       
  		
       int nVars = PyTuple_Size(args);
       
       double x;
       double y;
       double z;
       double t;
       double Ex;
       double Ey;
       double Ez;
       double Bx;
       double By;
       double Bz;
       

           //NO NEW OBJECT CREATED BY PyArg_ParseTuple! NO NEED OF Py_DECREF()
           if(!PyArg_ParseTuple(	args,"dddd:",&x,&y,&z,&t))
             error(" getFields(x,y,z,t) - parameters are needed");
           else
           cpp_fields->getElectricMagneticField(x,y,z,t,Ex,Ey,Ez,Bx,By,Bz);
           
           if (Ex==0.0&&Ey==0.0&&Ez==0.0) return Py_BuildValue("ddd",Bx,By,Bz);
           if (Bx==0.0&&By==0.0&&Bz==0.0) return Py_BuildValue("ddd",Ex,Ey,Ez);
  }
  
  


  //-----------------------------------------------------
  //destructor for python PyBaseFieldSource class (__del__ method).
  //-----------------------------------------------------
  static void RegularGridFS_del(pyORBIT_Object* self){
		//std::cerr<<"The RegularGridFS __del__ has been called!"<<std::endl;
		delete ((RegularGridFS*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }

	// defenition of the methods of the python PyBaseFieldSource wrapper class
	// they will be vailable from python level
  static PyMethodDef RegularGridFSClassMethods[] = {
		    { "setFieldOrientation",  RegularGridFS_setFieldOrientation,         METH_VARARGS,"Sets field orientation."},
			{ "getFields",         		RegularGridFS_getFields,         		 METH_VARARGS,"Returns tuple of EM fields."},
/*			{ "setLaserPower",         HermiteGaussianLFmode_setLaserPower,         METH_VARARGS,"Sets or returns the name of effects."},
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
	static PyMemberDef RegularGridFSClassMembers [] = {
		{NULL}
	};

	//new python PyBaseFieldSource wrapper type definition
	static PyTypeObject pyORBIT_RegularGridFS_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"RegularGridFS", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) RegularGridFS_del , /*tp_dealloc*/
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
		"The RegularGridFS python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		RegularGridFSClassMethods, /* tp_methods */
		RegularGridFSClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) RegularGridFS_init, /* tp_init */
		0, /* tp_alloc */
		RegularGridFS_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pyPyBaseFieldSource class
	//It will be called from Bunch wrapper initialization
	//--------------------------------------------------
  void initRegularGridFS(PyObject* module){
		if (PyType_Ready(&pyORBIT_RegularGridFS_Type) < 0) return;
		Py_INCREF(&pyORBIT_RegularGridFS_Type);
		PyModule_AddObject(module, "RegularGridFS", (PyObject *)&pyORBIT_RegularGridFS_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_utils_cpp_base_field_source
}
