#include "orbit_mpi.hh"

//#include "wrap_density_matrix.hh"

#include "wrap_nll.hh"

static PyMethodDef intSystemMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

  void initintsystem(){
    //create new module
    PyObject* module = Py_InitModule("intsystem",intSystemMethods);
	  wrap_utils_nll::initNll(module);

  }
	
#ifdef __cplusplus
}
#endif
