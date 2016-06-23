#ifndef WRAP_FUNCTIONS_HH_
#define WRAP_FUNCTIONS_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_functions{
    void initFunctions(PyObject* module);
  }

#ifdef __cplusplus
}
#endif


#endif /*WRAP_FUNCTIONS_HH_*/
