#ifndef WRAP_FUNCT_HH_
#define WRAP_FUNCT_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_funct{
    void initFunct(PyObject* module);
  }

#ifdef __cplusplus
}
#endif


#endif /*WRAP_FUNCTI_HH_*/
