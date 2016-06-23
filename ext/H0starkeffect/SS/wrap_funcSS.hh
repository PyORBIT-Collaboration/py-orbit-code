#ifndef WRAP_FUNCSS_HH_
#define WRAP_FUNCSS_HH_


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_funcSS{
    void initFuncSS(PyObject* module);
  }

#ifdef __cplusplus
}
#endif


#endif /*WRAP_FUNCSS_HH_*/
