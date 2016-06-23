#ifndef WRAP_SCHRODINGER_EQUATION_HH_
#define WRAP_SCHRODINGER_EQUATION_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_schrodinger_equation{
    void initSchrodingerEquation(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_SCHRODINGER_EQUATION_HH_*/



