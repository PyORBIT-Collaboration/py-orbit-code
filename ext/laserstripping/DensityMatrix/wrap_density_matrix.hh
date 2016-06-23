#ifndef WRAP_DENSITY_MATRIX_HH_
#define WRAP_DENSITY_MATRIX_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_density_matrix{
    void initDensityMatrix(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif /*WRAP_DENSITY_MATRIX_HH_*/



