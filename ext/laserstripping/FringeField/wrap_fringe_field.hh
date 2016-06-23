#ifndef WRAP_FRINGE_FIELD_HH_
#define WRAP_FRINGE_FIELD_HH_

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_fringe_field{
    void initFringeField(PyObject* module);
  }

#ifdef __cplusplus
}
#endif


#endif /*WRAP_FRINGE_FIELD_HH_*/

