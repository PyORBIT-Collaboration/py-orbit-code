#include "orbit_mpi.hh"
#include "wrap_mpi_comm.hh"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <set>
#include <map>
#include <vector>

#include "ParticleAttributes.hh"
#include "SyncPart.hh"

//from utils
#include "AttributesBucket.hh"
#include "CppPyWrapper.hh"

using namespace std;

#ifndef NLL_H
#define NLL_H

class Nll: public OrbitUtils::CppPyWrapper{
public:
  Nll();
  void TRACK_EXT(Bunch* bunch);
};

#endif
