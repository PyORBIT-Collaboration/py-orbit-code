## \namespace orbit::utils
## \brief Utility classes.
##
## Classes:
## - multiDimDoubleArray - Method. Generates multi-dimensional array with doubles.
## - multiDimIntArray    - Method. Generates multi-dimensional array with integers.
## - orbitFinalize    - Method. Finalizes ORBIT script execution.
## - NamedObject      - Class. Represents an object with a name.
## - TypedObject      - Class. Represents an object with a type.
## - ParamsDictObject - Class. Represents an object that has a parameters dictionary.

from orbit.utils.multiDimArray    import multiDimDoubleArray
from orbit.utils.multiDimArray    import multiDimIntArray
from orbit.utils.orbitFinalize    import orbitFinalize
from orbit.utils.NamedObject      import NamedObject
from orbit.utils.TypedObject      import TypedObject
from orbit.utils.ParamsDictObject import ParamsDictObject

from orbit.utils.phaseOperations  import phaseNearTargetPhase, phaseNearTargetPhaseDeg

__all__ = []
__all__.append("multiDimDoubleArray")
__all__.append("multiDimIntArray")
__all__.append("orbitFinalize")
__all__.append("NamedObject")
__all__.append("TypedObject")
__all__.append("ParamsDictObject")
__all__.append("phaseNearTargetPhase")
__all__.append("phaseNearTargetPhaseDeg")

