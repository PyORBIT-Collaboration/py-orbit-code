include ../../conf/make_root_config

DIRS  = $(patsubst %/, %, $(filter-out obj/,$(filter %/,$(shell ls -F))))
SRCS  = $(wildcard *.cc)
SRCS += $(foreach dir,$(DIRS),$(patsubst $(dir)/%.cc,%.cc,$(wildcard $(dir)/*.cc)))

OBJS = $(patsubst %.cc,./obj/%.o,$(SRCS))

#include files could be everywhere, we use only two levels
UPPER_DIRS = $(filter-out test%,$(patsubst %/, %,$(filter %/,$(shell ls -F ../../src))))
LOWER_DIRS = $(foreach dir,$(UPPER_DIRS),$(patsubst %/, ../../src/$(dir)/%,$(filter %/,$(shell ls -F ../../src/$(dir)))))

INCLUDES_LOCAL = $(patsubst %, -I../../src/%, $(UPPER_DIRS))
INCLUDES_LOCAL += $(filter-out %obj,$(patsubst %, -I%, $(LOWER_DIRS)))
INCLUDES_LOCAL += $(patsubst %, -I./%, $(filter %/,$(shell ls -F ./)))
INCLUDES_LOCAL += -I./

INC  = $(wildcard *.hh)
INC += $(wildcard *.h)
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.hh))
INC += $(foreach dir,$(DIRS),$(wildcard ./$(dir)/*.h))

#-------------------------------------------------------------------------------
# libraries files locations
#-------------------------------------------------------------------------------

LIBS +=  -L/home/tg4/tools/gmp-4.2.4/lib  -lgmp
LIBS +=  -L/home/tg4/tools/mpc-0.7/lib    -lmpc
LIBS +=  -L/home/tg4/tools/mpfr-2.4.1/lib -lmpfr
LIBS +=  -L/home/tg4/tools/gmpfrxx 		-lgmpfrxx



#-------------------------------------------------------------------------------
# include file locations
#-------------------------------------------------------------------------------


INCLUDES +=    -I/home/tg4/tools/mpc-0.7/include
INCLUDES +=    -I/home/tg4/tools/gmp-4.2.4/include
INCLUDES +=    -I/home/tg4/tools/mpfr-2.4.1/include
INCLUDES +=    -I/home/tg4/tools/gmpfrxx

#wrappers CC FLAGS
WRAPPER_FLAGS = -fno-strict-aliasing

#CXXFLAGS
CXXFLAGS += -fPIC

#shared library flags
SHARED_LIB = -shared

#tracker shared library
seffect_lib = starkeffect.so

#========rules=========================
compile: $(OBJS_WRAP) $(OBJS) $(INC)
	$(CXX) -fPIC $(SHARED_LIB) $(LIBS) $(LINKFLAGS) -o ../../lib/$(seffect_lib) $(OBJS)

./obj/wrap_%.o : wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/wrap_%.o : ./*/wrap_%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(WRAPPER_FLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

./obj/%.o : %.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;
	
./obj/%.o : ./*/%.cc $(INC)
	$(CXX) $(CXXFLAGS) $(INCLUDES_LOCAL) $(INCLUDES) -c $< -o $@;

clean:
	rm -rf ./obj/*.o
	rm -rf ./obj/*.os
	rm -rf ../../lib/$(seffect_lib)

