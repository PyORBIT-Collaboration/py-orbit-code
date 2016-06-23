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

#wrappers CC FLAGS
WRAPPER_FLAGS = -fno-strict-aliasing

#CXXFLAGS
CXXFLAGS += -fPIC

#shared library flags
#SHARED_LIB = -dynamiclib -undefined suppress -flat_namespace
SHARED_LIB = -shared

#tracker shared library
external_lib = intsystem.so

#LINKFLAGS += -L/System/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/config -lpython2.6

#========rules=========================
compile: $(OBJS_WRAP) $(OBJS) $(INC)
	$(CXX) -fPIC $(SHARED_LIB) $(LINKFLAGS) -o ../../lib/$(external_lib) $(OBJS)

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
	rm -rf ../../lib/$(external_lib)

