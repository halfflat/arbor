.PHONY: clean all
.SECONDARY:

topdir:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
srcdir:=$(topdir)src

VPATH:=$(srcdir)

all:: mocksim

sources:=mocksim.cc busywait.cc main.cc
objects:=$(patsubst %.cc,%.o,$(sources))
depends:=$(patsubst %.cc,%.d,$(sources))

OPTFLAGS?=-O3 -march=native
CXXFLAGS+=$(OPTFLAGS) -MMD -MP -std=c++14 -g
CPPFLAGS+=-I $(srcdir)
LDFLAGS+=-lhpx -lhpx_init -lhpx_iostream -lboost_program_options -lboost_system

-include $(depends)

mocksim: $(objects)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(objects)

realclean: clean
	rm -f mocksim $(depends)
