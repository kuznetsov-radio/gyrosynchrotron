# -D LINUX in sources means "not Windows" (OpenMP vs MSVC PPL). It is set on
# Darwin and Linux only; Windows builds use gscodes.vcxproj, not this makefile.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
BINDIR = ../binaries

CXX ?= g++
COMMON_COMPILE = -c -std=c++11 -O3 -fPIC

DEPS = Arr_DF.h Coulomb.h DF.h ExtMath.h FF.h GS.h IDLinterface.h Messages.h MWmain.h Plasma.h Std_DF.h Zeta.h

ifeq ($(UNAME_S),Darwin)
LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /opt/homebrew/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
LIBOMP := /usr/local/opt/libomp
endif
ifeq ($(wildcard $(LIBOMP)/include/omp.h),)
$(error Homebrew libomp not found (need omp.h). Install with: brew install libomp)
endif
ifeq ($(UNAME_M),arm64)
ONAME = $(BINDIR)/MWTransferArr_arm64.so
else
ONAME = $(BINDIR)/MWTransferArr_x86_64.so
endif
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -Xpreprocessor -fopenmp -I$(LIBOMP)/include
LINK_FLAGS = -shared
LDLIBS = -L$(LIBOMP)/lib -lomp -Wl,-rpath,$(LIBOMP)/lib
else ifeq ($(UNAME_S),Linux)
ONAME = $(BINDIR)/MWTransferArr.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LINK_FLAGS = -shared
LDLIBS = -fopenmp
else ifeq ($(OS),Windows_NT)
$(error Windows builds use gscodes.vcxproj; this makefile supports Darwin and Linux only)
else
ONAME = $(BINDIR)/MWTransferArr.so
CXXFLAGS = $(COMMON_COMPILE) -D LINUX -fopenmp
LINK_FLAGS = -shared
LDLIBS = -fopenmp
endif

.DEFAULT_GOAL := $(ONAME)

OBJS = Arr_DF.o Coulomb.o ExtMath.o FF.o getparms.o GS.o IDLinterface.o \
       PythonInterface.o Messages.o Plasma.o MWmain.o Std_DF.o Zeta.o

$(ONAME): $(OBJS) | $(BINDIR)
	$(CXX) $(LINK_FLAGS) -o $(ONAME) $(OBJS) $(LDLIBS)

Arr_DF.o: Arr_DF.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) Arr_DF.cpp
Coulomb.o: Coulomb.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) Coulomb.cpp
ExtMath.o: ExtMath.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) ExtMath.cpp
FF.o: FF.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) FF.cpp
getparms.o: getparms.cpp
	$(CXX) $(CXXFLAGS) getparms.cpp
GS.o: GS.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) GS.cpp
IDLinterface.o: IDLinterface.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) IDLinterface.cpp
PythonInterface.o: PythonInterface.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) PythonInterface.cpp
Messages.o: Messages.cpp
	$(CXX) $(CXXFLAGS) Messages.cpp
MWmain.o: MWmain.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) MWmain.cpp
Plasma.o: Plasma.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) Plasma.cpp
Std_DF.o: Std_DF.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) Std_DF.cpp
Zeta.o: Zeta.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) Zeta.cpp

$(BINDIR):
	mkdir -p $@

.PHONY: install clean MWTransferArr

MWTransferArr: $(ONAME)

install: $(ONAME)

clean:
	rm -f *.o $(ONAME)
