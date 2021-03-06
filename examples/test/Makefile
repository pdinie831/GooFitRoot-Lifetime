#------------------------------------------------------------------------------
CXX=nvcc
LD=g++  
OutPutOpt = -o

#CXXFLAGS     = -O3 -arch=sm_20 -g -std=c++11
CXXFLAGS     = -O3 
CUDALIBDIR=lib64

UNAME=$(shell uname)
ifeq ($(UNAME), Darwin)
CUDALIBDIR=lib
CXXFLAGS+=-m64
else
# Adds library symbols to debug info 
# - not strictly required, nice to have.
# Turn it off if it gives you trouble. 
LDFLAGS += -rdynamic 
endif

ifeq ($(CUDALOCATION), )
CUDALOCATION = /usr/local/cuda/
endif
CUDAHEADERS = $(CUDALOCATION)/include/
PWD = $(shell /bin/pwd)
GOOFITDIR = $(PWD)/../../
ROOTSYS=/usr/local/root
#INCLUDES += -I$(CUDAHEADERS) -I$(GOOFITDIR) -I$(GOOFITDIR)/rootstuff -I$(GOOFITDIR)/PDFs/ 
#LIBS += -L$(GOOFITDIR)/rootstuff -lRootUtils 
INCLUDES += -I$(CUDAHEADERS) -I$(GOOFITDIR)  -I$(GOOFITDIR)/PDFs/ 
LIBS += 
ifeq ($(TARGET_OMP),)
   LIBS += -L$(CUDALOCATION)/$(CUDALIBDIR) -lcudart
else
   LIBS += -fopenmp
endif

# These are for user-level programs that want access to the ROOT plotting stuff, 
# not just the fitting stuff included in the GooFit-local ripped library. 
ifneq ($(ROOTSYS), )
   CXXFLAGS += -DHAVE_ROOT=1
   ROOT_INCLUDES = -I$(shell root-config --incdir)
   ROOT_LIBS     = $(shell root-config --libs) -lMinuit
endif
WRKDIR = $(GOOFITDIR)/wrkdir/

ifneq ($(TARGET_OMP),)
CXXFLAGS    += -fopenmp -DTHRUST_SYSTEM_BACKEND=THRUST_DEVICE_BACKEND_OMP 
LDFLAGS     += -fopenmp -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP -lgomp
endif 


THRUSTO		= $(WRKDIR)/Variable.o $(WRKDIR)/FitManager.o $(WRKDIR)/GooPdfCUDA.o $(WRKDIR)/Faddeeva.o $(WRKDIR)/FitControl.o $(WRKDIR)/PdfBase.o $(WRKDIR)/DataSet.o $(WRKDIR)/BinnedDataSet.o $(WRKDIR)/UnbinnedDataSet.o $(WRKDIR)/FunctorWriter.o 
ROOTRIPDIR	= $(GOOFITDIR)/rootstuff
ROOTUTILLIB	= $(ROOTRIPDIR)/libRootUtils.so 
PROGRAMS	= testGoofit

.SUFFIXES: 

all:	$(PROGRAMS)

%.o:	%.cu
	$(CXX) $(INCLUDES) $(ROOT_INCLUDES) $(DEFINEFLAGS) $(CXXFLAGS) -c -o $@ $<

testGoofit:	testGoofit.o 
			$(LD) $(LDFLAGS) $^ $(THRUSTO) $(LIBS) $(ROOT_LIBS) $(OutPutOpt) $@
			@echo "$@ done"

clean:
		@rm -f *.o core $(PROGRAMS) 
