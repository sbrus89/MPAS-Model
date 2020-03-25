default: libs 

OS := $(shell uname -s)

CPPFLAGS += -DLMD94_OBL
DEBUG = -g
ifeq "$(Debug)" "true"
CPPFLAGS += -DDEBUG_OUTPUT
CXXFLAGS = -g -Wall -c $(CPPFLAGS)
else
CXXFLAGS = -O3 -Wall -c $(CPPFLAGS)
endif

F90=gfortran
F90FLAGS= -g -fcheck=all -fdefault-real-8 -Wno-character-truncation

LDFLAGS =
LINK = ${F90}

DEPFLAGS = -M

libs: libmixing.a

libmixing.a: cvmix_adc.o 
ifeq ($(OS),Linux)
	ar cr -o libmixing.a *.o 
endif
ifeq ($(OS),Darwin)
	libtool -static -o libmixing.a  *.o
endif
	ranlib libmixing.a

clean: 
	rm -f *.o *.mod libmixing.a

.cpp.o: 
	$(CXX) -c $(CXXFLAGS) $(EXTRA_INC) $<

%.o: %.f90
	$(F90) -c $(DEBUG) $(F90FLAGS) -I$(NETCDF)/include $<


