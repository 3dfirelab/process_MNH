
# The compiler
FC = gfortran 

HOST=$(shell hostname)
ifeq ($(HOST), ubu)
	FC = gfortran 
	FCFLAGS = -fPIC -fcray-pointer -fdefault-real-8 -fdefault-double-8  -O0 -finit-real=snan -g -fbounds-check -Wall -fbacktrace 
	INCFLAGS = -I./include/ -I/opt/netcdf-4.3.0/gnu/include # UBU
	LDFLAGS  = /opt/tec360/lib/libtecio.a -lstdc++ -L/opt/netcdf-4.3.0/gnu/lib -lnetcdf -lnetcdff    # UBU
endif
ifeq ($(HOST), kro)
	FC = gfortran 
	FCFLAGS = -fPIC -fcray-pointer -fdefault-real-8 -fdefault-double-8  -O0 -finit-real=snan -g -fbounds-check -Wall -fbacktrace 
	INCFLAGS = -I./include/ -I/opt/netcdf-4.3.0/gnu/include # UBU
	LDFLAGS  = /opt/tec360/lib/libtecio.a -lstdc++ -L/opt/netcdf-4.3.0/gnu/lib -lnetcdf -lnetcdff    # UBU
endif
ifeq ($(HOST),pyro)
	FC = gfortran 
	FCFLAGS = -fPIC -fcray-pointer -fdefault-real-8 -fdefault-double-8  -O0 -g -fbounds-check -Wall -fbacktrace 
	INCFLAGS = -I./include/ -I/opt/netcdf-4.3.0/gnu/include 
	LDFLAGS  = /opt/tec360/lib/libtecio.a -lstdc++ -L/opt/netcdf-4.3.0/gnu/lib -lnetcdf -lnetcdff 
endif
ifeq  ($(findstring kraken, $(HOST) ), kraken)
	FC = gfortran 
	FCFLAGS = -fPIC -fcray-pointer -fdefault-real-8 -fdefault-double-8  -O0 -g -fbounds-check -Wall -fbacktrace -fallow-invalid-boz
	INCFLAGS = -I./include/ -I/home/globc/paugam/anaconda3/envs/mypy3_fs/include
	LDFLAGS  = /home/globc/paugam/Src/Process_MNH/src/writePlt/lib/libtecio.a -lstdc++ -L/home/globc/paugam/anaconda3/envs/mypy3_fs/lib -lnetcdf
endif
ifeq  ($(HOST), moritz)
	FC = x86_64-conda-linux-gnu-gfortran
	FCFLAGS = -fPIC -fcray-pointer -fdefault-real-8 -fdefault-double-8  -O0 -g -fbounds-check -Wall -fbacktrace 
	INCFLAGS = -I./include/ -I/home/paugam/anaconda3/envs/mypy3/include/
	LDFLAGS  = /home/paugam/Src/Process_MNH/src/writePlt/lib/libtecio.a -lstdc++ -L/home/paugam/anaconda3/envs/mypy3/lib -lnetcdf
endif

SRC = ./src/

# List of executables to be built within the package
PROGRAMS = writePlt

LIBSRC_SOURCES = $(shell find $(SRC) -name "*.f90" '!' -name $(PROGRAMS).f90 -exec basename \{} .f90 \;)

LIBSRC_FILES_ = $(addsuffix .f90  , ${LIBSRC_SOURCES})
LIBSRC_FILES  = $(addprefix $(SRC), ${LIBSRC_FILES_} )
LIBSRC_OBJECTS = $(addsuffix .o,${LIBSRC_SOURCES})

$(PROGRAMS): save_tecplot.o cons.o lambda2.o tred1.o tqlrat.o r1mach.o pythag.o
save_tecplot.o : cons.o lambda2.o
lambda2.0: tred1.o tqlrat.o r1mach.o pythag.o

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) $(INCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: $(SRC)%.f90
	$(FC) $(FCFLAGS) $(INCFLAGS) -c $<

_$(PROGRAMS).so: ${LIBSRC_OBJECTS}	
	f90wrap -m $(PROGRAMS) ${LIBSRC_FILES} -k include/kind_map -v 
	f2py --build-dir . -c -m _$(PROGRAMS) ${LIBSRC_OBJECTS} f90wrap*.f90 --f90flags='$(FCFLAGS)' $(INCFLAGS) $(LDFLAGS)

python: $(PROGRAMS) _$(PROGRAMS).so ${LIBSRC_OBJECTS} 
	@if python ./test/test.py > test.out ; then echo "## test OK "; else echo "## test 1 failed"; exit 1; fi
	@echo "########"

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
	rm -rf src.linux-x*
