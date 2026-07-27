###
#   This is a simple makefile that builds the different
#   polyfuns libraries for the CHIDES library
###

FC  = gfortran
LD  = gfortran
FFLAGS = -c -O1 -fPIC
AR = ar

# gfortran -shared -fPIC Hermite.f90 -o test.so
#
# object files
#
OBJ = polymul.o polyfuns.o Hermite.o

# Name kept as .so (not .dylib) since Chides.jl loads "libchides.so"
LIB = libchides.so

ifeq ($(shell uname -s),Darwin)
    LDFLAGS += -L$(shell xcrun --show-sdk-path)/usr/lib
endif

.PHONY: chides_lib clean

chides_lib: $(OBJ)
	$(LD) -shared $(LDFLAGS) -o $(LIB) $(OBJ)
	mv $(LIB) ../
#
polymul.o: polymul.f90
	$(FC) $(FFLAGS) polymul.f90
#
Hermite.o: Hermite.f90
	$(FC) $(FFLAGS) Hermite.f90
#
polyfuns.o: polyfuns.f90
	$(FC) $(FFLAGS) polyfuns.f90
#
clean:
	/bin/rm -f *.o *.a *.mod ../$(LIB)
