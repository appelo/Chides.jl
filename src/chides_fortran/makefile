###
#   This is a simple makefile that builds the different
#   polyfuns libraries for the CHIDES library
###

FC  = gfortran
LD  = gfortran
FFLAGS = -c -O1 -shared -fPIC
AR = ar 

# gfortran -shared -fPIC Hermite.f90 -o test.so
#
# object files
#
OBJ = polymul.o polyfuns.o Hermite.o

#
chides_lib: $(OBJ)
	gfortran -shared -o libchides.so $(OBJ) 
	mv libchides.so ../
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
	/bin/rm -f *.o *.a
