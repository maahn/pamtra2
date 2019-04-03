# Only to create python signature file *pyf
# or to create fortran library (which is not required when using setup.py)

LIB=pyPamtraRadarSimulator

# Use these lines for Linux/g77
FC=g77
FC=gfortran
FFLAGS=-c -O2 -funroll-loops -fexpensive-optimizations -fPIC


all : 
	f2py radar_simulator.f90 radar_spectral_broadening.f90 random_module.f90 radar_spectrum.f90 rescale_spectra.f90 dia2vel.f90 report_module.f90 -m $(LIB)Lib -h $(LIB)Lib.pyf

OBJ=\
kinds.o\
report_module.o\
constants.o\
convolution.o\
dia2vel.o\
dsort.o\
rescale_spectra.o\
radar_simulator.o\
radar_spectral_broadening.o\
radar_spectrum.o\
random_module.o\
rho_air.o\
viscosity_air.o


%.o: %.f90
	$(FC) $(FFLAGS) $< -o $@


lib:	$(OBJ)
	ar -rcs lib$(LIB).a $(OBJ)

shared:$(OBJ)
	$(FC)  -shared -o lib$(LIB).so $(OBJ)

install:	lib$(LIB).a
	cp lib$(LIB).a /usr/local/lib
	#rm *.o

installshared:lib$(LIB).so
	cp lib$(LIB).so /usr/local/lib
	#rm *.o
	ldconfig 

clean: 
	 -rm -f -r *.o *.a *.so