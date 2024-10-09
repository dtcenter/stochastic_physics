.SUFFIXES:
.SUFFIXES: .F90 .o
.SUFFIXES: .F .o


FC = mpif90
FFLAGS_STOCH = -real-size 64 $(FFLAGS) 
FCINCLUDES_STOCH = $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90 
RM = rm
ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: stoch_physics_lib

dummy:
	echo "****** compiling stochastic physics ******"

OBJS = \
	kinddef.o \
	mpi_wrapper.o \
	stochy_internal_state_mod.o \
	stochy_namelist_def.o \
	spectral_transforms.o \
	mersenne_twister.o \
	stochy_patterngenerator.o \
	compns_stochy.o \
	stochy_data_mod.o \
	get_stochy_pattern.o \
        stochastic_physics.o \
        stochastic_physics_mpas.o	

stoch_physics_lib: $(OBJS)
	ar -ru libstochphys.a $(OBJS)

phys_interface: $(OBJS)

# DEPENDENCIES:
mpas_atmphys_camrad_init.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

clean:
	$(RM) *.o *.mod libstochphys.a

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS_STOCH) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES_STOCH) 

.F.o:
	echo $(FFLAGS_STOCH)
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS_STOCH) -c $*.F $(CPPINCLUDES) $(FCINCLUDES_STOCH) 
