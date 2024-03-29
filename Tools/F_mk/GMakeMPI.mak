# MPI reconfiguration...

ifndef MPI
$(error THIS FILE SHOULD ONLY BE LOADED WITH MPI defined)
endif

CCSE_MACHINES := angilas anguirus atragon baragon battra biollante dogora ebirah gabara
CCSE_MACHINES += gaira gamera garuda ghidorah gigan
CCSE_MACHINES += gimantis godzilla gojira hedorah kamacuras kiryu kong kraken kumonga manda
CCSE_MACHINES += megalon moguera mothra muto naphta orga rodan tindalos varan minilla winkle zilla
# CCSE's naphta orga are not included in CCSE_MACHINES

ifeq ($(HOST), $(findstring $(HOST), $(CCSE_MACHINES)))
    include $(FBOXLIB_HOME)/Tools/F_mk/Make.mpi.ccse
endif

DEFAULT_MACHINES := artoo naphta orga posse rob

ifeq ($(HOST), $(findstring $(HOST), $(DEFAULT_MACHINES)))
  MPIHOME=/usr/lib/mpich
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_libraries += -lmpich -lpthread   # -lmpichf90 might be needed
endif


# Architecture specific changes...
ifeq ($(ARCH),AIX)
  F90 = mpxlf95$(rsuf)
  FC  = mpxlf$(rsuf)
  CC  = mpcc$(rsuf)
  CXX = mpCC$(rsuf)
  mpi_libraries = -lmpi
endif

ifeq ($(ARCH),Darwin)

  # specifics for maudib
  ifeq ($(findstring rawk, $(UNAMEN)), rawk)
    FC = mpif90
    F90 = mpif90
    CC = mpicc
    CXX = mpicxx

    MPI_HOME=$(shell dirname `mpicc --showme:libdirs | cut -d" " -f2`)
    MPIHOME=$(MPI_HOME)
    LIBRARIES += lmpi_mpifh    

  # attempt at general support for Mac; only works if MPIHOME env var is set
  # and applies some general defaults
  else ifdef MPIHOME
    FC = $(MPIHOME)/bin/mpif90
    F90 = $(MPIHOME)/bin/mpif90
    CC = $(MPIHOME)/bin/mpicc
    CXX = $(MPIHOME)/bin/mpicxx

  # otherwise break
  else
    $(error SORRY, no MPI specification for your Darwin/Mac machine; check BoxLib/Tools/F_mk/GMakeMPI.mak and/or set MPIHOME environment variable)
  endif

  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = $(MPIHOME)/lib

endif

ifeq ($(ARCH),Linux)
  ifeq ($(COMP),g95)
    F90FLAGS += -fno-second-underscore
    FFLAGS   += -fno-second-underscore
    override F_C_LINK := UNDERSCORE
  endif
  ifeq ($(COMP),PathScale)
    FC = mpif90
    F90 = mpif90
    ifdef MPIHOME
      mpi_include_dir = $(MPIHOME)/include
    endif
  endif  # end of PathScale
endif
#
# Try and catch the Lawrencium cluster ...
#
ifeq ($(findstring .scs, $(HOSTNAMEF)), .scs)
    FC  = mpif90
    F90 = mpif90
endif
#
# Host changes ....
#
ifeq ($(findstring grace, $(HOST)), grace)
    #
    # grace.nersc.gov
    #
    ifdef MPI
        CXX := CC -target=linux
        CC  := cc -target=linux
        FC  := ftn -target=linux
        F90 := ftn -target=linux
    endif
endif
ifeq ($(findstring edison, $(HOST)), edison)
    #
    # edison.nersc.gov
    #
    ifdef MPI
        CXX := CC
        CC  := cc
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring cori, $(HOST)), cori)
    #
    # cori.nersc.gov
    #
    ifdef MPI
        CXX := CC
        CC  := cc
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring titan, $(HOST)), titan)
    #
    # titan (Oak Ridge, OLCF machine)
		#
		# Cray machines require you use their compiler wrappers
		# even if you aren't using Cray compiler
    #
    ifdef MPI
        CXX := CC 
        CC  := cc 
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring summitdev, $(HOST)), summitdev)
    #
    # summitdev (Oak Ridge, OLCF machine)
		#
		# Cray machines require you use their compiler wrappers
		# even if you aren't using Cray compiler
    #
    ifdef MPI
        CXX := mpicc
        CC  := mpicxx
        FC  := mpif77
        F90 := mpif90
    endif
endif
ifeq ($(findstring chester, $(HOST)), chester)
    #
    # titan (Oak Ridge, OLCF machine)
		#
		# Cray machines require you use their compiler wrappers
		# even if you aren't using Cray compiler
    #
    ifdef MPI
        CXX := CC 
        CC  := cc 
        FC  := ftn
        F90 := ftn
    endif
endif
ifeq ($(findstring h2o, $(UNAMEN)), h2o)
    #
    # Blue Waters
    #
    ifdef MPI
        CXX := CC
        FC  := ftn
        F90 := ftn
    endif

    ifeq ($(COMP),Cray)
      FFLAGS += -hpgas_runtime
      F90FLAGS += -hpgas_runtime
      CFLAGS += -hpgas_runtime
      CXXFLAGS += -hpgas_runtime
    endif
endif
ifeq ($(findstring mira, $(UNAMEN)), mira)
    #
    # The BlueGene/Q at ALCF
    #
    ifeq ($(COMP),IBM)
      ifdef MPI
        ifdef OMP
          CXX := mpixlcxx_r
          FC  := mpixlf95_r
          F90 := mpixlf95_r
        else
          CXX := mpixlcxx
          FC  := mpixlf95
          F90 := mpixlf95
        endif
      endif
    endif
endif


ifeq ($(HOST),cfe3)
  mpi_libraries += -lmpi
endif

ifeq ($(findstring cims.nyu.edu, $(HOSTNAME)), cims.nyu.edu)
   # Courant machines require "module load" statements for gcc and MPI
   # Check the latest list with "module avail" and then update here
   
   #ifeq ($(findstring ifort, $(FC)), ifort)
   #Compiler-specific MPI stuff
   #endif
   # Intel MPI libs:
   # Make sure to source /opt/pkg/intel/oneapi/setvars.sh
   # or run "module load intel-oneapi"
   # This exports variable $(I_MPI_ROOT)
   ifdef I_MPI_ROOT
      MPIHOME=$(I_MPI_ROOT)
      mpi_libraries += -lmpifort -lmpicxx
   else
      # OpenMPI v4.1, compiled with gcc 11.2
      MPIHOME=/usr/local/stow/openmpi-4.1
      mpi_libraries += -lmpi -lmpi_mpifh
   endif   

   # Generic stuff:
   mpi_include_dir = $(MPIHOME)/include
   mpi_lib_dir = $(MPIHOME)/lib       

else ifeq ($(findstring pop-os, $(HOSTNAME)), pop-os) 
   # A. Donev pop-os system76 laptop

   ifdef I_MPI_ROOT
      MPIHOME=$(I_MPI_ROOT)
      mpi_libraries += -lmpifort -lmpicxx
   else
      # OpenMPI on ubuntu:
      MPIHOME=/usr/lib/x86_64-linux-gnu/openmpi
      mpi_libraries += -lmpi -lmpi_mpifh
   endif   

   # Generic stuff:
   mpi_include_dir = $(MPIHOME)/include
   mpi_lib_dir = $(MPIHOME)/lib       

else 
   # OpenMPI on ubuntu
   MPIHOME=/usr/lib/x86_64-linux-gnu/openmpi
   mpi_libraries += -lmpi -lmpi_mpifh

   mpi_include_dir = $(MPIHOME)/include
   mpi_lib_dir = $(MPIHOME)/lib       
endif

ifeq ($(HOST),greenstreet)
  MPIHOME=/usr/local/anag/pkg/mpich-1.2.6-intel90
  mpi_lib_dir = $(MPIHOME)/lib
  mpi_include_dir = $(MPIHOME)/include
  mpi_libraries += -lmpich
endif

# generic linux install with MPICH wrappers -- set the 
# BOXLIB_USE_MPI_WRAPPERS environment variable for this
ifdef BOXLIB_USE_MPI_WRAPPERS
    F90 = mpif90
    CXX = mpicxx
endif

ifeq ($(HOST),lookfar)
  MPIHOME=/usr/local
  mpi_include_dir = $(MPIHOME)/include
  mpi_lib_dir = -L$(MPIHOME)/lib
  mpi_libraries += -lmpich
  ifeq ($(COMP),Intel)
    FFLAGS += -assume 2underscores
    F90FLAGS += -assume 2underscores
    CFLAGS += -DBL_FORT_USE_DBL_UNDERSCORE
    CFLAGS += -UBL_FORT_USE_UNDERSCORE
  endif
endif

ifeq ($(HOSTNAME),hyades.ucsc.edu)
  F90 := mpiifort
  FC := mpiifort
  fC := mpiifort
endif


ifeq ($(findstring lanl, $(UNAMEN)), lanl)
  F90 := mpif90
  FC := mpif90
  fc := mpif90
  CXX := mpic++
  CC := mpicc
endif
