# FBOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set FBOXLIB_HOME as an environment variable, this line will be ignored
FBOXLIB_HOME ?= ../../..

NDEBUG    := t
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
MKVERBOSE := t

include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
include ../../Source/GPackage.mak

VPATH_LOCATIONS += .	../../Source

include $(FBOXLIB_HOME)/Src/BaseLib/GPackage.mak

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak
