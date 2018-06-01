FBOXLIB_HOME := ../..

NDEBUG := t
PROF   := 
MPI    := 
OMP    := t

COMP := gfortran

MKVERBOSE := t

include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
VPATH_LOCATIONS += .
FINCLUDE_LOCATIONS += .

# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := Src/BaseLib

#-----------------------------------------------------------------------------
# Fmpack is the list of all the GPackage.mak files that we need to
# include into the build system to define the list of source files.
#
# Fmlocs is the list of all the directories that we want to search
# for the source files in -- this is usually going to be the
# same as the list of directories containing GPackage.mak defined
# above.
#
# Fincs is the list of directories that have include files that
# we need to tell the compiler about.

# BoxLib
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(FBOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(FBOXLIB_HOME)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files

#  Note: GMakerules.mak will include '.' at the start of the
#  VPATH_LOCATIONS to first search in the problem directory
VPATH_LOCATIONS += $(Fmlocs)


# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#-----------------------------------------------------------------------------
# define the build instructions for the executable
main.$(suf).exe: $(objects)
ifdef MKVERBOSE
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS
else
	@echo "Linking $@ ..."
	@$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS
endif


#-----------------------------------------------------------------------------
# runtime parameter stuff (probin.f90)

# template used by write_probin.py to build probin.f90
PROBIN_TEMPLATE := ./probin.template

# list of the directories to search for _parameters files
PROBIN_PARAMETER_DIRS = ./

# list of all valid _parameters files for probin
PROBIN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
ifdef MKVERBOSE
	$(FBOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)"
else
	@$(FBOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)"
endif
	@echo " "


#-----------------------------------------------------------------------------
# include the BoxLib Fortran Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak

ifeq ($(findstring bint01.nersc.gov, $(HOSTNAMEF)), bint01.nersc.gov)
#  F90FLAGS += -no-prec-div -fno-alias -fimf-precision=low -fimf-domain-exclusion=15 -align array64byte -opt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0
  F90FLAGS += -fno-alias -align array64byte -opt-assume-safe-padding -opt-streaming-stores always -opt-streaming-cache-evict=0
endif

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 -- this is where the '::' in make comes
# in handy
clean::
	$(RM) probin.f90

