f90sources += boxlib_f.f90
f90sources += f2kcli$(f2kcli_suf).f90

f90sources += BLProfiler_f90.f90

f90sources += bl_constants.f90
f90sources += bl_IO.f90
ifdef PROF
  f90sources += bl_prof.f90
else
  f90sources += bl_prof_stubs.f90
endif
f90sources += bl_error.f90
f90sources += bl_mem_stat.f90
f90sources += bl_parmparse.f90
f90sources += bl_space.f90
f90sources += bl_stream.f90
f90sources += bl_string.f90
f90sources += bl_system.f90
f90sources += bl_timer.f90
f90sources += bl_timer_c.f90
f90sources += bl_types.f90

f90sources += box_f.f90
f90sources += knapsack.f90
f90sources += layout.f90
f90sources += boxarray_f.f90
f90sources += box_util.f90
f90sources += fab.f90
f90sources += multifab_f.f90
f90sources += fabio.f90
f90sources += bl_fabio_c.f90
f90sources += plotfile.f90
f90sources += filler.f90
f90sources += cluster_f.f90
f90sources += interp.f90
f90sources += fourth_order_interp_coeffs.f90
f90sources += bc.f90
f90sources += bc_functions.f90
f90sources += bndry_reg.f90

f90sources += ml_boxarray.f90
f90sources += ml_layout.f90
f90sources += ml_multifab.f90
f90sources += ml_restrict_fill.f90

f90sources +=    cc_restriction.f90
f90sources += nodal_restriction.f90
f90sources += ml_cc_restriction.f90
f90sources += ml_nd_restriction.f90
f90sources += nodal_neumann_bcs.f90
f90sources += nodal_stencil_bc.f90

f90sources += create_umac_grown.f90
f90sources += define_bc_tower.f90
f90sources += fillpatch.f90
f90sources += multifab_fill_ghost_cells.f90
f90sources += multifab_physbc_edgevel.f90
f90sources += multifab_physbc.f90

f90sources += list_box.f90
f90sources += sort_box.f90
f90sources += vector_i.f90
f90sources += sort_d.f90
f90sources += sort_i.f90

f90sources += ppm_util.f90

f90sources += cutcells.f90

f90sources += regrid.f90
f90sources += make_new_grids.f90
f90sources += tag_boxes.f90

ifdef PARTICLES
  f90sources += particles_f.f90
endif

ifndef MPI
  f90sources += parallel_stubs.f90
else
  f90sources += parallel.f90
  ifeq ($(findstring gfortran, $(FC)), gfortran)
     F90FLAGS += --std=legacy
  endif
endif

ifdef OMP
  f90sources += omp.f90
else
  f90sources += omp_stubs.f90
endif

csources += fabio_c.c
csources += timer_c.c
csources += ppm_util_c.c
csources += system_util_c.c

ifdef RANDOM
  f90sources += mt19937ar.f90
endif

f90sources += backtrace_f.f90
cxxsources += backtrace_c.cpp

f90sources += bl_random_f.f90
cxxsources += bl_random_c.cpp

f90sources += bl_fort_mod.f90

cxxsources += BL_MemPool.cpp
cxxsources += BL_CArena.cpp
cxxsources += BL_Arena.cpp

f90sources += bl_mempool_f.f90

VPATH_LOCATIONS += $(FBOXLIB_HOME)/Src/BaseLib
INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/Src/BaseLib

