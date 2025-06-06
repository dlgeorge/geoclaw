
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FFLAGS += -DNETCDF -lnetcdff -I/usr/include -O2 -fopenmp -g -fcheck=all -Wall#-L$(NETCDF4_DIR)#--g -fbacktrace -fbounds-check#-O0 -W -Wall -fbounds-check -fcheck=all -pedantic-errors -Wunderflow -fbacktrace -ffpe-trap=invalid,zero,overflow -g
LFLAGS += $(FFLAGS) -lnetcdff

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

NEWAMR = ../../new_fortran_amrclaw
NEWGEO = ../../new_fortran_geoclaw
# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

RIEMANN = $(CLAW)/riemann/src

MODULES = \

SOURCES = \
  $(RIEMANN)/rpn2_geoclaw.f \
  $(RIEMANN)/rpt2_geoclaw.f \
  $(RIEMANN)/geoclaw_riemann_utils.f \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make optionse
include $(CLAWMAKE)
# Construct the topography data
.PHONY: input all fgmax_plots

input:
	python make_input_files.py

fgmax_plots:
	python process_fgmax.py

all:
	$(MAKE) .exe
	$(MAKE) input
	$(MAKE) output
	$(MAKE) plots
	$(MAKE) fgmax_plots
### DO NOT remove this line - make depends on it ###
