# Makefile for the df2d program.
# Author: Shayan Hoshyari
# Date: June 2015
#
# Available commands:
# make df2d: Creates the df2d executable file.

# Compilation flags---------------------------------------------------

# CPP linker and preprocessor flags are included here.
# If you have the include files  for your armadillo library somewhere
# that g++ cannot find, add -I/path/to/armadillo/hpp/files to them.
# -DARMA_NO_DEBUG means that armadillo does not do bound checking for
# the A(i,j) operators.
CFLAGS= -DARMA_NO_DEBUG
CPPFLAGS= -pedantic -DARMA_NO_DEBUG

# This flag contains the address of libraries used.
# If your armadillo library is located somewhere that g++ can find,
# simply set it to -larmadillo. If not, set it to
# -L/path/to/your/armadillo/(.so or .a)/file
MY_LIB = -larmadillo

# This flag shows the address where you want to put the df2d executable.
BINDIR = ../bin/

# Petsc build system.---------------------------------------------------
# You should have petsc installed and the
# environment variables PETSC_ARCH and PETSC_DIR configured
# correctly in order to use this makefile.
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Rules ------------.---------------------------------------------------

df2d: region.o jkfunc.o df2d.o error.o node.o formula.o \
mdata.o  bvertex.o mesh.o  driver.o visit_writer.o asciifile.o \
geom.o chkopts
	${CLINKER} -o ${BINDIR}$@  region.o jkfunc.o df2d.o error.o \
node.o formula.o mdata.o bvertex.o mesh.o driver.o  visit_writer.o \
asciifile.o geom.o ${PETSC_LIB} ${MY_LIB}
	${RM}  $@.o

clear:
	rm -rf *o *~
