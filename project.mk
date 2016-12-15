################################################################
# project name
################################################################

# This name is used for the distribution tar file.

project_name := vmc

################################################################
# modules -- list of directories in which to search
# for module.mk include files
################################################################

# Caution: The order in which modules are declares is important, since
# it is also used in linking.  The object/library file for the
# *caller* must precede the object/library file for the *callee* for
# successful linking.

################
# programs
################

modules += programs

################
# libraries
################

modules += libraries

################################################################
# extras -- list of extra files to be included
# in distribution tar file
################################################################

extras :=

################################################################
# additional project-specific make settings and rules
################################################################

# Intel MKL (Math Kernel Library)
CPPFLAGS += -mkl
LDLIBS += -mkl -liomp5 -lm -ldl

# Gnu Scientific Library
LDLIBS += -lgsl -lgslcblas
CPPFLAGS += -DHAVE_INLINE

# optimization mode
CPPFLAGS += -O3

# debugging mode
CXXFLAGS += -g -traceback
