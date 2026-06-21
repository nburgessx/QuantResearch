###########################################################################
# Settings for Linux using gcc
###########################################################################
#
#  ISDA CDS Standard Model
#
#  Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
#  Developed and supported in collaboration with Markit
#
#  This program is free software: you can redistribute it and/or modify it
#  under the terms of the ISDA CDS Standard Model Public License.
#
###########################################################################

BUILD = linux

# Tools
CC         = /usr/bin/gcc
CXX        = /usr/bin/g++
DLL_LINK   = /usr/bin/g++
LIB_LINK   = /usr/bin/ar
DIFF       = cmp

# Standard system include directories and library path
SYS_INCLUDES = -I/usr/include
SYS_LIB_PATH = 

###########################################################################
# Well-known file extension macros and o/s specifc settings
###########################################################################
OBJ     = o
DLL_PFX = lib
DLL_EXT = so
LIB_PFX =
LIB_EXT = a
EXE_EXT = 
EXPORTS = exp
MAPFILE = mapfile
CLEAN   =

############################################################################
# Standard system libraries
############################################################################
SYS_LIBS = -lc

############################################################################
# Compiler flags
############################################################################
DEFINES      = -DUNIX -DLINUX -DVERSION="$(VERSION)"
C_CFLAGS     = -g -c $(DEFINES)
CXX_CFLAGS   = $(C_CFLAGS)

############################################################################
# Linker flags
############################################################################
LIB_LFLAGS   = -q $@ 
DLL_LFLAGS   = -shared -o $@
EXE_LFLAGS   = -lm -o $@

############################################################################
# System libraries
############################################################################
C_SYS_LIBS   = $(SYS_LIB_PATH) $(SYS_LIBS)
CXX_SYS_LIBS = $(SYS_LIB_PATH) $(SYS_LIBS)

