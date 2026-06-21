###########################################################################
# Generic settings for Win32 using MSVC
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

###########################################################################
# MSVC_HOME is expected not to have spaces in the path, or the spaces must 
# be escaped with a backslash.
###########################################################################
MSVC_HOME  = C:/Program\ Files/Microsoft\ Visual\ Studio/VC98
MSVC_BIN   = $(MSVC_HOME)/bin

# Tools
CC         = $(MSVC_BIN)/cl
CXX        = $(MSVC_BIN)/cl
DLL_LINK   = $(MSVC_BIN)/link
LIB_LINK   = $(MSVC_BIN)/lib
DIFF       = diff -qw

# Standard system include directories and library path
SYS_INCLUDES = -I$(MSVC_HOME)/include -I$(MSVC_HOME)/PlatformSDK/include
SYS_LIB_PATH = /libpath:$(MSVC_HOME)/PlatformSDK/lib /libpath:$(MSVC_HOME)/lib

###########################################################################
# Well-known file extension macros and o/s specifc settings
###########################################################################
OBJ     = obj
DLL_PFX =
DLL_EXT = dll
LIB_PFX =
LIB_EXT = lib
EXE_EXT = .exe
EXPORTS = exp
CLEAN   =

############################################################################
# Standard system libraries
############################################################################
SYS_LIBS = kernel32.lib wsock32.lib netapi32.lib advapi32.lib oldnames.lib comdlg32.lib comctl32.lib user32.lib
