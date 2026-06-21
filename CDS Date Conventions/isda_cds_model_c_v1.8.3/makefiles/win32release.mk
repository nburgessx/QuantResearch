###########################################################################
# Generic settings for MSVC release build
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

include $(TOP)/makefiles/win32.mk
BUILD = win32release

############################################################################
# Compiler flags
############################################################################
# OPTIONS
# /MTd = Debug, multi-threaded, static, runtime library
# /W3  = Warning level
# /Od  = No optimisation
# /Z7  = No program database for debugging - stored in object file instead
# /GX  = Enable synchronous exception handling
# /G5  = Code optimisation for pentium
# /X   = Ignore standard include dirs
# /GR  = Enable C++ RTTI
# /c   = Produce object file only (do not link)
# /Zm1200 = Increase maximum compiler heap size by 1200%
#
# DEFINES
#	_MBCS = Microsoft multi-byte character set - See Using Generic-Text Mappings in MSDN.
#
DEFINES      = /DWIN32 /D_MBCS /DVERSION="$(VERSION)"
C_CFLAGS     = /c /MT /W3 /GX /Od /G5 /Z7 /X /Zm1200 /nologo $(DEFINES)
CXX_CFLAGS   = $(C_CFLAGS) /GR

############################################################################
# Linker flags
############################################################################
LIB_LFLAGS   = /out:"$@"
DLL_LFLAGS   = /out:"$@" /pdb:none /incremental:no /machine:i386 /nologo /subsystem:windows /dll 
EXE_LFLAGS   = /out:"$@" /pdb:none /incremental:no /machine:i386 /nologo /subsystem:console

############################################################################
# System libraries
############################################################################
C_SYS_LIBS   = $(SYS_LIB_PATH) libcmt.lib $(SYS_LIBS)
CXX_SYS_LIBS = $(SYS_LIB_PATH) libcpmt.lib libcmt.lib $(SYS_LIBS)
