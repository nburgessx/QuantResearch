###########################################################################
# Contains list of objects for the C example
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

INCLUDES = -I$(TOP)/lib/include/isda
SRC_DIR  = $(TOP)/examples/c/src

LIBS     = $(TOP)/lib/build/lib/$(BUILD)/$(LIB_PFX)$(LIBRARY).$(LIB_EXT)

OBJECTS  = main.$(OBJ)
