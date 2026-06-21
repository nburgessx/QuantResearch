############################################################################
# Platform and build-type independent rule for building a static library.
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

include $(TOP)/makefiles/objects.mk

# Convenience macros.
TARGET = $(LIB_PFX)$(LIBRARY).$(LIB_EXT)

############################################################################
# Default target.
############################################################################
default: $(TARGET)

############################################################################
# Build target static library.
############################################################################
$(TARGET): $(OBJECTS)
	$(LIB_LINK) $(LIB_LFLAGS) $(OBJECTS)

############################################################################
# Well-known 'clean' target - Remove all components built for target.
############################################################################
clean:
	rm -rf $(TARGET) $(OBJECTS) $(CLEAN)
	