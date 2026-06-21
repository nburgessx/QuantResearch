############################################################################
# Platform and build-type independent rule for building a dynamic library.
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
TARGET = $(DLL_PFX)$(LIBRARY).$(DLL_EXT)


############################################################################
# Default target.
############################################################################
default: $(TARGET)

############################################################################
# Build target dynamic library.
# Note that this command creates $(LIBRARY).lib as a by-product
############################################################################

$(TARGET): $(OBJECTS)
	$(DLL_LINK) $(DLL_LFLAGS) $(OBJECTS) $(SYS_LIB_PATH) $(SYS_LIBS)

############################################################################
# Well-known 'clean' target - Remove all components built for target.
############################################################################
clean:
	rm -rf $(OBJECTS) $(TARGET) $(LIB_PRX)$(LIBRARY).$(LIB_EXT) $(LIB_PRX)$(LIBRARY).$(EXPORTS) $(CLEAN)
