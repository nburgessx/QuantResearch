###########################################################################
# Makefile for the Excel add-in.
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
TARGET = $(LIBRARY).xll

############################################################################
# Default target.
############################################################################
default: $(TARGET)

###########################################################################
# Rule for linking xll
# Note that this command creates $(LIBRARY).lib as a by-product
###########################################################################
$(TARGET): $(OBJECTS) $(LIBS)
	$(DLL_LINK) $(DLL_LFLAGS) $(OBJECTS) $(LIBS) $(C_SYS_LIBS)

############################################################################
# Well-known 'clean' target - Remove all components built for target.
############################################################################
clean:
	rm -rf $(OBJECTS) $(TARGET) $(LIB_PRX)$(LIBRARY).$(LIB_EXT) $(LIB_PRX)$(LIBRARY).$(EXPORTS) $(CLEAN)

############################################################################
# Clean xll + main library.
############################################################################
cleanall:
	cd $(TOP)/lib/build/lib/$(BUILD); make clean
	rm -rf $(OBJECTS) $(TARGET) $(LIB_PRX)$(LIBRARY).$(LIB_EXT) $(LIB_PRX)$(LIBRARY).$(EXPORTS) $(CLEAN)

###########################################################################
# Rule for building main library
###########################################################################
$(TOP)/lib/build/lib/$(BUILD)/$(LIB_PRX)$(LIBRARY).$(LIB_EXT):
	cd $(TOP)/lib/build/lib/$(BUILD); make
	