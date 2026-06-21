############################################################################
# Platform and build-type independent rule for building an executable.
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
TARGET = $(LIBRARY)$(EXE_EXT)

############################################################################
# Default target.
############################################################################
default: $(TARGET)

############################################################################
# Build target executable.
############################################################################
$(TARGET): $(OBJECTS) $(LIBS)
	$(DLL_LINK) $(EXE_LFLAGS) $(OBJECTS) $(LIBS) $(C_SYS_LIBS)
	chmod u+x $(TARGET)

############################################################################
# Build purified executable.
############################################################################
purify: $(OBJECTS) $(LIBS)
	purify $(DLL_LINK) $(EXE_LFLAGS) $(OBJECTS) $(LIBS) $(C_SYS_LIBS)
	mv $@ $(TARGET)
	chmod u+x $(TARGET)

############################################################################
# Well-known 'clean' target - Remove all components built for target.
############################################################################
clean:
	rm -rf $(TARGET) $(OBJECTS) $(CLEAN) *_pure_*.o
	
############################################################################
# Clean exe + main library.
############################################################################
cleanall:
	cd $(TOP)/lib/build/lib/$(BUILD); make clean
	rm -rf $(OBJECTS) $(TARGET) $(LIB_PRX)$(LIBRARY).$(LIB_EXT) $(LIB_PRX)$(LIBRARY).$(EXPORTS) $(CLEAN)

###########################################################################
# Rule for building main library
###########################################################################
$(TOP)/lib/build/lib/$(BUILD)/$(LIB_PRX)$(LIBRARY).$(LIB_EXT):
	cd $(TOP)/lib/build/lib/$(BUILD); make

