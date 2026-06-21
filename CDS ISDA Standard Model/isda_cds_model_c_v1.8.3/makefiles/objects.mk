############################################################################
# Generic rules for building object files.
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

# Rule to compile C source files with '.c' extension in SRC_DIR.
%.$(OBJ): $(SRC_DIR)/%.c
	$(CC) -o $@ $(C_CFLAGS) $(INCLUDES) $(SYS_INCLUDES) "$<"

# Rule to compile C++ source files with '.cpp' extension in SRC_DIR.
%.$(OBJ): $(SRC_DIR)/%.cpp
	$(CXX) -o $@ $(CXX_CFLAGS) $(INCLUDES) $(SYS_INCLUDES) "$<"

# Rule to compile C source files with '.c' extension in current directory.
%.$(OBJ): %.c
	$(CC) -o $@ $(C_CFLAGS) $(INCLUDES) $(SYS_INCLUDES) "$<"

