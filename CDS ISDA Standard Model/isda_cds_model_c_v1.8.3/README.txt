/*
 * ISDA CDS Standard Model
 *
 * Copyright (C) 2009 International Swaps and Derivatives Association, Inc.
 * Developed and supported in collaboration with Markit
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the ISDA CDS Standard Model Public License.
 */


Package contents
----------------

The 'examples' directory contains an example Excel workbook demonstrating the library, and a simple program showing how to
call it from 'C'.

The 'excel' directory contains code for the Excel addin.

The 'lib' directory contains the actual analytics code.


Building
--------

To build the library go to the 'build' directory within the target's directory, then to the directory for the desired platform.
Typing make will invoke the build.

Note that for the library itself there is another level within the build directory depending on whether a static or dynamic
library is to be built.

Within the tests directory 'make retest' will cleanup after a previous test run, and 'make tests' will rerun all of the
regression tests.
