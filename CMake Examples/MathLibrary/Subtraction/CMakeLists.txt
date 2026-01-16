add_library(Subtraction STATIC
    subtract.h
	subtract.cpp
)

target_include_directories(Subtraction
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
