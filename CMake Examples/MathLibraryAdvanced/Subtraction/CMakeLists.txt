# Create Subtraction Project as a Static Library
add_library(Subtraction
    subtract.h
	subtract.cpp
)

# Declare Header File Location
target_include_directories(Subtraction PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Inherit all shared build settings via the interface target 'my_build_options'
target_link_libraries(Subtraction PRIVATE my_build_options)