# Create Addition Project as a Static Library
add_library(Addition STATIC
     add.h
	 add.cpp
)

# Declare Header File Location
target_include_directories(Addition PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Inherit all shared build settings via the interface target 'my_build_options'
target_link_libraries(Addition PRIVATE my_build_options)
