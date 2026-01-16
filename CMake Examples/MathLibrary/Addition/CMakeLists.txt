add_library(Addition STATIC
     add.h
	 add.cpp
)

target_include_directories(Addition
    PUBLIC 
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
