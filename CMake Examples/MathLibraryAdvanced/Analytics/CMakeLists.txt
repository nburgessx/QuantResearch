add_executable(Analytics
    main.cpp
)

# Inherit all shared build settings via the interface target 'my_build_options'
target_link_libraries(Analytics
	PRIVATE
		Addition
		Subtraction
		my_build_options
)
