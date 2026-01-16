add_executable(Analytics
    main.cpp
)

target_link_libraries(Analytics
	PRIVATE
		Addition
		Subtraction
)

