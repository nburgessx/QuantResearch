cmake_minimum_required(VERSION 3.20)

# Solution name as shown in Visual Studio
project(MathLibraryAdvanced LANGUAGES CXX)

# ------------------------------------------------------------
# Interface target holding shared compiler & linker behaviour.
# They only propagates flags on targets that specify the INTERFACE when linking
# ------------------------------------------------------------
add_library(my_build_options INTERFACE)

# C++20 → generates /std:c++20 on MSVC and equivalent elsewhere
target_compile_features(my_build_options INTERFACE cxx_std_20)

# OpenMP is a compiler feature; CMake probes the active toolchain
find_package(OpenMP REQUIRED)
target_link_libraries(my_build_options INTERFACE OpenMP::OpenMP_CXX)

# Windows Only - MSVC specific flags reflected directly in Project Properties
if (MSVC)

    # /MP  → parallel compilation
    # /EHa → C++ exceptions + SEH
    # /Qpar → auto-parallelised loops (requires optimisation)
    # /GL- → explicitly disable Whole Program Optimization
    target_compile_options(my_build_options INTERFACE
        /MP
        /EHa
        /Qpar
        /GL-
    )

    # /LTCG:OFF maps to "Enable Link Time Code Generation = No"
    target_link_options(my_build_options INTERFACE
        /LTCG:OFF
    )

endif()

# Add library targets (each becomes a .vcxproj)
add_subdirectory(Analytics)
add_subdirectory(Addition)
add_subdirectory(Subtraction)

# Set Root directory so linker paths can be relative
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})