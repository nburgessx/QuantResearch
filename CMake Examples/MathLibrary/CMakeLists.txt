cmake_minimum_required(VERSION 3.20)

project(MathLibrary LANGUAGES CXX)

# ---- Language standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Targets ----
add_subdirectory(Addition)
add_subdirectory(Subtraction)
add_subdirectory(Analytics)

