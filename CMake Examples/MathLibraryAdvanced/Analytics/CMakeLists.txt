cmake_minimum_required(VERSION 3.15)
project(Analytics LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===============================
# Include directories
# ===============================
include_directories(
    ${CMAKE_SOURCE_DIR}/external/xlw/include
)

# ===============================
# Prebuilt xlw libraries
# ===============================
if(MSVC)
    # Select correct xlw lib based on build type
    set(XLW_LIB_RELEASE "${CMAKE_SOURCE_DIR}/external/xlw/lib/x64/xlw-vc143-mt-0_0_0.lib")
    set(XLW_LIB_DEBUG   "${CMAKE_SOURCE_DIR}/external/xlw/lib/x64/xlw-vc143-mt-gd-0_0_0.lib")

    # Select correct xlcall32 lib based on build type
    set(XLCALL32_LIB_RELEASE "${CMAKE_SOURCE_DIR}/external/xlw/xlcall32stub/Release/xlcall32.lib")
    set(XLCALL32_LIB_DEBUG   "${CMAKE_SOURCE_DIR}/external/xlw/xlcall32stub/Debug/xlcall32.lib")
endif()

# ===============================
# Analytics Excel Add-in (.xll)
# ===============================
add_library(Analytics MODULE
    ExcelFunctions.cpp
)

target_compile_definitions(Analytics PRIVATE EXCEL12)

# ===============================
# Link libraries
# ===============================
target_link_libraries(Analytics PRIVATE
    Addition
    Subtraction
    user32
    optimized ${XLW_LIB_RELEASE}
    debug     ${XLW_LIB_DEBUG}
    optimized ${XLCALL32_LIB_RELEASE}
    debug     ${XLCALL32_LIB_DEBUG}
)

# ===============================
# Output as .xll
# ===============================
set_target_properties(Analytics PROPERTIES
    SUFFIX ".xll"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
