cmake_minimum_required(VERSION 3.20)
project(beam_profiler_sdk VERSION 1.0.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_STATIC_LIBS "Build static libraries" OFF)
option(BUILD_C_API "Build C API" ON)
option(BUILD_CPP_API "Build C++ API" ON)
option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS "Build unit tests" ON)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Windows-specific settings
if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Add subdirectories
add_subdirectory(src)
if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif() 