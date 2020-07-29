# Set minimun required CMake version
cmake_minimum_required (VERSION 2.8.12)

# Set project title
project(Tests)

# Include Absparac directory
include_directories(..)

# Set directory to put executables in
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY tests)

# Create executables
add_executable(test_1 test_1.cpp)