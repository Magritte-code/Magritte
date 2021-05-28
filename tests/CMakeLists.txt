# Set minimun required CMake version
cmake_minimum_required (VERSION 3.11)

# Set project title
project (Tests)

#Include pybind11
include_directories (../dependencies/pybind11)
include_directories (../dependencies/pybind11/include)

# Add Python package separately
find_package (PythonLibs REQUIRED)
# Add python dirs separately
include_directories (${PYTHON_INCLUDE_DIRS})

# Setup GoogleTest
enable_testing()
include(GoogleTest)
add_subdirectory("${PROJECT_SOURCE_DIR}/../dependencies/googletest" "dependencies/googletest")

# Marcro for creating tests
macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()


set (CMAKE_BUILD_TYPE Release)
# set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_CUDA_FLAGS " -g -G -x cu -O5") WARNING -G cfr issue:#112
#set (CMAKE_CUDA_FLAGS " -g")
#set (CMAKE_CUDA_FLAGS "-O5 -use_fast_math")
set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -w --extended-lambda --use_fast_math --expt-relaxed-constexpr")

# Setting to test for
set (PARACABS_USE_MESSAGE_PASSING false)
set (PARACABS_USE_MPI             false)
set (PARACABS_USE_MULTI_THREADING true)
set (PARACABS_USE_OPENMP          true)
set (PARACABS_USE_ACCELERATOR     true)
set (PARACABS_USE_CUDA            true)
set (PARACABS_USE_SYCL            false)
set (PARACABS_USE_SIMD            false)

if (PARACABS_USE_OPENMP)
    if (CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C                "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS          "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES      "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY   ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY  ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -latomic")
    else ()
        find_package (OpenMP REQUIRED)
        # Set compiler and linker flags
        set (CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX             "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS       "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES   "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY   ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY  ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -latomic")
    else ()
        find_package (OpenMP REQUIRED)
        # Set compiler and linker flags
        set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()


set (SOURCE_FILES
    test.cpp  
    test_datatypes.cpp
    test_multi_threading.cpp
)


if    (PARACABS_USE_CUDA)
    enable_language (CUDA)
    # Specify the accelerator architecture (currently we don't)
    set (CMAKE_CUDA_ARCHITECTURES OFF)
    # Set source file language to CUDA
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (PARACABS_USE_CUDA)


# Include Paracabs directory
include             (../CMakeLists.txt)
include_directories (../Paracabs      )

# Set directory to put executables in
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY tests)

package_add_test (test_multi_threading test_multi_threading.cpp)
# package_add_test (test_datatypes       test_datatypes.cpp      )

# Create executables
#      add_executable (test_message_passing test_message_passing.cpp)
#      add_executable (test_datatypes       test_datatypes.cpp      )
add_executable(test_accelerated_loops test.cpp)

#cuda_add_executable (test_accelerator test_accelerator.cu)
#cuda_add_executable (test_allocator   test_allocator.cu  )
#cuda_add_executable (test             test.cu            )

if    (PARACABS_USE_MPI)
#    target_link_libraries (test_message_passing ${MPI_C_LIBRARIES})
#    target_link_libraries (test_multi_threading ${MPI_C_LIBRARIES})
#    target_link_libraries (test_accelerator     ${MPI_C_LIBRARIES})
#    target_link_libraries (test_datatypes       ${MPI_C_LIBRARIES})
endif (PARACABS_USE_MPI)

if    (OPENMP_CXX_FOUND)
    target_link_libraries (test_accelerated_loops OpenMP::OpenMP_CXX)
endif (OPENMP_CXX_FOUND)


