# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)
include_directories (../dependencies/eigen)
include_directories (../dependencies/pybind11)
include_directories (../dependencies/pybind11/include)
include_directories (../dependencies/voro++/src)

# Add Python package separately
find_package (PythonLibs REQUIRED)
# Add python dirs separately
include_directories (${PYTHON_INCLUDE_DIRS})

# Setup GoogleTest
enable_testing()
include(GoogleTest)
add_subdirectory("${PROJECT_SOURCE_DIR}/dependencies/googletest" "dependencies/googletest")

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


# set (SOURCE_FILES
#     test_raytracer.cpp
#     test_multigrid.cpp
#     test_levelpops.cpp
#     test_shortchar_order_0.cpp
#     test_feautrier_order_2.cpp
#     test_solver_lambda.cpp
#     test_imager.cpp
# )

set (SOURCE_FILES test_*.cpp)


if (OMP_PARALLEL)
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
        set (CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS}        ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif (OMP_PARALLEL)

if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)

target_link_libraries(Magritte voro++)

add_executable        (test_raytracer test_raytracer.cpp)
target_link_libraries (test_raytracer Magritte)

add_executable        (test_multigrid test_multigrid.cpp)
target_link_libraries (test_multigrid Magritte)

add_executable        (test_levelpops test_levelpops.cpp)
target_link_libraries (test_levelpops Magritte)

add_executable        (test_shortchar_order_0 test_shortchar_order_0.cpp)
target_link_libraries (test_shortchar_order_0 Magritte)

add_executable        (test_feautrier_order_2 test_feautrier_order_2.cpp)
target_link_libraries (test_feautrier_order_2 Magritte)

add_executable        (test_imager test_imager.cpp)
target_link_libraries (test_imager Magritte)

package_add_test      (test_solver_lambda test_solver_lambda.cpp)
target_link_libraries (test_solver_lambda Magritte)

if (OpenMP_CXX_FOUND)
    target_link_libraries (test_raytracer         OpenMP::OpenMP_CXX)
    target_link_libraries (test_multigrid         OpenMP::OpenMP_CXX)
    target_link_libraries (test_levelpops         OpenMP::OpenMP_CXX)
    target_link_libraries (test_shortchar_order_0 OpenMP::OpenMP_CXX)
    target_link_libraries (test_feautrier_order_2 OpenMP::OpenMP_CXX)
    target_link_libraries (test_solver_lambda     OpenMP::OpenMP_CXX)
    target_link_libraries (test_imager            OpenMP::OpenMP_CXX)
endif()

if (OMP_PARALLEL)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_link_libraries (test_raytracer         atomic)
        target_link_libraries (test_multigrid         atomic)
        target_link_libraries (test_levelpops         atomic)
        target_link_libraries (test_shortchar_order_0 atomic)
        target_link_libraries (test_feautrier_order_2 atomic)
        target_link_libraries (test_solver_lambda     atomic)
        target_link_libraries (test_imager            atomic)
    else ()
        target_link_libraries (test_raytracer         OpenMP::OpenMP_CXX)
        target_link_libraries (test_multigrid         OpenMP::OpenMP_CXX)
        target_link_libraries (test_levelpops         OpenMP::OpenMP_CXX)
        target_link_libraries (test_shortchar_order_0 OpenMP::OpenMP_CXX)
        target_link_libraries (test_feautrier_order_2 OpenMP::OpenMP_CXX)
        target_link_libraries (test_solver_lambda     OpenMP::OpenMP_CXX)
        target_link_libraries (test_imager            OpenMP::OpenMP_CXX)
    endif ()
endif ()
