# Include directory containing source
include_directories (../src/)
include_directories (../dependencies)
include_directories (../dependencies/eigen)
include_directories (../dependencies/pybind11)
include_directories (../dependencies/pybind11/include)

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


set (SOURCE_FILES
    test_raytracer.cpp
    test_multigrid.cpp
    test_levelpops.cpp
    test_shortchar_order_0.cpp
    test_feautrier_order_2.cpp
    test_solver_lambda.cpp
    test_imager.cpp
)


if    (OMP_PARALLEL)
    find_package (OpenMP)
endif (OMP_PARALLEL)

if    (GPU_ACCELERATION)
    set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
endif (GPU_ACCELERATION)

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
