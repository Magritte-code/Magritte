cmake_minimum_required(VERSION 2.4)
project(CheckFortran Fortran)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/result.cmake"
  "set(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
  "set(CMAKE_Fortran_COMPILER_ID \"${CMAKE_Fortran_COMPILER_ID}\")\n"
  "set(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
  "set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 \"${CMAKE_Fortran_COMPILER_SUPPORTS_F90}\")\n"
  )
