cmake_minimum_required(VERSION 2.8)
project(Tournament_Pivoting C)

find_package(PkgConfig REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
pkg_search_module(PLASMA QUIET plasma)

if ( NOT PLASMA_FOUND) 
  MESSAGE(FATAL_ERROR "PLASMA is not found. We use pkg to find PLAsMA. Please install PLASMA and include the path to the plasma.pc file int PKG_CONFIG_PATH environment variable")
endif()
# Set the output folder where the program and libraries will be created
set(CMAKE_BINARY_DIR bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH lib)
SET(CMAKE_BUILD_TYPE DEBUG) 
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall -O3 ")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall  -ggdb ")

# The following folder will be included
include_directories(include)

# Plasma
pkg_check_modules(Plasma REQUIRED plasma )

# Verbosity
set(TP_VERBOSITY 1 CACHE STRING "Controls the verbosity. Default is 1.")

if (TP_VERBOSITY GREATER 0)
  add_definitions(-DTP_VERBOSITY=${TP_VERBOSITY})
endif()

file(GLOB_RECURSE SOURCES control/*.c pivots/*.c smatrix/*.c )

add_library(Tournament_Pivoting STATIC ${SOURCES})
target_link_libraries(Tournament_Pivoting ${PLASMA_LDFLAGS})
target_link_libraries(Tournament_Pivoting gfortran)
target_link_libraries(Tournament_Pivoting pthread)

target_include_directories(Tournament_Pivoting PUBLIC ${PLASMA_INCLUDE_DIRS})
target_compile_options(Tournament_Pivoting PUBLIC ${PLASMA_CFLAGS_OTHER})

if( SPRAL_DIR )
  target_include_directories(Tournament_Pivoting PUBLIC ${SPRAL_DIR}/include)
  link_directories(${SPRAL_DIR} ${SPRAL_DIR}/lib/)
  target_link_libraries(Tournament_Pivoting spral)
  else ()
  message(FATAL_ERROR "SPRAL_DIR is not defined. Use -D SPRAL_DIR to provide the root directory of the SPRAL install")
endif()

# Test driver: tp_test
add_executable(tp_test ${PROJECT_SOURCE_DIR}/test/test.c)
add_executable(rewrite ${PROJECT_SOURCE_DIR}/test/rewrite.c)
add_executable(test_paje ${PROJECT_SOURCE_DIR}/test/test_paje.c)

target_link_libraries(tp_test Tournament_Pivoting)
target_link_libraries(rewrite Tournament_Pivoting)
target_link_libraries(test_paje Tournament_Pivoting)
