cmake_minimum_required(VERSION 2.8)
project(Tournament_Pivoting C)

find_package(PkgConfig REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
pkg_search_module(PLASMA plasma)

# Set the output folder where the program and libraries will be created
set(CMAKE_BINARY_DIR bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH lib)
SET(CMAKE_BUILD_TYPE DEBUG) 
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -lm -Wall")

# The following folder will be included
include_directories(include)

# Spral
set(SPRAL_DIR "/home/cseg/numanlys/scarf523/install/spral" CACHE PATH "Root directory containing Spral")

if( SPRAL_DIR )

  set(ENV{PKG_CONFIG_PATH} "${SPRAL_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

  include_directories(${SPRAL_DIR}/include)
  link_directories(${SPRAL_DIR} ${SPRAL_DIR}/lib/)
endif()


# Plasma
pkg_check_modules(Plasma plasma )

file(GLOB_RECURSE SOURCES control/*.c pivots/*.c smatrix/*.c )

add_library(Tournament_Pivoting STATIC ${SOURCES})
target_link_libraries(Tournament_Pivoting spral)
target_link_libraries(Tournament_Pivoting gfortran)
target_link_libraries(Tournament_Pivoting pthread)
target_link_libraries(Tournament_Pivoting )

# Verbosity
set(TP_VERBOSITY 0 CACHE STRING "Controls the verbosity. Default is 0.")

if (TP_VERBOSITY GREATER 0)
  add_definitions(-DTP_VERBOSITY=${TP_VERBOSITY})
endif()
 
message(STATUS " libs     ${PLASMA_LIBRARIES} ")
message(STATUS " includes ${PLASMA_INCLUDE_DIRS} ")

target_link_libraries(Tournament_Pivoting ${PLASMA_LIBRARIES})
target_include_directories(Tournament_Pivoting PUBLIC ${PLASMA_INCLUDE_DIRS})
target_compile_options(Tournament_Pivoting PUBLIC ${PLASMA_CFLAGS_OTHER})
link_directories(${PLASMA_LIBRARY_DIRS})

# Test driver: tp_test
add_executable(tp_test ${PROJECT_SOURCE_DIR}/test/test.c)

target_link_libraries(tp_test Tournament_Pivoting)
