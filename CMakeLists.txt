cmake_minimum_required(VERSION 2.8)
project(ParSHUM C)
FIND_PACKAGE(OpenMP REQUIRED)

set(CMAKE_BINARY_DIR bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH lib)
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall -O3 -march=native ${OpenMP_C_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall -ggdb3 ${OpenMP_C_FLAGS}")

# The following folder will be included
include_directories(include)

file(GLOB_RECURSE SOURCES control/*.c pivots/*.c smatrix/*.c)

add_library(ParSHUM STATIC ${SOURCES})

if( SPRAL_DIR )
  target_include_directories(ParSHUM PUBLIC ${SPRAL_DIR}/include)
  link_directories( ${SPRAL_DIR}/lib/)
  target_link_libraries(ParSHUM spral)
  target_link_libraries(ParSHUM gfortran)
  add_definitions(-DHAVE_SPRAL)
endif()


if(MKL_LDFLAGS  AND MKL_CFLAGS)
  target_compile_options(ParSHUM PUBLIC ${MKL_CFLAGS})
  target_link_libraries(ParSHUM ${MKL_LDFLAGS})
else ()
  message(FATAL_ERROR "MKL was not provided . Use -DMKL_LDFLAGS and -DMKL_CFLAGS  to provide these flags. (provided by the MKL LINK ADVISOR: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) ")
endif()

add_executable(ParSHUM_simple ${PROJECT_SOURCE_DIR}/test/test.c)

target_link_libraries(ParSHUM_simple ParSHUM)
