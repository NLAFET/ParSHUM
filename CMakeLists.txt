cmake_minimum_required(VERSION 2.8)
project(ParSHUM C)

find_package(PkgConfig REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
pkg_search_module(PLASMA QUIET plasma)
find_package(MPI REQUIRED)

if ( NOT PLASMA_FOUND) 
  MESSAGE(FATAL_ERROR "PLASMA is not found. We use pkg to find PLAsMA. Please install PLASMA and include the path to the plasma.pc file int PKG_CONFIG_PATH environment variable")
endif()
# Set the output folder where the program and libraries will be created
set(CMAKE_BINARY_DIR bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH lib)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall -O3")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -lm -Wall  -ggdb3")

# The following folder will be included
include_directories(include)
# include_directories(/home/s/snakov/install/zoltan/include)

# Plasma
pkg_check_modules(Plasma REQUIRED plasma )

# Verbosity
set(ParSHUM_VERBOSITY 1 CACHE STRING "Controls the verbosity. Default is 1.")

if (ParSHUM_VERBOSITY GREATER 0)
  add_definitions(-DParSHUM_VERBOSITY=${ParSHUM_VERBOSITY})
endif()

file(GLOB_RECURSE SOURCES control/*.c pivots/*.c smatrix/*.c )

add_library(ParSHUM STATIC ${SOURCES})
target_link_libraries(ParSHUM ${PLASMA_LDFLAGS})
target_link_libraries(ParSHUM gfortran)
target_link_libraries(ParSHUM pthread)

target_include_directories(ParSHUM PUBLIC ${PLASMA_INCLUDE_DIRS})
target_compile_options(ParSHUM PUBLIC ${PLASMA_CFLAGS_OTHER})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

if( SPRAL_DIR )
  target_include_directories(ParSHUM PUBLIC ${SPRAL_DIR}/include)
  link_directories(${SPRAL_DIR} ${SPRAL_DIR}/lib/)
  target_link_libraries(ParSHUM spral)
  add_definitions(-DHAVE_SPRAL)
endif()

# if( ZOLTAN_DIR )
#   target_include_directories(ParSHUM PUBLIC ${ZOLTAN_DIR}/include)
#   link_directories(${ZOLTAN_DIR} ${ZOLTAN_DIR}/lib/)
#   target_link_libraries(ParSHUM zoltan)
#   else ()
#   message(FATAL_ERROR "ZOLTAN_DIR is not defined. Use -D ZOLTAN_DIR to provide the root directory of the Zoltan install")
# endif()

# if( PATOH_DIR )
#   target_include_directories(ParSHUM PUBLIC ${PATOH_DIR})
#   link_directories(${PATOH_DIR} ${PATOH_DIR})
#   target_link_libraries(ParSHUM patoh)
#   else ()
#   message(FATAL_ERROR "PATOH_DIR is not defined. Use -D PATOH_DIR to provide the root directory of the PaToH install")
# endif()

# Test driver: tp_test
add_executable(ParSHUM_simple ${PROJECT_SOURCE_DIR}/test/test.c)
add_executable(rewrite ${PROJECT_SOURCE_DIR}/test/rewrite.c)
#add_executable(test_paje ${PROJECT_SOURCE_DIR}/test/test_paje.c)
#add_executable(zoltan_test ${PROJECT_SOURCE_DIR}/test/zoltan.c)
#add_executable(ParSHUM_SBBD ${PROJECT_SOURCE_DIR}/test/SBBD.c)

target_link_libraries(ParSHUM_simple ParSHUM)
target_link_libraries(rewrite ParSHUM)
#target_link_libraries(test_paje ParSHUM)
#target_link_libraries(ParSHUM_SBBD ParSHUM)
#target_link_libraries(zoltan_test ParSHUM)
#target_link_libraries(zoltan_test ${MPI_C_LIBRARIES})
#target_link_libraries(ParSHUM_SBBD ${MPI_C_LIBRARIES})
#target_link_libraries(ParSHUM_SBBD gfortran)
