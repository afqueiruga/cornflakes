cmake_minimum_required(VERSION 2.8.9)

project( proj_cornflakes )

# The source files
#file(GLOB SOURCES "*.c")
set(SOURCES
  kernel.c
  dofmap.c
  cfdata.c
  cfmat.c
  target.c
  
  cfdata_default.c
  cfmat_default.c
  cfdata_petsc.c
  #cfmat_petsc.c
  
  hypergraph.c
  assemble.c
  spatialhash.c
  graphers.c
  util.c
  )
# Add in kernels to be included in the distribution
add_subdirectory(kernels)

option(USE_PETSC "Enable PETSc backend for assembly" ON)
option(USE_LIS "Enable LIS backend for assembly" ON)

# Make the static library
add_library( lib_cornflakes STATIC ${SOURCES} )


if(USE_PETSC)
  set(PETSC_DIR /home/afq/Downloads/petsc)
  set(PETSC_ARCH arch-linux2-c-debug )
  target_compile_definitions(lib_cornflakes PRIVATE -DUSE_PETSC)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
  find_package(PETSc REQUIRED)
  include_directories(${PETSC_INCLUDES})
  add_definitions(${PETSC_DEFINITIONS})
  
  target_link_libraries(lib_cornflakes ${PPETSC_LIBRARIES})
endif(USE_PETSC)
if(USE_LIS)
  set(lis_path "${PROJECT_SOURCE_DIR}/../lis/" CACHE STRING "Path to lis")
  set(lis_lib_location "${lis_path}/src/.libs/liblis.a" CACHE STRING "Path to liblis.a for this architecture")
  add_library(lib_lis STATIC IMPORTED)
  set_property(TARGET lib_lis PROPERTY IMPORTED_LOCATION ${lis_lib_location})
  target_compile_definitions(lib_cornflakes PRIVATE -DUSE_LIS)
  include_directories(${lis_path}/include/)
endif(USE_LIS)


set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Ofast")

target_link_libraries( lib_cornflakes packaged_kernels m )


# Configure swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "")

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(cornflakes_library_swig.i PROPERTIES SWIG_FLAGS "-includeall")
swig_add_module(cornflakes_library python cornflakes_library_swig.i ${SOURCES})

if(USE_PETSC)
  #SWIG_LINK_LIBRARIES(cornflakes_library ${PETSC_LIBRARIES})
endif(USE_PETSC)
swig_link_libraries(cornflakes_library packaged_kernels m ${PYTHON_LIBRARIES})
if(USE_LIS)
  #set_source_files_properties( ${swig_generated_file_fullname}
  #    PROPERTIES COMPILE_FLAGS "-DUSE_LIS")
endif(USE_LIS)



# Put the swig files up a level into the python library
install(TARGETS ${SWIG_MODULE_cornflakes_library_REAL_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/..)
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/cornflakes_library.py DESTINATION ${PROJECT_SOURCE_DIR}/..)