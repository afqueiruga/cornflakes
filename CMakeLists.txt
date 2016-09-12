cmake_minimum_required(VERSION 2.8.9)

project( proj_cornflakes )
enable_testing()

# The source files
set(SOURCES
  src/kernel.c
  src/dofmap.c
  src/indexmap.c
  src/indexset.c
  src/sparsity_pattern.c

  src/cfdata.c
  src/cfmat.c
  src/target.c

  src/toctic.c

  src/cfdata_bc.c
  src/cfmat_bc.c
  src/cfdata_default.c
  src/cfmat_default.c
  src/cfmat_csr.c
  src/cfdata_petsc.c
  src/cfmat_petsc.c
  src/cfdata_lis.c
  src/cfmat_lis.c
  src/cfmat_lis_2.c
  
  src/hypergraph.c
  src/assemble.c
  src/assemble_omp.c
  src/filter.c
  src/fill_sparsity.c
  src/spatialhash.c
  src/graphers.c
  src/util.c
  )
set( HEADERS
  src/cornflakes.h
  src/kernel.h
  src/dofmap.h
  src/indexmap.h
  src/indexset.h
  src/sparsity_pattern.h
  src/cfdata.h
  src/cfmat.h
  src/target.h
  
  #src/toctic.h
  
  src/cfdata_bc.h
  src/cfmat_bc.h
  src/cfdata_default.h
  src/cfmat_default.h
  src/cfmat_csr.h
  src/cfdata_petsc.h
  src/cfmat_petsc.h
  src/cfdata_lis.h
  src/cfmat_lis.h
  src/cfmat_lis_2.h

  src/hypergraph.h
  src/assemble.h
  src/assemble_omp.h
  src/filter.h
  src/fill_sparsity.h
  src/spatialhash.h
  src/graphers.h
  src/util.h
  )
set( PY_FILES
  cornflakes/__init__.py
  cornflakes/Assemble.py
  cornflakes/Dofmap.py
  cornflakes/GraphIO.py
  cornflakes/Hypergraph.py
  cornflakes/ParticlePlacers.py
  cornflakes/SpatialHash.py
  )
  
# Add in kernels to be included in the distribution
add_subdirectory(src/kernels)
include_directories(src/)



set(EXTRA_INCLUDES "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/" CACHE STRING "Extra includes; e.g. location to numpy/core/include/")
include_directories(${EXTRA_INCLUDES})
option(USE_PETSC "Enable PETSc backend for assembly" ON)
option(USE_LIS "Enable LIS backend for assembly" ON)

# Make the static library
add_library( cornflakes STATIC ${SOURCES} )
set(CORNFLAKES_LIBRARIES m)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(USE_PETSC)
  set(PETSC_DIR $ENV{PETSC_DIR})
  set(PETSC_ARCH ) #$ENV{PETSC_ARCH})

  #set(PETSC_DIR /home/afq/Downloads/petsc)
  #set(PETSC_ARCH arch-linux2-c-debug )
  target_compile_definitions(cornflakes PRIVATE -DUSE_PETSC)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
  find_package(PETSc REQUIRED)
  list(APPEND CORNFLAKES_INCLUDES ${PETSC_INCLUDES})
  list(APPEND CORNFLAKES_DEFINITIONS ${PETSC_DEFINITIONS} -DUSE_PETSC)
  list(APPEND CORNFLAKES_LIBRARIES ${PETSC_LIBRARIES})
endif(USE_PETSC)

if(USE_LIS)
  set(lis_path "/home/afq/Research/Distribute/Stonedist" CACHE STRING "Path to lis")
  set(lis_lib_location "${lis_path}/lib/liblis.a")
  add_library(lib_lis STATIC IMPORTED)
  set_property(TARGET lib_lis PROPERTY IMPORTED_LOCATION ${lis_lib_location})
  #target_compile_definitions(cornflakes PRIVATE -DUSE_LIS)
  list(APPEND CORNFLAKES_DEFINITIONS -DUSE_LIS)
  list(APPEND CORNFLAKES_INCLUDES ${lis_path}/include/)
  list(APPEND CORNFLAKES_LIBRARIES lib_lis)
endif(USE_LIS)


target_link_libraries(cornflakes ${CORNFLAKES_LIBRARIES} )
include_directories( ${CORNFLAKES_INCLUDES} )
target_compile_definitions(cornflakes PRIVATE ${CORNFLAKES_DEFINITIONS})

set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Ofast")

target_link_libraries( cornflakes packaged_kernels m )


# Configure swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "")

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(src/cornflakes_library_swig.i PROPERTIES SWIG_FLAGS "-includeall")
swig_add_module(cornflakes_library python src/cornflakes_library_swig.i ${SOURCES})

if(USE_PETSC)
  #SWIG_LINK_LIBRARIES(cornflakes_library ${PETSC_LIBRARIES})
endif(USE_PETSC)
swig_link_libraries(cornflakes_library packaged_kernels m ${PYTHON_LIBRARIES})
if(USE_LIS)
  #set_source_files_properties( ${swig_generated_file_fullname}
  #    PROPERTIES COMPILE_FLAGS "-DUSE_LIS")
endif(USE_LIS)

#
# The tests
#
add_subdirectory(units)

#
# Install stuffs
#
install(TARGETS cornflakes DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

set( PY_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/pythonlib/cornflakes)
install(TARGETS ${SWIG_MODULE_cornflakes_library_REAL_NAME} DESTINATION ${PY_INSTALL_DIR} OPTIONAL)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cornflakes_library.py DESTINATION ${PY_INSTALL_DIR} OPTIONAL)
install(FILES ${PY_FILES} DESTINATION ${PY_INSTALL_DIR})

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
# I don't like it, but I'll just do it
install(FILES src/numpy.i DESTINATION ${CMAKE_INSTALL_PREFIX}/include )

configure_file(cornflakes.cmake.in cornflakes.cmake @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cornflakes.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)
configure_file(env.conf.in env.conf @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/env.conf DESTINATION ${CMAKE_INSTALL_PREFIX})