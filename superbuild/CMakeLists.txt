cmake_minimum_required(VERSION 2.8.9)
include(ExternalProject)

project(SUPERBUILD NONE)
enable_language(C)

# CITATIONS:
# Sarcasm, https://github.com/Sarcasm/cmake-superbuild/blob/master/cmake/SuperBuild.cmake#L16
# mfreiholz, https://github.com/mfreiholz/cmake-example-external-project/blob/master/CMakeLists-ExternalProjects.txt
# Nektar project, https://gitlab.nektar.info/nektar/nektar/commit/d27126171289eeba610eb7a1e97deede8ae53d9f

option(USE_LIS "Use LIS" OFF)
option(DOWNLOAD_LIS "Download LIS" ON)
option(USE_PETSC "Use PETSc" OFF)
option(DOWNLOAD_PETSC "Download PETSc" ON)

option(DOWNLOAD_CORNFLAKES "Pull cornflakes from github" OFF)
option(DOWNLOAD_POPCORN "Pull popcorn from github" ON)

set(PATH_TO_CORNFLAKES_SRC ${CMAKE_SOURCE_DIR}/../
  CACHE STRING "Location of cornflakes src")
set(PATH_TO_POPCORN_SRC ${CMAKE_SOURCE_DIR}/../../popcorn
  CACHE STRING "Location of popcorn src")

set (DEPENDENCIES)
set (EXTRA_CMAKE_ARGS)


if(CMAKE_TOOLCHAIN_FILE)
  set(TOOLARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DSWIG_EXECUTABLE=${SWIG_EXECUTABLE}
  -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
  -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
  )
endif(CMAKE_TOOLCHAIN_FILE)

message(STATUS "CC = $ENV{CC}")
message(STATUS "CMAKE = ${CMAKE_C_COMPILER}")

if(USE_LIS)
  ExternalProject_Add (lis
    #GIT_REPOSITORY "https://github.com/anishida/lis"
    URL http://www.ssisc.org/lis/dl/lis-1.7.13.tar.gz
    #SOURCE_DIR "${CMAKE_SOURCE_DIR}/3rdparty/lis"
    
    CONFIGURE_COMMAND ./configure
    CC=${CMAKE_C_COMPILER}
    --enable-fortran
    --enable-omp
    --enable-shared
    --prefix=${CMAKE_INSTALL_PREFIX}

    BUILD_COMMAND make
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 1
    EXCLUDE_FROM_ALL 1
    )
  list(APPEND DEPENDENCIES lis)
  list(APPEND EXTRA_CMAKE_ARGS -Dlis_path=${CMAKE_INSTALL_PREFIX})
endif(USE_LIS)
if(USE_PETSC)
  ExternalProject_Add(petsc
    URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.4.tar.gz
    BUILD_IN_SOURCE 1
    
    CONFIGURE_COMMAND ./configure
    CC=${CMAKE_C_COMPILER}
    --prefix=${CMAKE_INSTALL_PREFIX}
    --with-petsc-arch=c-opt
    --with-fc=gfortran
    --with-mpi=0
    --download-superlu
    
    #--with-cc=${CMAKE_C_COMPILER}
    #--with-cxx=${CMAKE_CXX_COMPILER}
    )
  list(APPEND DEPENDENCIES petsc)
  list(APPEND EXTRA_CMAKE_ARGS -DPETSC_DIR=${CMAKE_INSTALL_PREFIX} -DPETSC_ARCH=${c-opt})
endif(USE_PETSC)



ExternalProject_Add(cornflakes
  #GIT_REPOSITORY "https://github.com/afqueiruga/cornflakes"
  #GIT_TAG "develop"
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${PATH_TO_CORNFLAKES_SRC}

  DEPENDS ${DEPENDENCIES}
  
  CMAKE_ARGS ${TOOLARGS}
  -DEXTRA_INCLUDES=${EXTRA_INCLUDES}
  -DUSE_LIS=${USE_LIS} -DUSE_PETSC=${USE_PETSC}
  -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${EXTRA_CMAKE_ARGS}
  
  INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )
ExternalProject_Add(popcorn
  GIT_REPOSITORY "https://github.com/afqueiruga/popcorn"
  #SOURCE_DIR ${PATH_TO_POPCORN_SRC}
  #DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND sh lazyinstall.sh ${CMAKE_INSTALL_PREFIX}/pythonlib/
  BUILD_IN_SOURCE 1
  )
