cmake_minimum_required(VERSION 2.8.9)

project( proj_cornflakes )

# The source files
file(GLOB SOURCES "*.c")
# Add in kernels to be included in the distribution
add_subdirectory(kernels)


# Make the static library
add_library( lib_cornflakes STATIC ${SOURCES} )
target_link_libraries( lib_cornflakes packaged_kernels m )


# Configure swig
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
SET(CMAKE_SWIG_FLAGS "")

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET_SOURCE_FILES_PROPERTIES(cornflakes_library_swig.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(cornflakes_library python cornflakes_library_swig.i ${SOURCES})
SWIG_LINK_LIBRARIES(cornflakes_library packaged_kernels m ${PYTHON_LIBRARIES})


# Put the swig files up a level into the python library
install(TARGETS ${SWIG_MODULE_cornflakes_library_REAL_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/..)
install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/cornflakes_library.py DESTINATION ${PROJECT_SOURCE_DIR}/..)