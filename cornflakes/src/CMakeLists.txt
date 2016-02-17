cmake_minimum_required(VERSION 2.8.9)

project( cornflakes )


FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

file(GLOB SOURCES "*.c")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "")

# Add in kernels to be included in the distribution
add_subdirectory(kernels)


# Configure swig
SET_SOURCE_FILES_PROPERTIES(cornflakes_library_swig.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(cornflakes_library python cornflakes_library_swig.i ${SOURCES})
SWIG_LINK_LIBRARIES(cornflakes_library packaged_kernels m ${PYTHON_LIBRARIES})


# Put the swig files up a level into the python library
install(TARGETS ${SWIG_MODULE_cornflakes_library_REAL_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/..)
install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/cornflakes_library.py DESTINATION ${PROJECT_SOURCE_DIR}/..)