file(GLOB PACKAGED_KERNELS "*.c")
set_property(SOURCE ${PACKAGED_KERNELS} APPEND PROPERTY COMPILE_FLAGS -fPIC)
add_library(packaged_kernels STATIC ${PACKAGED_KERNELS})
target_link_libraries(packaged_kernels m)