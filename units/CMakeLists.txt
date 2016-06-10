

add_executable(test_indexmap test_indexmap.c)
target_link_libraries(test_indexmap cornflakes ${CORNFLAKES_LIBRARIES})
add_test(test_indexmap test_indexmap)