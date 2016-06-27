
set(TESTS
  test_indexmap
  test_cfmat_bc
  test_index_set
  )

foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.c)
  if(USE_PETSC)
        target_compile_definitions(${TEST} PRIVATE -DUSE_PETSC)
  endif(USE_PETSC)
  target_link_libraries(${TEST} cornflakes ${CORNFLAKES_LIBRARIES})
  add_test(${TEST} ${TEST})
endforeach(TEST)
