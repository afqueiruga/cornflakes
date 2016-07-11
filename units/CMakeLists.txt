
set(TESTS
  test_indexmap
  test_cfmat_bc
  test_index_set
  test_sparsity_pattern
  test_fill_sparsity
  test_assemble_csr
  )

foreach(TEST ${TESTS})
  add_executable(${TEST} ${TEST}.c)
  if(USE_PETSC)
        target_compile_definitions(${TEST} PRIVATE -DUSE_PETSC)
  endif(USE_PETSC)
  target_link_libraries(${TEST} cornflakes ${CORNFLAKES_LIBRARIES})
  add_test(${TEST} ${TEST})
endforeach(TEST)
