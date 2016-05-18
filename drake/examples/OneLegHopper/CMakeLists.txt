if (LONG_RUNNING_TESTS)
  add_matlab_test(NAME examples/OneLegHopper/runPassive COMMAND runPassive PROPERTIES TIMEOUT 3000)
endif()
