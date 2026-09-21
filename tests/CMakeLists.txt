# Test files are appended here as each task lands.
set(MOLDYN_TEST_SOURCES
  test_main.cpp
  test_vec3.cpp
  test_box.cpp
  test_system.cpp
  test_nist_config.cpp
  test_nist_lj.cpp
  test_forces.cpp
  test_neighbor.cpp
  test_integrator.cpp
)

add_executable(moldyn_tests ${MOLDYN_TEST_SOURCES})
target_include_directories(moldyn_tests PRIVATE ${CMAKE_SOURCE_DIR}/third_party)
target_link_libraries(moldyn_tests PRIVATE moldyn moldyn_warnings)
target_compile_definitions(moldyn_tests PRIVATE
  MOLDYN_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

add_test(NAME moldyn_tests COMMAND moldyn_tests)
