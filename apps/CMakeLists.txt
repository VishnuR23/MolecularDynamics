add_executable(moldyn_run moldyn_run.cpp)
target_link_libraries(moldyn_run PRIVATE moldyn moldyn_warnings)

# Record the git commit the binary was built from, so every CSV moldyn_run
# writes can be traced back to the exact code that produced it.
execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE MOLDYN_GIT_SHA
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT MOLDYN_GIT_SHA)
  set(MOLDYN_GIT_SHA "unknown")
endif()
target_compile_definitions(moldyn_run PRIVATE MOLDYN_GIT_SHA="${MOLDYN_GIT_SHA}")
