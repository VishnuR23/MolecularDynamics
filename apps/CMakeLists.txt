add_executable(moldyn_run moldyn_run.cpp)
target_link_libraries(moldyn_run PRIVATE moldyn moldyn_warnings)

# Record the git revision the binary was built from, so every CSV
# moldyn_run writes can be traced back to the exact code that produced it.
#
# This is deliberately NOT an execute_process at configure time: that runs
# once, when the build tree is created, and then never again, so a binary
# built weeks and many commits later still stamps the revision the build
# directory was configured at. (That is exactly what happened before this
# was fixed -- every committed result carried a stale sha.)
#
# Instead a custom target regenerates a header on every build. Custom
# targets are always considered out of date, so the git query really does
# re-run each time; GitSha.cmake writes through configure_file, which only
# touches the header when the revision actually changed, so this costs one
# `git describe` per build and no extra recompilation.
set(MOLDYN_GIT_SHA_HEADER ${CMAKE_CURRENT_BINARY_DIR}/generated/moldyn_git_sha.hpp)

add_custom_target(moldyn_git_sha
  BYPRODUCTS ${MOLDYN_GIT_SHA_HEADER}
  COMMAND ${CMAKE_COMMAND}
          -DMOLDYN_SOURCE_DIR=${CMAKE_SOURCE_DIR}
          -DMOLDYN_OUTPUT=${MOLDYN_GIT_SHA_HEADER}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/GitSha.cmake
  COMMENT "Recording the git revision for the provenance stamp"
  VERBATIM)

add_dependencies(moldyn_run moldyn_git_sha)
target_include_directories(moldyn_run PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)

# The water structure experiment. Same provenance wiring as moldyn_run.
add_executable(water_rdf water_rdf.cpp)
target_link_libraries(water_rdf PRIVATE moldyn moldyn_warnings)
add_dependencies(water_rdf moldyn_git_sha)
target_include_directories(water_rdf PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
