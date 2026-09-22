# Generates moldyn_git_sha.hpp, the provenance stamp every CSV moldyn_run
# writes carries. Run as a script (cmake -P) by a custom target that is
# deliberately always out of date, so the revision is re-read on every
# build rather than frozen at configure time.
#
# Requires: MOLDYN_SOURCE_DIR (the repository root), MOLDYN_OUTPUT (the
# header to write).
#
# `git describe --always --dirty` -- not `rev-parse --short HEAD` -- so a
# build from a modified working tree is stamped `<sha>-dirty` and cannot
# masquerade as the clean commit it was derived from. A result stamped
# with a bare sha is a promise that the tree was exactly that commit.

execute_process(
  COMMAND git describe --always --dirty
  WORKING_DIRECTORY ${MOLDYN_SOURCE_DIR}
  OUTPUT_VARIABLE MOLDYN_GIT_SHA
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE git_result)

if(NOT git_result EQUAL 0 OR NOT MOLDYN_GIT_SHA)
  set(MOLDYN_GIT_SHA "unknown")
endif()

# configure_file only rewrites the output when its content actually
# changes, so an unchanged revision leaves the header's timestamp alone
# and moldyn_run.cpp is not recompiled on every build.
configure_file(${CMAKE_CURRENT_LIST_DIR}/moldyn_git_sha.hpp.in ${MOLDYN_OUTPUT} @ONLY)
