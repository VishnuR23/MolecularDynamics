# Generates moldyn_git_sha.hpp, the provenance stamp every CSV moldyn_run
# writes carries. Run as a script (cmake -P) by a custom target that is
# deliberately always out of date, so the revision is re-read on every
# build rather than frozen at configure time.
#
# Requires: MOLDYN_SOURCE_DIR (the repository root), MOLDYN_OUTPUT (the
# header to write).
#
# The stamp is the revision of the CODE, so a build from a modified tree
# is marked `<sha>-dirty` and cannot masquerade as the clean commit it came
# from. Dirtiness deliberately ignores results/: the experiments write into
# the repository, so regenerating one would otherwise mark every later run
# dirty for a reason unrelated to the code. This matches
# experiments/git_sha.sh and plot_common.py, so all three agree.

execute_process(
  COMMAND git describe --always
  WORKING_DIRECTORY ${MOLDYN_SOURCE_DIR}
  OUTPUT_VARIABLE MOLDYN_GIT_SHA
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE git_result)

execute_process(
  COMMAND git status --porcelain -- . ":(exclude)results"
  WORKING_DIRECTORY ${MOLDYN_SOURCE_DIR}
  OUTPUT_VARIABLE MOLDYN_GIT_DIRTY
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)

if(NOT git_result EQUAL 0 OR NOT MOLDYN_GIT_SHA)
  set(MOLDYN_GIT_SHA "unknown")
elseif(MOLDYN_GIT_DIRTY)
  set(MOLDYN_GIT_SHA "${MOLDYN_GIT_SHA}-dirty")
endif()

# configure_file only rewrites the output when its content actually
# changes, so an unchanged revision leaves the header's timestamp alone
# and moldyn_run.cpp is not recompiled on every build.
configure_file(${CMAKE_CURRENT_LIST_DIR}/moldyn_git_sha.hpp.in ${MOLDYN_OUTPUT} @ONLY)
