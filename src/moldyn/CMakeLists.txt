# Source files are appended here as each task lands.
set(MOLDYN_SOURCES
  core/box.cpp
  core/system.cpp
  io/nist_config.cpp
  forcefield/lennard_jones.cpp
  neighbor/cell_list.cpp
  neighbor/verlet_list.cpp
  integrate/velocity_verlet.cpp
  integrate/langevin.cpp
  integrate/nose_hoover.cpp
  util/lattice.cpp
  analyze/rdf.cpp
  analyze/msd.cpp
  analyze/vacf.cpp
  analyze/block_average.cpp
)

add_library(moldyn STATIC ${MOLDYN_SOURCES})
target_include_directories(moldyn PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(moldyn PRIVATE moldyn_warnings)
