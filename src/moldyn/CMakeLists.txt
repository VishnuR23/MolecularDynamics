# Source files are appended here as each task lands.
set(MOLDYN_SOURCES
  core/placeholder.cpp
)

add_library(moldyn STATIC ${MOLDYN_SOURCES})
target_include_directories(moldyn PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(moldyn PRIVATE moldyn_warnings)
