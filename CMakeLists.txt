cmake_minimum_required(VERSION 3.16)
project(moldyn LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(MOLDYN_BUILD_TESTS  "Build the test suite" ON)
option(MOLDYN_BUILD_PYTHON "Build the pybind11 bindings" OFF)

add_library(moldyn_warnings INTERFACE)
target_compile_options(moldyn_warnings INTERFACE
  -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion)

add_subdirectory(src/moldyn)

if(MOLDYN_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
