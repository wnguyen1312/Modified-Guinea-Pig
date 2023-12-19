cmake_minimum_required(VERSION 2.6)

project(guinea-pig)

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

include(guinea_definitions)

include(compilers)

# With this true we automatically include current src/binary directory
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

add_subdirectory(src)

# Some commands in testing require cmake 2.8
if(CMAKE_VERSION GREATER 2.7)
  enable_testing()
  include(CTest)
  add_subdirectory(testing)
endif()
