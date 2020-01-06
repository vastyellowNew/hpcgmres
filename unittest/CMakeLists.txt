cmake_minimum_required(VERSION 3.1)
project(hpcgmres_unittest C CXX)
set(CMAKE_CXX_STANDARD 11)
enable_testing()

include(cmake/gtest.cmake)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/../src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../include)

add_subdirectory(
    ${CMAKE_SOURCE_DIR}/common
)