cmake_minimum_required(VERSION 3.1)
project(hpcgmres_unittest C CXX)
set(CMAKE_CXX_STANDARD 11)
enable_testing()

include(cmake/gtest.cmake)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/../src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../include)
set(NMPC_MODEL_DIR ${CMAKE_SOURCE_DIR}/../examples/cartpole/model)
set(BLASFEO_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/../blasfeo/include)
set(BLASFEO_LIB_PATH ${CMAKE_SOURCE_DIR}/../blasfeo/lib)

link_directories(${BLASFEO_LIB_PATH})

add_subdirectory(${CMAKE_SOURCE_DIR}/common)
# add_subdirectory(${CMAKE_SOURCE_DIR}/cgmres_initializer)