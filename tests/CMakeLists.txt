cmake_minimum_required(VERSION 3.1)
project(hpcgmres_unittest C CXX)
set(CMAKE_CXX_STANDARD 11)

enable_testing()
include(cmake/gtest.cmake)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/../src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../include)
set(NMPC_MODEL_DIR ${CMAKE_SOURCE_DIR}/../examples/cartpole/model)

add_library(coverage_config INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/nmpc_model)
add_subdirectory(${CMAKE_SOURCE_DIR}/common)
add_subdirectory(${CMAKE_SOURCE_DIR}/cgmres_initializer)
add_subdirectory(${CMAKE_SOURCE_DIR}/single_shooting_cgmres)