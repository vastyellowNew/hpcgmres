cmake_minimum_required(VERSION 3.1)
project(cartpolesim C CXX)


set(SRC_DIR ${CMAKE_SOURCE_DIR}/../../../src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../../../include)
set(NMPC_MODEL_DIR ${CMAKE_SOURCE_DIR}/../model)


add_library(
    d_cgmres_initializer
    STATIC
    ${SRC_DIR}/cgmres_initializer/d_cgmres_initializer.c
    ${SRC_DIR}/cgmres_initializer/d_mfgmres_for_cgmres_initializer.c
    ${SRC_DIR}/cgmres_initializer/d_inexact_newton_for_zero_horizon_ocp.c
    ${SRC_DIR}/cgmres_initializer/d_zero_horizon_ocp.c
    ${NMPC_MODEL_DIR}/d_nmpc_model.c
    ${SRC_DIR}/common/d_memory_manager.c
    ${SRC_DIR}/common/d_linear_algebra.c
    ${SRC_DIR}/common/d_givens_rotation.c
)
target_include_directories(
    d_cgmres_initializer
    PRIVATE
    ${NMPC_MODEL_DIR}
    ${INCLUDE_DIR}
)

add_library(
    d_single_shooting_cgmres
    STATIC
    ${SRC_DIR}/single_shooting_cgmres/d_single_shooting_cgmres.c
    ${SRC_DIR}/single_shooting_cgmres/d_mfgmres_for_single_shooting_cgmres.c
    ${SRC_DIR}/single_shooting_cgmres/d_single_shooting_continuation.c
    ${SRC_DIR}/single_shooting_cgmres/d_single_shooting_ocp.c
    ${NMPC_MODEL_DIR}/d_nmpc_model.c
    ${SRC_DIR}/common/d_time_varying_smooth_horizon.c
    ${SRC_DIR}/common/d_memory_manager.c
    ${SRC_DIR}/common/d_linear_algebra.c
    ${SRC_DIR}/common/d_givens_rotation.c
)
target_include_directories(
    d_single_shooting_cgmres
    PRIVATE
    ${NMPC_MODEL_DIR}
    ${INCLUDE_DIR}
)
target_link_libraries(
    d_single_shooting_cgmres
    PRIVATE
    d_cgmres_initializer
)

add_executable(
    dtest
    main.cpp
    ${CMAKE_SOURCE_DIR}/../../simulator/d_simulator.cpp
    ${CMAKE_SOURCE_DIR}/../../simulator/numerical_integrator.cpp
    ${CMAKE_SOURCE_DIR}/../../simulator/save_simulation_data.cpp
)
target_link_libraries(
    dtest
    PRIVATE
    d_single_shooting_cgmres
)
target_include_directories(
    dtest
    PRIVATE
    ${INCLUDE_DIR}
    ${NMPC_MODEL_DIR}
    ${CMAKE_SOURCE_DIR}/../../simulator
)