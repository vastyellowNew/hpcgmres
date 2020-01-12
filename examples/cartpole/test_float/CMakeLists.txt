cmake_minimum_required(VERSION 3.1)
project(cartpolesim C)


set(SRC_DIR ${CMAKE_SOURCE_DIR}/../../../src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../../../include)
set(NMPC_MODEL_DIR ${CMAKE_SOURCE_DIR}/../model)


add_library(
    s_cgmres_initializer
    STATIC
    ${SRC_DIR}/cgmres_initializer/s_cgmres_initializer.c
    ${SRC_DIR}/cgmres_initializer/s_mfgmres_for_cgmres_initializer.c
    ${SRC_DIR}/cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp.c
    ${SRC_DIR}/cgmres_initializer/s_zero_horizon_ocp.c
    ${NMPC_MODEL_DIR}/s_nmpc_model.c
    ${SRC_DIR}/common/s_memory_manager.c
    ${SRC_DIR}/common/s_linear_algebra.c
    ${SRC_DIR}/common/s_givens_rotation.c
)
target_include_directories(
    s_cgmres_initializer
    PRIVATE
    ${NMPC_MODEL_DIR}
    ${INCLUDE_DIR}
)

add_library(
    s_single_shooting_cgmres
    STATIC
    ${SRC_DIR}/single_shooting_cgmres/s_single_shooting_cgmres.c
    ${SRC_DIR}/single_shooting_cgmres/s_mfgmres_for_single_shooting_cgmres.c
    ${SRC_DIR}/single_shooting_cgmres/s_single_shooting_continuation.c
    ${SRC_DIR}/single_shooting_cgmres/s_single_shooting_ocp.c
    ${NMPC_MODEL_DIR}/s_nmpc_model.c
    ${SRC_DIR}/common/s_time_varying_smooth_horizon.c
    ${SRC_DIR}/common/s_memory_manager.c
    ${SRC_DIR}/common/s_linear_algebra.c
    ${SRC_DIR}/common/s_givens_rotation.c
)
target_include_directories(
    s_single_shooting_cgmres
    PRIVATE
    ${NMPC_MODEL_DIR}
    ${INCLUDE_DIR}
)
target_link_libraries(
    s_single_shooting_cgmres
    PRIVATE
    s_cgmres_initializer
)

add_executable(
    s_main
    s_main.c
    ${SRC_DIR}/common/s_memory_manager.c
    ${SRC_DIR}/common/s_linear_algebra.c
)
target_link_libraries(
    s_main
    PRIVATE
    s_single_shooting_cgmres
)
target_include_directories(
    s_main
    PRIVATE
    ${NMPC_MODEL_DIR}
    ${INCLUDE_DIR}
)