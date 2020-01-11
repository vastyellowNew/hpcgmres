add_executable(
    d_time_varying_smooth_horizon_test 
    d_time_varying_smooth_horizon_test.cpp
    ${SRC_DIR}/common/d_time_varying_smooth_horizon.c
)
target_include_directories(
    d_time_varying_smooth_horizon_test 
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    d_time_varying_smooth_horizon_test 
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    d_time_varying_smooth_horizon_test 
    googletest
)
add_test(
    NAME d_time_varying_smooth_horizon_test 
    COMMAND $<TARGET_FILE:d_time_varying_smooth_horizon_test>
)

add_executable(
    s_time_varying_smooth_horizon_test 
    s_time_varying_smooth_horizon_test.cpp
    ${SRC_DIR}/common/s_time_varying_smooth_horizon.c
)
target_include_directories(
    s_time_varying_smooth_horizon_test 
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    s_time_varying_smooth_horizon_test 
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    s_time_varying_smooth_horizon_test 
    googletest
)
add_test(
    NAME s_time_varying_smooth_horizon_test 
    COMMAND $<TARGET_FILE:s_time_varying_smooth_horizon_test>
)


add_executable(
    d_givens_rotation_test
    d_givens_rotation_test.cpp
    ${SRC_DIR}/common/d_givens_rotation.c
)
target_include_directories(
    d_givens_rotation_test
    PUBLIC
    ${BLASFEO_INCLUDE_PATH}
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    d_givens_rotation_test
    PUBLIC
    blasfeo
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    d_givens_rotation_test
    googletest
)
add_test(
    NAME d_givens_rotation_test 
    COMMAND $<TARGET_FILE:d_givens_rotation_test>
)


add_executable(
    s_givens_rotation_test
    s_givens_rotation_test.cpp
    ${SRC_DIR}/common/s_givens_rotation.c
)
target_include_directories(
    s_givens_rotation_test
    PUBLIC
    ${BLASFEO_INCLUDE_PATH}
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    s_givens_rotation_test
    PUBLIC
    blasfeo
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    s_givens_rotation_test
    googletest
)
add_test(
    NAME s_givens_rotation_test 
    COMMAND $<TARGET_FILE:s_givens_rotation_test>
)


add_executable(
    d_memory_manager_test
    d_memory_manager_test.cpp
    ${SRC_DIR}/common/d_memory_manager.c
)
target_include_directories(
    d_memory_manager_test
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    d_memory_manager_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    d_memory_manager_test
    googletest
)
add_test(
    NAME d_memory_manager_test
    COMMAND $<TARGET_FILE:d_memory_manager_test>
)

add_executable(
    s_memory_manager_test
    s_memory_manager_test.cpp
    ${SRC_DIR}/common/s_memory_manager.c
)
target_include_directories(
    s_memory_manager_test
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    s_memory_manager_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    s_memory_manager_test
    googletest
)
add_test(
    NAME s_memory_manager_test
    COMMAND $<TARGET_FILE:s_memory_manager_test>
)

add_executable(
    d_linear_algebra_test
    d_linear_algebra_test.cpp
    ${SRC_DIR}/common/d_linear_algebra.c
    ${SRC_DIR}/common/d_memory_manager.c
)
target_include_directories(
    d_linear_algebra_test
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    d_linear_algebra_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    d_linear_algebra_test
    googletest
)
add_test(
    NAME d_linear_algebra_test
    COMMAND $<TARGET_FILE:d_linear_algebra_test>
)

add_executable(
    s_linear_algebra_test
    s_linear_algebra_test.cpp
    ${SRC_DIR}/common/s_linear_algebra.c
    ${SRC_DIR}/common/s_memory_manager.c
)
target_include_directories(
    s_linear_algebra_test
    PUBLIC
    ${INCLUDE_DIR}/common
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    s_linear_algebra_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    s_linear_algebra_test
    googletest
)
add_test(
    NAME s_linear_algebra_test
    COMMAND $<TARGET_FILE:s_linear_algebra_test>
)
