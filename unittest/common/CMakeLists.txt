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
