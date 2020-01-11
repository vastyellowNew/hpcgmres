add_executable(
    d_nmpc_model_test
    d_nmpc_model_test.cpp
    ${NMPC_MODEL_DIR}/d_nmpc_model.c
)
target_include_directories(
    d_nmpc_model_test
    PUBLIC
    ${NMPC_MODEL_DIR}
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    d_nmpc_model_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    d_nmpc_model_test
    googletest
)
add_test(
    NAME d_nmpc_model_test
    COMMAND $<TARGET_FILE:d_nmpc_model_test>
)

add_executable(
    s_nmpc_model_test
    s_nmpc_model_test.cpp
    ${NMPC_MODEL_DIR}/s_nmpc_model.c
)
target_include_directories(
    s_nmpc_model_test
    PUBLIC
    ${NMPC_MODEL_DIR}
    ${GTEST_INCLUDE_PATH}
    ${GMOCK_INCLUDE_PATH}
    ${GTEST_LIBRARY_PATH}
)
target_link_libraries(
    s_nmpc_model_test
    PUBLIC
    GTest::GTest
    GTest::GMock
)
add_dependencies(
    s_nmpc_model_test
    googletest
)
add_test(
    NAME s_nmpc_model_test
    COMMAND $<TARGET_FILE:s_nmpc_model_test>
)