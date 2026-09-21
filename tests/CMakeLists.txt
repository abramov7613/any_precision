find_package(Boost 1.82 CONFIG REQUIRED)
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(mpz_tests
    test_mpz.cpp
)

target_include_directories(mpz_tests PRIVATE ${Boost_INCLUDE_DIRS})

target_compile_features(mpz_tests PRIVATE cxx_std_20)

target_link_libraries(mpz_tests PRIVATE
    z3mp
    GTest::gtest_main
)

gtest_discover_tests(mpz_tests)
