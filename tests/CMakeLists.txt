find_package(GTest 1.17 QUIET)
if(NOT GTest_FOUND)
  message(STATUS "Google Test not found or version too old, fetching GTest")
  include(FetchContent)
  FetchContent_Declare(googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        v1.17.0
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()
if(NOT TARGET GTest::gtest_main)
    message(FATAL_ERROR "Google Test not found")
endif()
include(GoogleTest)

add_executable(mpz_tests
    test_mpz.cpp
    test_mpz1.cpp
)

target_compile_features(mpz_tests PRIVATE cxx_std_20)

if(Boost_FOUND)
  target_include_directories(mpz_tests PRIVATE ${Boost_INCLUDE_DIRS})
else()
  target_link_libraries(mpz_tests PRIVATE boost_multiprecision)
endif()

target_link_libraries(mpz_tests PRIVATE ${PROJECT_NAME} GTest::gtest_main)

gtest_discover_tests(mpz_tests)
