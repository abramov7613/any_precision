cmake_minimum_required(VERSION 3.20)

project(z3mp LANGUAGES CXX)

file(GLOB SRCS LIST_DIRECTORIES FALSE "${CMAKE_CURRENT_SOURCE_DIR}/src/util/*.cpp")

add_library(z3mp STATIC ${SRCS})

target_compile_features(z3mp PUBLIC cxx_std_20)

target_compile_definitions(z3mp
    PRIVATE
        _MP_INTERNAL
)

target_include_directories(z3mp
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
