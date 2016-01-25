cmake_minimum_required(VERSION 3.3)
project(LSModelCrypto)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Fundamentals.c
    Fundamentals.h
    LSbox.c
    LSbox.h
    Test.c)

add_executable(LSModelCrypto ${SOURCE_FILES})