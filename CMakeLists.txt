cmake_minimum_required(VERSION 3.29)
project(lp)

set(CMAKE_CXX_STANDARD 20)

add_executable(lp
        src/include/matrix.hpp
        src/main.cpp
        src/include/simplex.hpp
        src/include/parsing.hpp)
target_include_directories(lp PUBLIC src/include/)
