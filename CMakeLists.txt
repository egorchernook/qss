cmake_minimum_required(VERSION 3.16)
project(qss VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(src)
add_subdirectory(src)

set(WARNINGS "-Wall -Wextra -Wshadow -Wconversion -Wpedantic")
add_compile_options(${WARNINGS} "-O3")
