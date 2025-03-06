cmake_minimum_required(VERSION 3.10)
project(heartrate)

SET(CMAKE_CXX_FLAGS_DEBUG "-O0")
SET(CMAKE_C_FLAGS_DEBUG "-O0")

add_subdirectory(algorithms/dummy)

add_executable(heartrate main.c)

target_link_libraries(heartrate dummy)
