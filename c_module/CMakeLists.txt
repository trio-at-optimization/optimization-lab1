cmake_minimum_required(VERSION 3.25)
project(untitled2)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

add_executable(untitled2 main.c)
