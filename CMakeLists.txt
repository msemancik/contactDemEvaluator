cmake_minimum_required(VERSION 3.31.6)
project(contactDem)

set(CMAKE_CXX_STANDARD 20)

add_executable(contactDem main.cpp
)
set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++") # For better distribution
