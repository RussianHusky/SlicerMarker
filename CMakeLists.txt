cmake_minimum_required(VERSION 3.16)

project(SlicerMarker LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(hatch_generator main.cpp)

if (MINGW)
    set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
endif()

if (MSVC)
    target_compile_options(hatch_generator PRIVATE /W4 /permissive-)
else()
    target_compile_options(hatch_generator PRIVATE -Wall -Wextra -Wpedantic)
endif()
