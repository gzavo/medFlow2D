cmake_minimum_required(VERSION 3.2)
project(tracer2D)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../bin/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    utils/interpolation.h
    utils/miniz.h
    ini/ini.h
    ini/ini.c
    ini/INIReader.h
    ini/INIReader.cpp
#    utils/tinycthread.c
#    utils/tinycthread.h
    main.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


add_executable(tracer2D ${SOURCE_FILES})