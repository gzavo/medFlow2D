cmake_minimum_required(VERSION 3.0)
project(medFlow2D)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/../bin/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

set(SOURCE_FILES
    inih/ini.c
    inih/ini.h
    lbm/boundaries.c
    lbm/boundaries.h
    lbm/lb.c
    lbm/lb.h
    platelet/cell.c
    platelet/cell.h
    utils/arraylist.c
    utils/arraylist.h
    utils/queue.c
    utils/queue.h
    utils/utils.c
    utils/utils.h
    utils/slog.c
    utils/slog.h
    utils/miniz.h
    utils/tinycthread.c
    utils/tinycthread.h
    common.h
    main.c)

#list( APPEND CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fopenmp")
#list( APPEND CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp")

list( APPEND CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
list( APPEND CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

add_executable(medFlow2D ${SOURCE_FILES})

if(APPLE)
    #set(CMAKE_C_COMPILER "/usr/local/bin/gcc")
    #set(CMAKE_CXX_COMPILER "/usr/local/bin/g++")
    #set(CMAKE_LINKER "/usr/local/bin/ld")
elseif(UNIX)
    target_link_libraries(medFlow2D m pthread)
endif()


#find_package(OpenMP)
#if (OPENMP_FOUND)
#    list (APPEND CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    list (APPEND CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    message("OpenMP found successfully!")
#endif()

