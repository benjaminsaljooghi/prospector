cmake_minimum_required (VERSION 3.19)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
SET(CMAKE_CUDA_ARCHITECTURES 52)
project(prospector LANGUAGES CUDA CXX)

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Find Boost at /usr/local/boost_1_76_0
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.76.0)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(include)
include_directories(lib)


add_library(util src/cas.cpp src/array_discovery.cpp src/cas_profiles.cpp src/crispr.cpp src/debug.cpp src/util.cpp src/config.cpp src/cas_helpers.cpp)
set_property(TARGET util PROPERTY CXX_STANDARD 20)

add_library(prospector src/cu/prospector.cu src/cu/cas.cu src/cu/cuda_helpers.cu src/cu/cas_helpers.cu)
set_property(TARGET prospector PROPERTY CXX_STANDARD 20)

add_executable(exec src/main.cpp include/config.h include/system.h)
set_property(TARGET exec PROPERTY CXX_STANDARD 20)
set_target_properties(exec PROPERTIES OUTPUT_NAME "prospector")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    SET(OPENMP OpenMP::OpenMP_CXX)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(exec util prospector ${Boost_LIBRARIES} ${OPENMP})
