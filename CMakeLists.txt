project(prospector)

cmake_minimum_required (VERSION 3.19)

set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
SET(CMAKE_CUDA_ARCHITECTURES 52)

enable_language("CUDA")

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


add_library(util "src/cas.cpp" "src/array_discovery.cpp" "src/cas_profiles.cpp" "src/crispr.cpp" "src/debug.cpp" "src/util.cpp" src/config.cpp)
set_property(TARGET util PROPERTY CXX_STANDARD 20)

add_library(prospector "src/prospector.cu")
set_property(TARGET prospector PROPERTY CXX_STANDARD 17)

add_executable(exec "src/main.cpp" include/config.h include/system.h)
set_property(TARGET exec PROPERTY CXX_STANDARD 20)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    SET(OPENMP OpenMP::OpenMP_CXX)
endif()


target_link_libraries(exec util prospector ${Boost_LIBRARIES} ${OPENMP})
