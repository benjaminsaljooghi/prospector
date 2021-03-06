project(prospector)

cmake_minimum_required (VERSION 3.19)

enable_language("CUDA")

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories(include)
include_directories(lib)

add_library(util "src/cas.cpp" "src/array_discovery.cpp" "src/cas_profiles.cpp" "src/crispr.cpp" "src/debug.cpp" "src/util.cpp")
set_property(TARGET util PROPERTY CXX_STANDARD 20)

add_library(prospector "src/prospector.cu")
set_property(TARGET prospector PROPERTY CUDA_ARCHITECTURES 75)
set_property(TARGET prospector PROPERTY CXX_STANDARD 17)

add_executable(exec "src/main.cpp")
set_property(TARGET exec PROPERTY CXX_STANDARD 20)
target_link_libraries(exec util prospector)
