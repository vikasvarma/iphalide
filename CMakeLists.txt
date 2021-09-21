# CMake preferences:
cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Project definition:
project(
    IPHALIDE 
    VERSION 0.1.0
)

find_package(Halide REQUIRED)

# Find dependent libraries:
find_library(HALIDE_LIB Halide /usr/local/lib)
find_library(PNG_LIB png /usr/local/libpng/lib)
find_library(JPEG_LIB jpeg /usr/local/libjpeg/lib)

message(STATUS "Halide library: ${HALIDE_LIB}")
message(STATUS "PNG library: ${PNG_LIB}")
message(STATUS "JPEG library: ${JPEG_LIB}")

# Create conv generator executable:
add_executable(conv_halide src/conv.cpp)
target_link_libraries(conv_halide PRIVATE Halide::Generator)

# Use genreator to create a library:
add_halide_library(conv FROM conv_halide
                   GENERATOR conv_generator
                   PARAMS auto_schedule=true
                   AUTOSCHEDULER Halide::Mullapudi2016)

# Create ssim generator executable:
add_executable(ssim_halide src/ssim.cpp)
target_link_libraries(ssim_halide PRIVATE Halide::Generator)

# Use genreator to create a library:
add_halide_library(ssim FROM ssim_halide
                   GENERATOR ssim_generator
                   PARAMS auto_schedule=true
                   AUTOSCHEDULER Halide::Mullapudi2016)


# Testing setup ----------------------------------------------------------------
add_executable(tests
    test/main.cpp
)

target_link_libraries(tests PRIVATE
    conv
    ssim
    ${HALIDE_LIB}
    ${PNG_LIB}
    ${JPEG_LIB}
    Halide::Tools
)