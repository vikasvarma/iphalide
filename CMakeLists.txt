# cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr/local/Cellar/llvm/12.0.1/lib/cmake/llvm -install ..

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
find_library(HALIDE_LIB Halide)
find_library(PNG_LIB png)
find_library(JPEG_LIB jpeg)

message(STATUS "Halide library: ${HALIDE_LIB}")
message(STATUS "PNG library: ${PNG_LIB}")
message(STATUS "JPEG library: ${JPEG_LIB}")

# Create utils library:
add_library(extern_utils STATIC src/util.cpp)
set_target_properties(extern_utils PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(extern_utils 
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    /usr/local/include
)
target_link_libraries(extern_utils PRIVATE ${HALIDE_LIB})

find_library(UTILS_LIB extern_utils)

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

# Create 2D DCT executable and library:
add_executable(dct2d_halide src/dct2d.cpp)
target_link_libraries(dct2d_halide PRIVATE Halide::Generator)
add_halide_library(dct2d FROM dct2d_halide
                   GENERATOR dct2d_generator
                   PARAMS auto_schedule=true
                   AUTOSCHEDULER Halide::Mullapudi2016)

add_executable(idct2d_halide src/idct2d.cpp)
target_link_libraries(idct2d_halide PRIVATE Halide::Generator)
add_halide_library(idct2d FROM idct2d_halide
                   GENERATOR idct2d_generator
                   PARAMS auto_schedule=false
                   AUTOSCHEDULER Halide::Mullapudi2016)

#-------------------------------------------------------------------------------
# Create bm3d generator executable:
add_executable(bm3d_halide src/bm3d/bm3d_halide.cpp)
target_link_libraries(bm3d_halide PRIVATE
     Halide::Generator
     extern_utils
)

# Use genreator to create a library:
add_halide_library(
    blockdct FROM bm3d_halide
    GENERATOR blockdct_generator
    PARAMS auto_schedule=true
    AUTOSCHEDULER Halide::Mullapudi2016
)

add_halide_library(
    blockdist FROM bm3d_halide
    GENERATOR blockdist_generator
    PARAMS auto_schedule=true
    AUTOSCHEDULER Halide::Mullapudi2016
)

add_halide_library(
    argsort FROM bm3d_halide
    GENERATOR argsort_generator
    PARAMS auto_schedule=true
    AUTOSCHEDULER Halide::Mullapudi2016
)

# Testing setup ----------------------------------------------------------------
add_executable(tests
    test/main.cpp
)

target_link_libraries(tests PRIVATE
    conv
    ssim
    dct2d
    idct2d
    blockdct
    blockdist
    argsort
    extern_utils
    ${HALIDE_LIB}
    ${PNG_LIB}
    ${JPEG_LIB}
    Halide::Tools
)