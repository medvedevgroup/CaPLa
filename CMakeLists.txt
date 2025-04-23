cmake_minimum_required(VERSION 3.16)

project(CaPLa VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS  "-O3 -mbmi2 -msse4.2 -DNDEBUG -march=native ")

# set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build 64-bit version of libdivsufsort")
add_subdirectory(libdivsufsort)

add_library(cmd STATIC 
  src/cmdline.cpp
)
target_link_libraries(cmd PUBLIC)



add_executable(count_segments 
    src/count_segments.cpp    
)
target_link_libraries(count_segments PUBLIC cmd)

add_executable(process_fasta 
    src/process_fasta.cpp    
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(count_segments PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(create_random_genome
    src/create_random_genome.cpp    
)

add_executable(convert_kmer_types
    src/convert_kmer_types.cpp    
)
target_link_libraries(convert_kmer_types PUBLIC cmd)