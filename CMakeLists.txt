cmake_minimum_required(VERSION 3.16)

project(CaPLa VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS  "-O3 -DNDEBUG -march=native")

add_subdirectory(libdivsufsort)
set(DIVSUFSORT_HEADER "${CMAKE_BINARY_DIR}/libdivsufsort/include/divsufsort.h")
file(READ "${DIVSUFSORT_HEADER}" file_content)
string(REPLACE "int32_t" "int64_t" file_content "${file_content}")
file(WRITE "${DIVSUFSORT_HEADER}" "${file_content}")
set_target_properties(mksary PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

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
