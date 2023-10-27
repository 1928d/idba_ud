
cmake_minimum_required(VERSION 3.20)
project(idba_ud)

set(MAXK 4)
set(MAXSS 768)
add_compile_options(-Wall -pedantic -O3 -march=native -DMAXSS=${MAXSS} -DMAXK=${MAXK})

find_package(ZLIB)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(idba_lib STATIC
    misc/log.cpp
    misc/utils.cpp
    misc/options_description.cpp
    misc/hash_aligner.cpp
    assembly/assembly_utility.cpp
    assembly/local_assembler.cpp
    graph/contig_info.cpp
    graph/hash_graph.cpp
    graph/hash_graph_branch_group.cpp
    graph/scaffold_graph.cpp
    graph/contig_graph.cpp
    graph/contig_graph_branch_group.cpp
    basic/math.cpp
    sequence/sequence_io.cpp
    sequence/sequence_reader.cpp
    sequence/sequence_writer.cpp
    sequence/sequence.cpp
    sequence/compact_sequence.cpp
    sequence/short_sequence.cpp
)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(idba_ud idba_ud.cpp)
target_compile_features(idba_ud PRIVATE cxx_std_11)
target_link_libraries(idba_ud idba_lib ZLIB::ZLIB)
