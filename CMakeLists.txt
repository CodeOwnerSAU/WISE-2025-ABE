cmake_minimum_required(VERSION 3.22)
project(ABE)

set(CMAKE_CXX_STANDARD 17)

add_executable(ABE ABE_main.cpp tools.h tools.cpp graph.h
        heap.h Point.cpp Point.h graph.cpp IG-Tree.cpp MBR.cpp MBR.h H2H.cpp Grid.cpp SKORP.cpp )
target_link_libraries(ABE metis)
find_package(Threads REQUIRED)
target_link_libraries(ABE Threads::Threads)