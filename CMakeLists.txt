cmake_minimum_required(VERSION 3.22)
project(SY_SKORP)

set(CMAKE_CXX_STANDARD 11)

add_executable(SY_SKORP main.cpp tools.h tools.cpp graph.h
        heap.h Point.cpp Point.h graph.cpp IG-Tree.cpp MBR.cpp MBR.h H2H.cpp Grid.cpp SKORP.cpp Top_K_SKORP.cpp CBIE.cpp)
target_link_libraries(SY_SKORP metis)
find_package(Threads REQUIRED)
target_link_libraries(SY_SKORP Threads::Threads)