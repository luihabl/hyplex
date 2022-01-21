cmake_minimum_required(VERSION 3.11)
project(hyplex)

# C++ version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native -Wall")

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

add_compile_definitions(GIT_VERSION="?.?")
link_directories(
    ${MPI_CXX_LINK_DIRS}
    ${CMAKE_SOURCE_DIR}/lib/yaml
    ${CMAKE_SOURCE_DIR}/lib/hypre/src/hypre/lib
)

add_executable(hyplex 
src/main.cpp
src/configuration.cpp
src/diagnostics.cpp
src/dsmc.cpp
src/fields.cpp
src/input-output.cpp
src/mcc.cpp
src/num-tools.cpp
src/particles-in-mesh.cpp
src/particles.cpp
src/random-numbers.cpp
src/rsolver.cpp
src/simulation.cpp
)

include_directories(
lib
lib/yaml
lib/hypre/src/hypre/include
)

target_link_libraries(hyplex
yaml-cpp
HYPRE
MPI::MPI_CXX
)