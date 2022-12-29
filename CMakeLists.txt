cmake_minimum_required(VERSION 3.11)
project(hyplex)

include(FetchContent)


# Fetching external libraries
FetchContent_Declare(yamlcpp                             # name of the content
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git  # the repository
  GIT_TAG        yaml-cpp-0.7.0                          # the tag
)
set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(yamlcpp)


FetchContent_Declare(hypre                                   # name of the content
  GIT_REPOSITORY https://github.com/hypre-space/hypre.git    # the repository
  GIT_TAG        v2.27.0                                     # the tag
  SOURCE_DIR src
)

FetchContent_GetProperties(hypre)
if(NOT hypre_POPULATED)
  FetchContent_Populate(hypre)
  add_subdirectory(${hypre_SOURCE_DIR}/src ${hypre_BINARY_DIR})
endif()

# C++ version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wall")

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

add_compile_definitions(GIT_VERSION="?.?")
link_directories(
    ${MPI_CXX_LINK_DIRS}
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
    ${CMAKE_SOURCE_DIR}/src/external
)

target_link_libraries(hyplex
yaml-cpp
HYPRE
MPI::MPI_CXX
)