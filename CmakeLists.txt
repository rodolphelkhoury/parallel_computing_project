cmake_minimum_required(VERSION 3.16)

project(ParallelHeatSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

find_package(MPI REQUIRED)
find_package(OpenGL REQUIRED)

add_subdirectory(simulation)
add_subdirectory(visualization)

add_executable(heat_solver main.cpp)

target_link_libraries(heat_solver
    PRIVATE simulation_module
    PRIVATE visualization_module
    PRIVATE MPI::MPI_CXX     # <-- link MPI
)

# Include module headers
target_include_directories(heat_solver
    PRIVATE
        simulation/include
        visualization/include
)
