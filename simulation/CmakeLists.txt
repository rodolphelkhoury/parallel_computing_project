add_library(simulation_module
    src/Grid.cpp
    src/SubGridInformation.cpp
)

target_include_directories(simulation_module
    PUBLIC include
)

target_link_libraries(simulation_module
    PRIVATE MPI::MPI_CXX
)
