#include <iostream>
#include "simulation/include/Grid.hpp"

int main() {
    try {
        simulation::Grid grid(
            10, 10,            // totalX, totalY
            0.1, 0.1,          // dx, dy
            0.01,              // dt
            0.5,               // alpha
            100,               // total iterations
            simulation::EdgeType::DIRICHLET, // North
            simulation::EdgeType::NEUMANN,   // South
            simulation::EdgeType::DIRICHLET, // East
            simulation::EdgeType::NEUMANN,   // West
            100.0              // Dirichlet value
        );

        std::cout << "Grid created successfully!\n";
        std::cout << "Size: " << grid.getTotalGridX() << " x " << grid.getTotalGridY() << "\n";
        std::cout << "Cell size: " << grid.getCellSizeX() << " x " << grid.getCellSizeY() << "\n";
        std::cout << "Time step: " << grid.getTimeStep() << "\n";
        std::cout << "Thermal diffusivity: " << grid.getThermalDiffusivity() << "\n";
        std::cout << "Dirichlet value: " << grid.getDirichletValue() << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error creating grid: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
