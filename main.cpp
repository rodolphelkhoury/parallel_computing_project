#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

#include "Grid.hpp"
#include "Window.hpp"
#include "GridRenderer.hpp"
#include "ColorMap.hpp"

int main() {
    // 1) Create a Grid to verify simulation code compiles/links.
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
    } catch (const std::exception& e) {
        std::cerr << "Error creating grid: " << e.what() << "\n";
    }

    // 2) Create a Window and render a test image continuously
    visualization::Window* window = nullptr;
    try {
        window = new visualization::Window(640, 480, "Simulation Window");
        if (!window->isValid()) {
            std::cout << "Window created but not valid; exiting.\n";
            delete window;
            return 1;
        }

        visualization::GridRenderer renderer(*window);

        const int nx = 16, ny = 12;
        std::vector<double> data(nx * ny);

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                data[j * nx + i] = (double(i) / (nx - 1)) * 100.0;
            }
        }

        while (window->isValid() && !glfwWindowShouldClose(window->get())) {
            renderer.render(data, nx, ny, 0.0, 100.0);
            std::this_thread::sleep_for(std::chrono::milliseconds(16)); // ~60 FPS
        }

    } catch (const std::exception& e) {
        std::cerr << "Window/Renderer error: " << e.what() << "\n";
    }

    delete window;
    return 0;
}
