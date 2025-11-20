#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

#include <mpi.h>

#include "Grid.hpp"
#include "SubGridInformation.hpp"

#include "Window.hpp"
#include "GridRenderer.hpp"
#include "ColorMap.hpp"

int main(int argc, char** argv) {
    // ----------------------------------------------------
    // 1) Initialize MPI
    // ----------------------------------------------------
    MPI_Init(&argc, &argv);

    int worldRank = -1;
    int worldSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    if (worldRank == 0) {
        std::cout << "Running with " << worldSize << " MPI processes.\n";
    }

    // ----------------------------------------------------
    // 2) Test Grid
    // ----------------------------------------------------
    try {
        simulation::Grid grid(
            10, 10,            // totalX, totalY
            0.1, 0.1,          // dx, dy
            0.01,              // dt
            0.5,               // alpha
            100,               // iterations
            simulation::EdgeType::DIRICHLET, // North
            simulation::EdgeType::NEUMANN,   // South
            simulation::EdgeType::DIRICHLET, // East
            simulation::EdgeType::NEUMANN,   // West
            100.0              // Dirichlet
        );

        if (worldRank == 0) {
            std::cout << "Grid created successfully: "
                      << grid.getTotalGridX() << " x " << grid.getTotalGridY()
                      << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "[Rank " << worldRank << "] Grid creation error: "
                  << e.what() << "\n";
    }

    // ----------------------------------------------------
    // 3) Test SubGridInformation
    // ----------------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    simulation::SubGridInformation pos(MPI_COMM_WORLD, 10, 10);

    // Print info per rank
    std::cout << "-------------------------------------------\n";
    std::cout << "Rank " << pos.getRank() << " / " << pos.getNumberOfProcesses() << "\n";
    std::cout << "Process grid dims: ["
              << pos.getProcessGridDimensions()[0] << ", "
              << pos.getProcessGridDimensions()[1] << "]\n";
    std::cout << "My coords: ["
              << pos.getProcessCoordinates()[0] << ", "
              << pos.getProcessCoordinates()[1] << "]\n";
    std::cout << "Neighbors (N,S,E,W): "
              << pos.getNorthNeighborRank() << ", "
              << pos.getSouthNeighborRank() << ", "
              << pos.getEastNeighborRank()  << ", "
              << pos.getWestNeighborRank()  << "\n";

    MPI_Barrier(MPI_COMM_WORLD);

    // ----------------------------------------------------
    // 4) Visualization test (ONLY rank 0 to avoid many windows)
    // ----------------------------------------------------
    if (worldRank == 0) {
        std::cout << "Starting window on Rank 0...\n";

        visualization::Window* window = nullptr;
        try {
            window = new visualization::Window(640, 480, "Subgrid Visualization (Rank 0)");

            visualization::GridRenderer renderer(*window);

            const int nx = 16, ny = 12;
            std::vector<double> data(nx * ny);

            // Simple gradient
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    data[j * nx + i] = (double(i) / (nx - 1)) * 100.0;

            while (window->isValid() && !glfwWindowShouldClose(window->get())) {
                renderer.render(data, nx, ny, 0.0, 100.0);
                std::this_thread::sleep_for(std::chrono::milliseconds(16));
            }

        } catch (const std::exception& e) {
            std::cerr << "[Rank 0] Window/Renderer error: " << e.what() << "\n";
        }

        delete window;
    }

    MPI_Finalize();
    return 0;
}
