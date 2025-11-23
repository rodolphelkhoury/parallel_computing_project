#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

#include <mpi.h>

#include "Grid.hpp"
#include "SubGridInformation.hpp"
#include "SubGrid.hpp"

#include "Window.hpp"
#include "GridRenderer.hpp"
#include "ColorMap.hpp"

int main(int argc, char** argv) {
    int worldRank = -1;
    int worldSize = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    try {
        int totalX = 129;
        int totalY = 129;

        double lengthX = 20.0;
        double lengthY = 20.0;

        double dx = lengthX / (totalX - 1);
        double dy = lengthY / (totalY - 1);

        double alpha = 1.0;
        double dt = 0.2 * dx * dy / alpha;
        int iterations = static_cast<int>(200 / dt);

        double dirichletValue = 100.0;

        simulation::Grid grid(
            totalX, totalY,
            dx, dy,
            dt,
            alpha,
            iterations,
            simulation::EdgeType::DIRICHLET, // North
            simulation::EdgeType::DIRICHLET, // South
            simulation::EdgeType::NEUMANN,   // East
            simulation::EdgeType::NEUMANN,   // West
            dirichletValue
        );

        if (worldRank == 0) {
            std::cout << "Grid created successfully: "
                      << grid.getTotalGridCellsCountX() << " x "
                      << grid.getTotalGridCellsCountY() << "\n";
        }

        simulation::SubGridInformation subGridInfo(MPI_COMM_WORLD, totalX, totalY);
        simulation::SubGrid subGrid(grid, subGridInfo);

        visualization::Window* window = nullptr;
        visualization::GridRenderer* renderer = nullptr;
        if (worldRank == 0) {
            window = new visualization::Window(800, 800, "Heat Simulation");
            renderer = new visualization::GridRenderer(*window);
        }

        for (int iter = 0; iter < grid.getTotalIterations(); ++iter) {
            subGrid.exchangeGhostCells();
            subGrid.applyBoundaryConditions();
            subGrid.updateCellTemp();

            if (iter % 10 == 0) {
                auto localInterior = subGrid.getInteriorCells();
                int procCountX = subGridInfo.getNumberOfProcessesOnX();
                int procCountY = subGridInfo.getNumberOfProcessesOnY();
                int procCoordX = subGridInfo.getProcessCoordinatesX();
                int procCoordY = subGridInfo.getProcessCoordinatesY();

                int baseX = totalX / procCountX;
                int leftoverX = totalX % procCountX;
                int startX = procCoordX * baseX + std::min(procCoordX, leftoverX);

                int baseY = totalY / procCountY;
                int leftoverY = totalY % procCountY;
                int startY = procCoordY * baseY + std::min(procCoordY, leftoverY);

                if (worldRank == 0) {
                    std::vector<double> globalGrid(totalX * totalY, 0.0);

                    // Copy master local data
                    for (int i = 0; i < subGrid.getCellCountX(); ++i) {
                        for (int j = 0; j < subGrid.getCellCountY(); ++j) {
                            int localIdx = i * subGrid.getCellCountY() + j;
                            int globalIdx = (startY + j) * totalX + (startX + i);
                            globalGrid[globalIdx] = localInterior[localIdx];
                        }
                    }

                    // Receive from other processes
                    for (int rank = 1; rank < worldSize; ++rank) {
                        int remoteCellsX, remoteCellsY, remoteStartX, remoteStartY;
                        MPI_Recv(&remoteCellsX, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&remoteCellsY, 1, MPI_INT, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&remoteStartX, 1, MPI_INT, rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(&remoteStartY, 1, MPI_INT, rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        std::vector<double> remoteData(remoteCellsX * remoteCellsY);
                        MPI_Recv(remoteData.data(), remoteCellsX * remoteCellsY, MPI_DOUBLE, rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        for (int i = 0; i < remoteCellsX; ++i) {
                            for (int j = 0; j < remoteCellsY; ++j) {
                                int remoteIdx = i * remoteCellsY + j;
                                int globalIdx = (remoteStartY + j) * totalX + (remoteStartX + i);
                                globalGrid[globalIdx] = remoteData[remoteIdx];
                            }
                        }
                    }

                    double minTemp = 0.0;
                    double maxTemp = dirichletValue;
                    renderer->render(globalGrid, totalX, totalY, minTemp, maxTemp);

                    if (glfwWindowShouldClose(window->get())) {
                        break;
                    }

                } else {
                    // Send local data to master
                    int cellsX = subGrid.getCellCountX();
                    int cellsY = subGrid.getCellCountY();

                    MPI_Send(&cellsX, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                    MPI_Send(&cellsY, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                    MPI_Send(&startX, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
                    MPI_Send(&startY, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
                    MPI_Send(localInterior.data(), cellsX * cellsY, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
                }
            }
        }

        if (worldRank == 0) {
            delete renderer;
            delete window;
        }

    } catch (const std::exception& e) {
        std::cerr << "[Rank " << worldRank << "] Error: " << e.what() << "\n";
    }

    MPI_Finalize();
    return 0;
}
