#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <string>

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

    // lightweight timestamped logging helpers

    auto current_time_str = []() -> std::string {
        using namespace std::chrono;
        auto now = system_clock::now();
        auto t = system_clock::to_time_t(now);
        std::tm tm;
        localtime_r(&t, &tm);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    };

    auto log_info = [&](int rank, const std::string& msg) {
        std::cout << "[" << current_time_str() << "] [Rank " << rank << "] INFO: " << msg << std::endl;
    };

    auto log_error = [&](int rank, const std::string& msg) {
        std::cerr << "[" << current_time_str() << "] [Rank " << rank << "] ERROR: " << msg << std::endl;
    };

    log_info(worldRank, std::string("MPI initialized; world size = ") + std::to_string(worldSize));

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

        simulation::SubGridInformation subGridInfo(MPI_COMM_WORLD, totalX, totalY);
        simulation::SubGrid subGrid(grid, subGridInfo);

        log_info(worldRank, std::string("Created grid ") + std::to_string(totalX) + "x" + std::to_string(totalY));
        log_info(worldRank, std::string("Subgrid cells: ") + std::to_string(subGrid.getCellCountX()) + "x" + std::to_string(subGrid.getCellCountY()));

        visualization::Window* window = nullptr;
        visualization::GridRenderer* renderer = nullptr;
        if (worldRank == 0) {
            window = new visualization::Window(800, 800, "Heat Simulation");
            renderer = new visualization::GridRenderer(*window);
            log_info(0, "Created window and renderer");
        }

        log_info(worldRank, std::string("Starting iterations: ") + std::to_string(grid.getTotalIterations()));

        for (int iter = 0; iter < grid.getTotalIterations(); ++iter) {
            if (iter % 50 == 0) {
                log_info(worldRank, std::string("Iteration ") + std::to_string(iter));
            }
            subGrid.exchangeGhostCells();
            subGrid.applyBoundaryConditions();
            subGrid.updateCellTemp();

            if (iter % 10 == 0) {
                log_info(worldRank, std::string("Iteration ") + std::to_string(iter) + ": preparing to exchange/collect data");
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
                    log_info(0, std::string("Master assembling global grid for iter ") + std::to_string(iter));
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
                    log_info(0, std::string("Rendering iteration ") + std::to_string(iter));
                    renderer->render(globalGrid, totalX, totalY, minTemp, maxTemp);

                    if (glfwWindowShouldClose(window->get())) {
                        break;
                    }

                } else {
                    log_info(worldRank, std::string("Sending local data to master at iter ") + std::to_string(iter));
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
            log_info(0, "Shutting down visualization");
            delete renderer;
            delete window;
        }

    } catch (const std::exception& e) {
        log_error(worldRank, std::string("Exception caught: ") + e.what());
    }

    log_info(worldRank, "Finalizing MPI");
    MPI_Finalize();
    return 0;
}
