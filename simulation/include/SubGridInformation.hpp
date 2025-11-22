#ifndef SIMULATION_INCLUDE_SUBGRIDINFORMATION_HPP
#define SIMULATION_INCLUDE_SUBGRIDINFORMATION_HPP

#include <mpi.h>
#include <array>

namespace simulation {

/**
 * @brief Represents the MPI sub-grid placement of a process within the
 * overall decomposition.
 *
 * Tracks the communicator, rank and total number of processes, the
 * 2D process-grid dimensions and the coordinates of this process within
 * that grid. Provides the ranks of the neighboring processes (north,
 * south, east, west) for halo exchanges.
 */
class SubGridInformation {
private:
    MPI_Comm m_communicator;
    int m_rank;
    int m_numberOfProcesses;
    std::array<int, 2> m_processGridDimensions; // number of processes on each dimension
    std::array<int, 2> m_processCoordinates; // coordinates of this process in the process grid
    int m_northNeighborRank;
    int m_southNeighborRank;
    int m_eastNeighborRank;
    int m_westNeighborRank;
public:
    explicit SubGridInformation(MPI_Comm communicator, int totalGridX, int totalGridY);

    [[nodiscard]] MPI_Comm getCommunicator() const;
    [[nodiscard]] int getRank() const;
    [[nodiscard]] int getNumberOfProcesses() const;
    [[nodiscard]] std::array<int, 2> getProcessGridDimensions() const;
    [[nodiscard]] int getNumberOfProcessesOnX() const;
    [[nodiscard]] int getNumberOfProcessesOnY() const;
    [[nodiscard]] std::array<int, 2> getProcessCoordinates() const;
    [[nodiscard]] int getProcessCoordinatesX() const;
    [[nodiscard]] int getProcessCoordinatesY() const;
    [[nodiscard]] int getNorthNeighborRank() const;
    [[nodiscard]] int getSouthNeighborRank() const;
    [[nodiscard]] int getEastNeighborRank() const;
    [[nodiscard]] int getWestNeighborRank() const;
    [[nodiscard]] bool hasNorthNeighbor() const;
    [[nodiscard]] bool hasSouthNeighbor() const;
    [[nodiscard]] bool hasEastNeighbor() const;
    [[nodiscard]] bool hasWestNeighbor() const;
};

} // namespace simulation

#endif // SIMULATION_INCLUDE_SUBGRIDINFORMATION_HPP
