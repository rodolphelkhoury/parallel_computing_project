#include "SubGridInformation.hpp"

namespace simulation {

SubGridInformation::SubGridInformation(MPI_Comm communicator, int totalGridX, int totalGridY)
    : m_communicator(MPI_COMM_NULL),
      m_rank(-1),
      m_numberOfProcesses(0),
      m_northNeighborRank(MPI_PROC_NULL),
      m_southNeighborRank(MPI_PROC_NULL),
      m_eastNeighborRank(MPI_PROC_NULL),
      m_westNeighborRank(MPI_PROC_NULL)
{
    MPI_Comm_size(communicator, &m_numberOfProcesses);
    MPI_Comm_rank(communicator, &m_rank);

    m_processGridDimensions[0] = 0; // should be initialized to 0 before calling MPI_Dims_create, so MPI can understand that it needs to decide the dimensions
    m_processGridDimensions[1] = 0;
    int periods[2] = {0, 0};  // non-periodic in both dimensions, means sub-grids on the edges does not have neighbors on the opposite side

    MPI_Dims_create(m_numberOfProcesses, 2, m_processGridDimensions.data()); // let MPI decide dimensions

    // takes an existing communicator and builds a new communicator that is better structured for grid-based communication.
    MPI_Cart_create(
        communicator, // existing communicator
        2,
        m_processGridDimensions.data(),
        periods,
        1, // reorder allowed
        &m_communicator // new communicator
    );

    MPI_Comm_rank(m_communicator, &m_rank); // get rank in new communicator

    MPI_Cart_coords( // get coordinates of this process in the process grid
        m_communicator,
        m_rank,
        2,
        m_processCoordinates.data()
    );

    MPI_Cart_shift(m_communicator, 0, 1, &m_westNeighborRank, &m_eastNeighborRank);
    MPI_Cart_shift(m_communicator, 1, 1, &m_southNeighborRank, &m_northNeighborRank);
}

MPI_Comm SubGridInformation::getCommunicator() const {
    return m_communicator;
}

int SubGridInformation::getRank() const {
    return m_rank;
}

int SubGridInformation::getNumberOfProcesses() const {
    return m_numberOfProcesses;
}

std::array<int, 2> SubGridInformation::getProcessGridDimensions() const {
    return m_processGridDimensions;
}

std::array<int, 2> SubGridInformation::getProcessCoordinates() const {
    return m_processCoordinates;
}

int SubGridInformation::getNorthNeighborRank() const {
    return m_northNeighborRank;
}

int SubGridInformation::getSouthNeighborRank() const {
    return m_southNeighborRank;
}

int SubGridInformation::getEastNeighborRank() const {
    return m_eastNeighborRank;
}

int SubGridInformation::getWestNeighborRank() const {
    return m_westNeighborRank;
}

} // namespace simulation
