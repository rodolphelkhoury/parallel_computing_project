#include "HeatEquation.hpp"
#include "SubGrid.hpp"
#include "Grid.hpp"

#include <mpi.h>
#include <vector>
#include <algorithm>

namespace simulation {

// Forward declaration of computeHeatUpdate
double computeHeatUpdate(double center, double west, double east, double north, double south,
                         double alpha, double dt, double dx);

SubGrid::SubGrid(const Grid& parentGrid, const SubGridInformation& subGridInfo)
    : m_parentGrid(parentGrid)
    , m_subGridInfo(subGridInfo)
{
    m_hasGhostWest  = subGridInfo.hasWestNeighbor();
    m_hasGhostEast  = subGridInfo.hasEastNeighbor();
    m_hasGhostNorth = subGridInfo.hasNorthNeighbor();
    m_hasGhostSouth = subGridInfo.hasSouthNeighbor();

    int totalCellsX = parentGrid.getTotalGridCellsCountX();
    int totalCellsY = parentGrid.getTotalGridCellsCountY();

    int procCountX = subGridInfo.getNumberOfProcessesOnX();
    int procCountY = subGridInfo.getNumberOfProcessesOnY();

    int processCoordinateX = subGridInfo.getProcessCoordinatesX();
    int processCoordinateY = subGridInfo.getProcessCoordinatesY();

    int baseCellsX     = totalCellsX / procCountX;
    int leftoverCellsX = totalCellsX % procCountX;

    int baseCellsY     = totalCellsY / procCountY;
    int leftoverCellsY = totalCellsY % procCountY;

    // this process gets 1 extra cell in X/Y if it's in the leftover region
    m_cellCountX = baseCellsX + (processCoordinateX < leftoverCellsX ? 1 : 0);
    m_cellCountY = baseCellsY + (processCoordinateY < leftoverCellsY ? 1 : 0);

    m_totalCellCountX = m_cellCountX + static_cast<int>(m_hasGhostWest) + static_cast<int>(m_hasGhostEast);
    m_totalCellCountY = m_cellCountY + static_cast<int>(m_hasGhostSouth) + static_cast<int>(m_hasGhostNorth);

    m_currentTemperature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);
    m_nextTemperature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);

    double dirichletValue = parentGrid.getDirichletValue();

    // Set Dirichlet boundaries using parentGrid edges
    if (processCoordinateY == procCountY - 1 && parentGrid.getEdgeNorth() == EdgeType::DIRICHLET) {
        int j = m_cellCountY;
        for (int i = 1; i <= m_cellCountX; ++i) {
            m_currentTemperature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (processCoordinateY == 0 && parentGrid.getEdgeSouth() == EdgeType::DIRICHLET) {
        int j = 1;
        for (int i = 1; i <= m_cellCountX; ++i) {
            m_currentTemperature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (processCoordinateX == 0 && parentGrid.getEdgeWest() == EdgeType::DIRICHLET) {
        int i = 1;
        for (int j = 1; j <= m_cellCountY; ++j) {
            m_currentTemperature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (processCoordinateX == procCountX - 1 && parentGrid.getEdgeEast() == EdgeType::DIRICHLET) {
        int i = m_cellCountX;
        for (int j = 1; j <= m_cellCountY; ++j) {
            m_currentTemperature[cellIndex(i, j)] = dirichletValue;
        }
    }
}

void SubGrid::updateCellTemp() {
    int processCoordinateX = m_subGridInfo.getProcessCoordinatesX();
    int processCoordinateY = m_subGridInfo.getProcessCoordinatesY();

    int procCountX = m_subGridInfo.getNumberOfProcessesOnX();
    int procCountY = m_subGridInfo.getNumberOfProcessesOnY();

    for (int i = 1; i <= m_cellCountX; ++i) {
        for (int j = 1; j <= m_cellCountY; ++j) {

            bool isNorth = (processCoordinateY == procCountY - 1) && (j == m_cellCountY);
            bool isSouth = (processCoordinateY == 0) && (j == 1);
            bool isWest  = (processCoordinateX == 0) && (i == 1);
            bool isEast  = (processCoordinateX == procCountX - 1) && (i == m_cellCountX);

            if (isNorth) {
                if (m_parentGrid.getEdgeNorth() == EdgeType::DIRICHLET) {
                    m_nextTemperature[cellIndex(i,j)] = m_parentGrid.getDirichletValue();
                } else {
                    int ghostY = m_totalCellCountY - 1;
                    m_currentTemperature[cellIndex(i, ghostY)] = m_currentTemperature[cellIndex(i, j)];
                    m_nextTemperature[cellIndex(i,j)] = m_currentTemperature[cellIndex(i,j)];
                }
                continue;
            }

            if (isSouth) {
                if (m_parentGrid.getEdgeSouth() == EdgeType::DIRICHLET) {
                    m_nextTemperature[cellIndex(i,j)] = m_parentGrid.getDirichletValue();
                } else {
                    int ghostY = 0;
                    m_currentTemperature[cellIndex(i, ghostY)] = m_currentTemperature[cellIndex(i,j)];
                    m_nextTemperature[cellIndex(i,j)] = m_currentTemperature[cellIndex(i,j)];
                }
                continue;
            }

            if (isWest) {
                if (m_parentGrid.getEdgeWest() == EdgeType::DIRICHLET) {
                    m_nextTemperature[cellIndex(i,j)] = m_parentGrid.getDirichletValue();
                } else {
                    int ghostX = 0;
                    m_currentTemperature[cellIndex(ghostX, j)] = m_currentTemperature[cellIndex(i,j)];
                    m_nextTemperature[cellIndex(i,j)] = m_currentTemperature[cellIndex(i,j)];
                }
                continue;
            }

            if (isEast) {
                if (m_parentGrid.getEdgeEast() == EdgeType::DIRICHLET) {
                    m_nextTemperature[cellIndex(i,j)] = m_parentGrid.getDirichletValue();
                } else {
                    int ghostX = m_totalCellCountX - 1;
                    m_currentTemperature[cellIndex(ghostX, j)] = m_currentTemperature[cellIndex(i,j)];
                    m_nextTemperature[cellIndex(i,j)] = m_currentTemperature[cellIndex(i,j)];
                }
                continue;
            }

            double center = m_currentTemperature[cellIndex(i, j)];
            double north  = m_currentTemperature[cellIndex(i, j + 1)];
            double south  = m_currentTemperature[cellIndex(i, j - 1)];
            double east   = m_currentTemperature[cellIndex(i + 1, j)];
            double west   = m_currentTemperature[cellIndex(i - 1, j)];

            double updated = computeHeatUpdate(
                center, west, east, north, south,
                m_parentGrid.getThermalDiffusivity(),
                m_parentGrid.getTimeStep(),
                m_parentGrid.getCellSizeX()
            );

            m_nextTemperature[cellIndex(i, j)] = updated;
        }
    }

    m_currentTemperature.swap(m_nextTemperature);
}

void SubGrid::exchangeGhostCells() {
    std::vector<MPI_Request> requests;

    std::vector<double> sendRowNorth(m_cellCountX), recvRowNorth(m_cellCountX);
    std::vector<double> sendRowSouth(m_cellCountX), recvRowSouth(m_cellCountX);

    std::vector<double> sendColEast(m_cellCountY), recvColEast(m_cellCountY);
    std::vector<double> sendColWest(m_cellCountY), recvColWest(m_cellCountY);

    for (int i = 1; i <= m_cellCountX; ++i) {
        sendRowNorth[i-1] = m_currentTemperature[cellIndex(i, m_cellCountY)];
        sendRowSouth[i-1] = m_currentTemperature[cellIndex(i, 1)];
    }

    for (int j = 1; j <= m_cellCountY; ++j) {
        sendColEast[j-1] = m_currentTemperature[cellIndex(m_cellCountX, j)];
        sendColWest[j-1] = m_currentTemperature[cellIndex(1, j)];
    }

    // Use m_subGridInfo instead of subGridInfo
    if (m_subGridInfo.hasWestNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Isend(sendColWest.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getWestNeighborRank(), 0, m_subGridInfo.getCommunicator(), &reqSend);
        MPI_Irecv(recvColWest.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getWestNeighborRank(), 1, m_subGridInfo.getCommunicator(), &reqRecv);
        requests.push_back(reqSend);
        requests.push_back(reqRecv);
    }

    if (m_subGridInfo.hasEastNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Isend(sendColEast.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getEastNeighborRank(), 1, m_subGridInfo.getCommunicator(), &reqSend);
        MPI_Irecv(recvColEast.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getEastNeighborRank(), 0, m_subGridInfo.getCommunicator(), &reqRecv);
        requests.push_back(reqSend);
        requests.push_back(reqRecv);
    }

    if (m_subGridInfo.hasNorthNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Isend(sendRowNorth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getNorthNeighborRank(), 2, m_subGridInfo.getCommunicator(), &reqSend);
        MPI_Irecv(recvRowNorth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getNorthNeighborRank(), 3, m_subGridInfo.getCommunicator(), &reqRecv);
        requests.push_back(reqSend);
        requests.push_back(reqRecv);
    }

    if (m_subGridInfo.hasSouthNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Isend(sendRowSouth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getSouthNeighborRank(), 3, m_subGridInfo.getCommunicator(), &reqSend);
        MPI_Irecv(recvRowSouth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getSouthNeighborRank(), 2, m_subGridInfo.getCommunicator(), &reqRecv);
        requests.push_back(reqSend);
        requests.push_back(reqRecv);
    }

    if (!requests.empty()) {
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }

    for (int i = 1; i <= m_cellCountX; ++i) {
        if (m_subGridInfo.hasNorthNeighbor()) {
            m_currentTemperature[cellIndex(i, m_totalCellCountY - 1)] = recvRowNorth[i-1];
        }
        if (m_subGridInfo.hasSouthNeighbor()) {
            m_currentTemperature[cellIndex(i, 0)] = recvRowSouth[i-1];
        }
    }

    for (int j = 1; j <= m_cellCountY; ++j) {
        if (m_subGridInfo.hasEastNeighbor()) {
            m_currentTemperature[cellIndex(m_totalCellCountX - 1, j)] = recvColEast[j-1];
        }
        if (m_subGridInfo.hasWestNeighbor()) {
            m_currentTemperature[cellIndex(0, j)] = recvColWest[j-1];
        }
    }
}

const std::vector<double>& SubGrid::getCurrentTemperature() const {
    return m_currentTemperature;
}

int SubGrid::getCellCountX() const { return m_cellCountX; }
int SubGrid::getCellCountY() const { return m_cellCountY; }
int SubGrid::getTotalCellCountX() const { return m_totalCellCountX; }
int SubGrid::getTotalCellCountY() const { return m_totalCellCountY; }

std::vector<double> SubGrid::getInteriorCells() const {
    std::vector<double> interior(m_cellCountX * m_cellCountY);
    for (int i = 1; i <= m_cellCountX; ++i) {
        for (int j = 1; j <= m_cellCountY; ++j) {
            int localIdx = (i - 1) * m_cellCountY + (j - 1);
            interior[localIdx] = m_currentTemperature[cellIndex(i, j)];
        }
    }
    return interior;
}

} // namespace simulation
