#include "HeatEquation.hpp"
#include "SubGrid.hpp"
#include "Grid.hpp"

#include <mpi.h>
#include <vector>
#include <algorithm>

namespace simulation {

SubGrid::SubGrid(const Grid& parentGrid, const SubGridInformation& info)
    : m_parentGrid(parentGrid),
      m_subGridInfo(info)
{
    m_hasGhostWest  = info.hasWestNeighbor();
    m_hasGhostEast  = info.hasEastNeighbor();
    m_hasGhostNorth = info.hasNorthNeighbor();
    m_hasGhostSouth = info.hasSouthNeighbor();

    const int totalX = parentGrid.getTotalGridCellsCountX();
    const int totalY = parentGrid.getTotalGridCellsCountY();

    const int numberOfProcessesSuivantX = info.getNumberOfProcessesOnX();
    const int numberOfProcessesSuivantY = info.getNumberOfProcessesOnY();
    const int processCoordinatesSuivantX = info.getProcessCoordinatesX();
    const int processCoordiantesSuivantY = info.getProcessCoordinatesY();

    const int baseX = totalX / numberOfProcessesSuivantX;
    const int remainingProcessesSuivantX  = totalX % numberOfProcessesSuivantX;
    const int baseY = totalY / numberOfProcessesSuivantY;
    const int remainingProcessesSuivantY  = totalY % numberOfProcessesSuivantY;

    m_cellCountX = baseX + (processCoordinatesSuivantX < remainingProcessesSuivantX);
    m_cellCountY = baseY + (processCoordiantesSuivantY < remainingProcessesSuivantY);

    m_totalCellCountX = m_cellCountX + 2;
    m_totalCellCountY = m_cellCountY + 2;

    m_currentTemperature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);
    m_nextTemperature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);
}

void SubGrid::exchangeGhostCells() {
    std::vector<MPI_Request> requests;

    int startI = 1;
    int startJ = 1;

    std::vector<double> sendColWest(m_cellCountY), recvColWest(m_cellCountY);
    std::vector<double> sendColEast(m_cellCountY), recvColEast(m_cellCountY);
    std::vector<double> sendRowSouth(m_cellCountX), recvRowSouth(m_cellCountX);
    std::vector<double> sendRowNorth(m_cellCountX), recvRowNorth(m_cellCountX);

    for (int j = 0; j < m_cellCountY; ++j) {
        sendColWest[j] = m_currentTemperature[cellIndex(startI, startJ + j)];
        sendColEast[j] = m_currentTemperature[cellIndex(startI + m_cellCountX - 1, startJ + j)];
    }

    for (int i = 0; i < m_cellCountX; ++i) {
        sendRowSouth[i] = m_currentTemperature[cellIndex(startI + i, startJ)];
        sendRowNorth[i] = m_currentTemperature[cellIndex(startI + i, startJ + m_cellCountY - 1)];
    }

    if (m_subGridInfo.hasWestNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Irecv(recvColWest.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getWestNeighborRank(), 1, m_subGridInfo.getCommunicator(), &reqRecv);
        MPI_Isend(sendColWest.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getWestNeighborRank(), 0, m_subGridInfo.getCommunicator(), &reqSend);
        requests.push_back(reqRecv);
        requests.push_back(reqSend);
    }

    if (m_subGridInfo.hasEastNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Irecv(recvColEast.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getEastNeighborRank(), 0, m_subGridInfo.getCommunicator(), &reqRecv);
        MPI_Isend(sendColEast.data(), m_cellCountY, MPI_DOUBLE, m_subGridInfo.getEastNeighborRank(), 1, m_subGridInfo.getCommunicator(), &reqSend);
        requests.push_back(reqRecv);
        requests.push_back(reqSend);
    }

    if (m_subGridInfo.hasSouthNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Irecv(recvRowSouth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getSouthNeighborRank(), 2, m_subGridInfo.getCommunicator(), &reqRecv);
        MPI_Isend(sendRowSouth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getSouthNeighborRank(), 3, m_subGridInfo.getCommunicator(), &reqSend);
        requests.push_back(reqRecv);
        requests.push_back(reqSend);
    }

    if (m_subGridInfo.hasNorthNeighbor()) {
        MPI_Request reqSend, reqRecv;
        MPI_Irecv(recvRowNorth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getNorthNeighborRank(), 3, m_subGridInfo.getCommunicator(), &reqRecv);
        MPI_Isend(sendRowNorth.data(), m_cellCountX, MPI_DOUBLE, m_subGridInfo.getNorthNeighborRank(), 2, m_subGridInfo.getCommunicator(), &reqSend);
        requests.push_back(reqRecv);
        requests.push_back(reqSend);
    }

    if (!requests.empty()) {
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }

    if (m_subGridInfo.hasWestNeighbor()) {
        for (int j = 0; j < m_cellCountY; ++j) {
            m_currentTemperature[cellIndex(0, startJ + j)] = recvColWest[j];
        }
    }

    if (m_subGridInfo.hasEastNeighbor()) {
        for (int j = 0; j < m_cellCountY; ++j) {
            m_currentTemperature[cellIndex(m_totalCellCountX - 1, startJ + j)] = recvColEast[j];
        }
    }

    if (m_subGridInfo.hasSouthNeighbor()) {
        for (int i = 0; i < m_cellCountX; ++i) {
            m_currentTemperature[cellIndex(startI + i, 0)] = recvRowSouth[i];
        }
    }

    if (m_subGridInfo.hasNorthNeighbor()) {
        for (int i = 0; i < m_cellCountX; ++i) {
            m_currentTemperature[cellIndex(startI + i, m_totalCellCountY - 1)] = recvRowNorth[i];
        }
    }
}


void SubGrid::updateCellTemp() {
    const double diffusivity = m_parentGrid.getThermalDiffusivity();
    const double dt = m_parentGrid.getTimeStep();
    const double dx = m_parentGrid.getCellSizeX();

    for (int i = 1; i <= m_cellCountX; ++i) {
        for (int j = 1; j <= m_cellCountY; ++j) {

            const double center = m_currentTemperature[cellIndex(i, j)];
            const double west = m_currentTemperature[cellIndex(i - 1, j)];
            const double east = m_currentTemperature[cellIndex(i + 1, j)];
            const double south = m_currentTemperature[cellIndex(i, j - 1)];
            const double north = m_currentTemperature[cellIndex(i, j + 1)];

            m_nextTemperature[cellIndex(i, j)] = computeHeatUpdate(
                center, west, east, north, south,
                diffusivity, dt, dx
            );
        }
    }

    m_currentTemperature.swap(m_nextTemperature);
}

void SubGrid::applyBoundaryConditions() {
    const int procX = m_subGridInfo.getProcessCoordinatesX();
    const int procY = m_subGridInfo.getProcessCoordinatesY();
    const int procCountX = m_subGridInfo.getNumberOfProcessesOnX();
    const int procCountY = m_subGridInfo.getNumberOfProcessesOnY();

    auto& temp = m_currentTemperature;
    const double dirichletValue = m_parentGrid.getDirichletValue();

    if (procY == procCountY - 1) {
        int j = 1 + m_cellCountY - 1;
        for (int i = 1; i <= m_cellCountX; ++i) {
            if (m_parentGrid.getEdgeNorth() == EdgeType::DIRICHLET) {
                temp[cellIndex(i, j)] = dirichletValue;
            } else {
                temp[cellIndex(i, m_totalCellCountY - 1)] = temp[cellIndex(i, j)];
            }
        }
    }

    if (procY == 0) {
        int j = 1;
        for (int i = 1; i <= m_cellCountX; ++i) {
            if (m_parentGrid.getEdgeSouth() == EdgeType::DIRICHLET) {
                temp[cellIndex(i, j)] = dirichletValue;
            } else {
                temp[cellIndex(i, 0)] = temp[cellIndex(i, j)];
            }
        }
    }

    if (procX == 0) {
        int i = 1;
        for (int j = 1; j <= m_cellCountY; ++j) {
            if (m_parentGrid.getEdgeWest() == EdgeType::DIRICHLET) {
                temp[cellIndex(i, j)] = dirichletValue;
            } else {
                temp[cellIndex(0, j)] = temp[cellIndex(i, j)];
            }
        }
    }

    if (procX == procCountX - 1) {
        int i = 1 + m_cellCountX - 1;
        for (int j = 1; j <= m_cellCountY; ++j) {
            if (m_parentGrid.getEdgeEast() == EdgeType::DIRICHLET) {
                temp[cellIndex(i, j)] = dirichletValue;
            } else {
                temp[cellIndex(m_totalCellCountX - 1, j)] = temp[cellIndex(i, j)];
            }
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
    std::vector<double> interior;
    interior.reserve(m_cellCountX * m_cellCountY);

    for (int i = 1; i <= m_cellCountX; ++i) {
        for (int j = 1; j <= m_cellCountY; ++j) {
            interior.push_back(m_currentTemperature[cellIndex(i, j)]);
        }
    }

    return interior;
}

} // namespace simulation