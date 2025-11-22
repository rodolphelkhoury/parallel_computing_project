#include "SubGrid.hpp"
#include "Grid.hpp"
#include <mpi.h>
#include <vector>

namespace simulation {

SubGrid::SubGrid(const Grid& parentGrid, const SubGridInformation& subGridInfo)
    : m_parentGrid(parentGrid)
    , m_subGridInfo(subGridInfo)
{
    m_haloWest  = subGridInfo.hasWestNeighbor();
    m_haloEast  = subGridInfo.hasEastNeighbor();
    m_haloNorth = subGridInfo.hasNorthNeighbor();
    m_haloSouth = subGridInfo.hasSouthNeighbor();

    int totalCellsX = parentGrid.getTotalGridCellsCountX();
    int totalCellsY = parentGrid.getTotalGridCellsCountY();

    int procCountX = subGridInfo.getNumberOfProcessesOnX();
    int procCountY = subGridInfo.getNumberOfProcessesOnY();

    int procX = subGridInfo.getProcessCoordinatesX();
    int procY = subGridInfo.getProcessCoordinatesY();

    int baseCellsX     = totalCellsX / procCountX;
    int leftoverCellsX = totalCellsX % procCountX;

    int baseCellsY     = totalCellsY / procCountY;
    int leftoverCellsY = totalCellsY % procCountY;

    // this process gets 1 extra cell in X if it's in the leftover region
    m_cellCountX = baseCellsX + (procX < leftoverCellsX ? 1 : 0);
    m_cellCountY = baseCellsY + (procY < leftoverCellsY ? 1 : 0);

    m_totalCellCountX = m_cellCountX + static_cast<int>(m_haloWest) + static_cast<int>(m_haloEast);

    m_totalCellCountY = m_cellCountY + static_cast<int>(m_haloSouth) + static_cast<int>(m_haloNorth);

    m_currectTemprature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);
    m_nextTemperature.assign(m_totalCellCountX * m_totalCellCountY, 0.0);

    double dirichletValue = parentGrid.getDirichletValue();
    if (procY == subGridInfo.getProcessGridDimensionsY() - 1 && parentGrid.getEdgeNorth() == EdgeType::DIRICHLET) {
        int j = m_cellCountY;
        for (int i = 1; i <= m_cellCountX; ++i) {
            m_currectTemprature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (procY == 0 && subGridInfo.getEdgeSouth() == EdgeType::DIRICHLET) {
        int j = 1;  // first interior row
        for (int i = 1; i <= m_cellCountX; ++i) {
            m_currectTemprature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (procX == 0 && subGridInfo.getEdgeWest() == EdgeType::DIRICHLET) {
        int i = 1;
        for (int j = 1; j <= m_cellCountY; ++j) {
            m_currectTemprature[cellIndex(i, j)] = dirichletValue;
        }
    }

    if (procX == subGridInfo.getProcessGridDimensionsX() - 1 && subGridInfo.getEdgeEast() == EdgeType::DIRICHLET) {
        int i = m_cellCountX;
        for (int j = 1; j <= m_cellCountY; ++j) {
            m_currectTemprature[cellIndex(i, j)] = dirichletValue;
        }
    }

}

void SubGrid::updatecelltemp() {
    for(size_t i = 1; i <= m_cellCountX; ++i) {
        for(size_t j = 1; j <= m_cellCountY; ++j) {
            // Example update logic (to be replaced with actual computation)
            double center = m_currectTemprature[cellIndex(i, j)];
            double north  = m_currectTemprature[cellIndex(i, j + 1)];
            double south  = m_currectTemprature[cellIndex(i, j - 1)];
            double east   = m_currectTemprature[cellIndex(i + 1, j)];
            double west   = m_currectTemprature[cellIndex(i - 1, j)];
            double cellNextValue=computeHeatUpdate(center, west, east,north, south,m_parentGrid.getThermalDiffusivity(),m_parentGrid.getTimeStep(),m_parentGrid.getCellSizeX());m_nextTemperature[cellIndex(i, j)] = cellNextValue;
        }
    }
    m_currectTemprature.swap(m_nextTemperature);
   
} // namespace simulation