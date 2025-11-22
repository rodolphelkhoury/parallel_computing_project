#include "Grid.hpp"

#include <stdexcept>

namespace simulation {

Grid::Grid(int totalX, int totalY,
           double dx, double dy,
           double dt, double alpha,
           int iterations,
           EdgeType north, EdgeType south, EdgeType east, EdgeType west,
           double dirichlet)
    : m_totalGridX(totalX), m_totalGridY(totalY),
      m_cellSizeX(dx), m_cellSizeY(dy),
      m_timeStep(dt), m_thermalDiffusivity(alpha),
      m_totalIterations(iterations),
      m_edgeNorth(north), m_edgeSouth(south), m_edgeEast(east), m_edgeWest(west),
      m_dirichletValue(dirichlet)
{
    // some validation
    if (m_totalGridX <= 0 || m_totalGridY <= 0)
        throw std::invalid_argument("Grid dimensions must be positive.");
    if (m_cellSizeX <= 0 || m_cellSizeY <= 0)
        throw std::invalid_argument("Cell sizes must be positive.");
    if (m_timeStep <= 0)
        throw std::invalid_argument("Time step must be positive.");
}

int Grid::getTotalGridCellsCountX() const { return m_totalGridX; }
int Grid::getTotalGridCellsCountY() const { return m_totalGridY; }
double Grid::getCellSizeX() const { return m_cellSizeX; }
double Grid::getCellSizeY() const { return m_cellSizeY; }
double Grid::getTimeStep() const { return m_timeStep; }
double Grid::getThermalDiffusivity() const { return m_thermalDiffusivity; }
int Grid::getTotalIterations() const { return m_totalIterations; }
EdgeType Grid::getEdgeNorth() const { return m_edgeNorth; }
EdgeType Grid::getEdgeSouth() const { return m_edgeSouth; }
EdgeType Grid::getEdgeEast() const { return m_edgeEast; }
EdgeType Grid::getEdgeWest() const { return m_edgeWest; }
double Grid::getDirichletValue() const { return m_dirichletValue; }

} // namespace simulation
