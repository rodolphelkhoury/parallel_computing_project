#ifndef SIMULATION_INCLUDE_GRID_HPP
#define SIMULATION_INCLUDE_GRID_HPP

#include <cstdint>

namespace simulation {

/**
 * @brief Enum for the type of boundary condition on each edge of the grid.
 * 
 * DIRICHLET: fixed value at the boundary.
 * NEUMANN: fixed gradient (flux) at the boundary.
 */
enum class EdgeType : std::uint8_t {
    DIRICHLET,
    NEUMANN
};

/**
 * @brief Class representing a 2D simulation grid with physical and numerical parameters.
 * 
 * Stores the size of the grid, cell spacing, simulation time step, thermal diffusivity,
 * boundary conditions, and optional Dirichlet values.
 */
class Grid {
private:
    int m_totalGridX;
    int m_totalGridY;

    double m_cellSizeX; // dx
    double m_cellSizeY; // dy

    double m_timeStep; // dt
    double m_thermalDiffusivity; // α

    int m_totalIterations;

    EdgeType m_edgeNorth;
    EdgeType m_edgeSouth;
    EdgeType m_edgeEast;
    EdgeType m_edgeWest;

    double m_dirichletValue;

public:
    explicit Grid(int totalX, int totalY,
        double dx, double dy,
        double dt, double alpha,
        int iterations,
        EdgeType north, EdgeType south, EdgeType east, EdgeType west,
        double dirichlet);

    [[nodiscard]] int getTotalGridX() const;
    [[nodiscard]] int getTotalGridY() const;
    [[nodiscard]] double getCellSizeX() const;
    [[nodiscard]] double getCellSizeY() const;
    [[nodiscard]] double getTimeStep() const;
    [[nodiscard]] double getThermalDiffusivity() const;
    [[nodiscard]] int getTotalIterations() const;
    [[nodiscard]] EdgeType getEdgeNorth() const;
    [[nodiscard]] EdgeType getEdgeSouth() const;
    [[nodiscard]] EdgeType getEdgeEast() const;
    [[nodiscard]] EdgeType getEdgeWest() const;
    [[nodiscard]] double getDirichletValue() const;
};

} // namespace simulation

#endif // SIMULATION_INCLUDE_GRID_HPP
