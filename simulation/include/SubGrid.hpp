#ifndef SIMULATION_INCLUDE_SUBGRID_HPP
#define SIMULATION_INCLUDE_SUBGRID_HPP

#include "SubgridInformation.hpp"
#include "Grid.hpp"

#include <vector>

namespace simulation {

class SubGrid {
public:
    SubGrid(const Grid& parentGrid, const SubGridInformation& subGridInfo);
    void updateCellTemp();
    void exchangeGhostCells();
    void applyBoundaryConditions();

    [[nodiscard]] const std::vector<double>& getCurrentTemperature() const;
    [[nodiscard]] int getCellCountX() const;
    [[nodiscard]] int getCellCountY() const;
    [[nodiscard]] int getTotalCellCountX() const;
    [[nodiscard]] int getTotalCellCountY() const;
    [[nodiscard]] std::vector<double> getInteriorCells() const;

private:
    const Grid& m_parentGrid;
    const SubGridInformation& m_subGridInfo;

    int m_cellCountX;
    int m_cellCountY;

    bool m_hasGhostWest;
    bool m_hasGhostEast;
    bool m_hasGhostNorth;
    bool m_hasGhostSouth;

    int m_totalCellCountX;  // interior + ghostLeft + ghostRight
    int m_totalCellCountY;  // interior + ghostBottom + ghostTop

    std::vector<double> m_currentTemperature;
    std::vector<double> m_nextTemperature;

    inline int cellIndex(int i, int j) const {
        return i * m_totalCellCountY + j;
    }
};

} // namespace simulation

#endif // SIMULATION_INCLUDE_SUBGRID_HPP