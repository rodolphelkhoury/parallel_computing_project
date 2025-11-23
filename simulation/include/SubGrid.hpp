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
        // If i/j already are raw coordinates (0..m_totalCellCountX-1 / 0..m_totalCellCountY-1),
        // return direct mapping.
        if (i >= 0 && i < m_totalCellCountX && j >= 0 && j < m_totalCellCountY) {
            return i * m_totalCellCountY + j;
        }

        // Otherwise assume i,j are 1-based interior coordinates (1..m_cellCountX / 1..m_cellCountY)
        // Convert to raw coordinates depending on whether a ghost cell exists on the low side.
        int rawI = (m_hasGhostWest ? i : i - 1);
        int rawJ = (m_hasGhostSouth ? j : j - 1);

        // Safety: clamp to valid range (optional)
        if (rawI < 0) rawI = 0;
        if (rawI >= m_totalCellCountX) rawI = m_totalCellCountX - 1;
        if (rawJ < 0) rawJ = 0;
        if (rawJ >= m_totalCellCountY) rawJ = m_totalCellCountY - 1;

        return rawI * m_totalCellCountY + rawJ;
    }

};

} // namespace simulation

#endif // SIMULATION_INCLUDE_SUBGRID_HPP
