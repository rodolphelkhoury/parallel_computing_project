#ifndef SIMULATION_INCLUDE_SubGrid_HPP
#define SIMULATION_INCLUDE_SubGrid_HPP

#include "SubgridInformation.hpp"
#include "Grid.hpp"

#include <vector>

namespace simulation {

class SubGrid {
public:
    SubGrid(const Grid& parentGrid, const SubGridInformation& subGridInfo);

private:
    const Grid& m_parentGrid;
    const SubGridInformation& m_subGridInfo;

    int m_cellCountX;
    int m_cellCountY;

    bool m_haloWest;
    bool m_haloEast;
    bool m_haloNorth;
    bool m_haloSouth;

    int m_totalCellCountX;  // interior + haloLeft + haloRight
    int m_totalCellCountY;  // interior + haloBottom + haloTop

    std::vector<double> m_currentTemperature;
    std::vector<double> m_nextTemperature;

    std::vector<double> m_sendBufferLeft;
    std::vector<double> m_sendBufferRight;
    std::vector<double> m_recvBufferLeft;
    std::vector<double> m_recvBufferRight;

    inline int cellIndex(int i, int j) const {
        return i * m_totalCellCountY + j;
    }
};

} // namespace simulation

#endif // SIMULATION_INCLUDE_SubGrid_HPP