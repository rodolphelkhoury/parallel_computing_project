#ifndef SIMULATION_INCLUDE_HEATEQUATION_HPP
#define SIMULATION_INCLUDE_HEATEQUATION_HPP

namespace simulation {

[[nodiscard]] inline double computeHeatUpdate(
    double centerValue,
    double leftValue,
    double rightValue,
    double topValue,
    double bottomValue,
    double thermalDiffusivity,
    double timeStep,
    double cellSizeX
) {
    const double laplacian = leftValue + rightValue + topValue + bottomValue - 4.0 * centerValue;
    return centerValue + (thermalDiffusivity * timeStep / (cellSizeX * cellSizeX)) * laplacian;
}

} // namespace simulation

#endif // SIMULATION_INCLUDE_HEATEQUATION_HPP
