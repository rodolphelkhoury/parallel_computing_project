#ifndef VISUALIZATION_INCLUDE_COLOR_MAP_HPP
#define VISUALIZATION_INCLUDE_COLOR_MAP_HPP

namespace visualization {

/**
 * @brief Utilities for mapping scalar temperature values to RGB colors.
 *
 * Provides a simple `Color` struct holding RGB components and a
 * `temperatureToColor` function that converts a temperature value in the
 * range `[minValue, maxValue]` into an RGB color suitable for visualization.
 */
struct Color {
    float r;
    float g;
    float b;
};

/**
 * @brief Convert a temperature value to an RGB color.
 *
 * @param value The temperature value to map.
 * @param minValue Minimum expected temperature.
 * @param maxValue Maximum expected temperature.
 * @return Color RGB color corresponding to `value`.
 */
Color temperatureToColor(double value, double minValue, double maxValue);

} // namespace visualization

#endif // VISUALIZATION_INCLUDE_COLOR_MAP_HPP
