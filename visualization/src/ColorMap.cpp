#include "ColorMap.hpp"

#include <algorithm>

namespace visualization {

Color temperatureToColor(double temp, double minValue, double maxValue) {
    double t = (temp - minValue) / (maxValue - minValue);
    t = std::clamp(t, 0.0, 1.0); // keep value between 0 and 1

    return Color{
        static_cast<float>(t),         // red increases with temperature
        static_cast<float>(0.5 * t),   // green increases slowly
        static_cast<float>(1.0 - t)    // blue decreases as temp increases
    };
}

} // namespace visualization
