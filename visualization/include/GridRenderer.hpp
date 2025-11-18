#ifndef VISUALIZATION_INCLUDE_GRID_RENDERER_HPP
#define VISUALIZATION_INCLUDE_GRID_RENDERER_HPP

#include "Window.hpp"
#include "ColorMap.hpp"

#include <vector>

namespace visualization {

/**
 * @brief Renderer for drawing a simulation grid into a `Window`.
 *
 * Converts scalar grid data into colors using the colormap utilities and
 * renders the resulting image into the associated OpenGL/GLFW `Window`.
 * The `render` method accepts a flat vector of temperature values along
 * with grid dimensions and min/max temperature for color mapping.
 */
class GridRenderer {
private:
    Window& m_window;

public:
    GridRenderer(Window& window);

    void render(const std::vector<double>& data, int nx, int ny, double minTemp, double maxTemp);
};

} // namespace visualization

#endif // VISUALIZATION_INCLUDE_GRID_RENDERER_HPP
