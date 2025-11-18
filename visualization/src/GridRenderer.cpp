#include "GridRenderer.hpp"
#include "ColorMap.hpp"
#include "Window.hpp"

#include <vector>
#include <GLFW/glfw3.h>

namespace visualization {

GridRenderer::GridRenderer(Window& window) : m_window(window) {}

void GridRenderer::render(const std::vector<double>& grid, int nx, int ny, double minTemp, double maxTemp) {
    GLFWwindow* win = m_window.get();
    if (!win) {
        return;
    }

    glClear(GL_COLOR_BUFFER_BIT);

    int winWidth, winHeight;
    glfwGetFramebufferSize(win, &winWidth, &winHeight);

    float cellW = float(winWidth) / nx;
    float cellH = float(winHeight) / ny;

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {

            double temp = grid[j * nx + i];
            Color c = temperatureToColor(temp, minTemp, maxTemp);

            glColor3f(c.r, c.g, c.b);

            float x0 = i * cellW;
            float y0 = j * cellH;
            float x1 = x0 + cellW;
            float y1 = y0 + cellH;

            glBegin(GL_QUADS);
            glVertex2f(x0, y0);
            glVertex2f(x1, y0);
            glVertex2f(x1, y1);
            glVertex2f(x0, y1);
            glEnd();
        }
    }

    glfwSwapBuffers(win);
    glfwPollEvents();
}

} // namespace visualization
