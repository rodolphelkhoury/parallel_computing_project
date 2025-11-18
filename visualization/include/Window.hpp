#ifndef VISUALIZATION_INCLUDE_WINDOW_HPP
#define VISUALIZATION_INCLUDE_WINDOW_HPP

#include <GLFW/glfw3.h>
#include <string>

namespace visualization {

/**
 * @brief Class representing a window for OpenGL/GLFW-based visualization.
 *
 * Manages a GLFW window used to render visualization output. Stores the
 * underlying `GLFWwindow*` and provides simple accessors to check
 * validity and retrieve the native window handle.
 */
class Window {
private:
    GLFWwindow* m_window = nullptr;

public:
    Window(int width, int height, const std::string& title);
    ~Window();

    bool isValid() const;
    GLFWwindow* get() const;
};

} // namespace visualization

#endif // VISUALIZATION_INCLUDE_WINDOW_HPP