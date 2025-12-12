find_package(OpenGL REQUIRED)

add_library(visualization_module
    src/Window.cpp
    src/ColorMap.cpp
    src/GridRenderer.cpp
)

target_include_directories(visualization_module
    PUBLIC include
)

target_link_libraries(visualization_module
    PUBLIC simulation_module
    PRIVATE ${OPENGL_LIBRARIES}
    PRIVATE glfw
)
