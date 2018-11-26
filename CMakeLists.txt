cmake_minimum_required(VERSION 3.1)
project(mcf_srep)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)

find_package(LIBIGL REQUIRED QUIET)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Add your project files
file(GLOB SRCFILES *.cxx)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw ${VTK_LIBRARIES})
