cmake_minimum_required(VERSION 2.6)
project(Tiling)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(libigl REQUIRED)

if (NOT LIBIGL_FOUND)
   message(FATAL_ERROR "libigl not found --- You can download it using: \n git clone --recursive https://github.com/libigl/libigl.git ${PROJECT_SOURCE_DIR}/../libigl")
endif()

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") #### Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltet -DTETLIBRARY")

### libIGL options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

### add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)

### libIGL options: choose your dependencies (by default everything is OFF except opengl) ###
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL_GLFW "Needed to use viewer" ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_GLFW        "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            OFF)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         ON)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       ON)
option(LIBIGL_WITH_XML         "Use XML"            ON)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_BOOLEAN     "Use Cork boolean"   OFF) # This option is not supported yet
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           ON)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

### Adding libIGL: choose the path to your local copy libIGL ###
### This is going to compile everything you requested ###
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
set(LIBIGL_CMAKE "${LIBIGL_INCLUDE_DIR}/../shared/cmake" CACHE PATH "Path for libigl cmake folder")
add_subdirectory("${LIBIGL_CMAKE}" "libigl")

### ligIGL information ###
message("libIGL includes: ${LIBIGL_INCLUDE_DIRS}")
message("libIGL libraries: ${LIBIGL_LIBRARIES}")
message("libIGL extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libIGL extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libIGL definitions: ${LIBIGL_DEFINITIONS}")

# # Also need ibtourtre
# find_package(libtourtre REQUIRED)
# 
# if (NOT LIBTOURTRE_FOUND)
#   message(FATAL_ERROR "libtourtre not found -- You can download it using: \n git clone https://github.com/sedillard/libtourtre.git")
# endif(NOT LIBTOURTRE_FOUND)
# message("libtourtre includes: ${LIBTOURTRE_INCLUDE_DIR}")
# message("libtourtre lib: ${LIBTOURTRE_LIB}")

### Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS} ${LIBTOURTRE_INCLUDE_DIR})
add_definitions(${LIBIGL_DEFINITIONS})

### Add your project files
SET(BINARY_FILES "${CMAKE_CURRENT_LIST_DIR}/main.cpp;${CMAKE_CURRENT_LIST_DIR}/test_harmonic.cpp;${CMAKE_CURRENT_LIST_DIR}/contourtree.cpp;${CMAKE_CURRENT_LIST_DIR}/tiler_main.cpp;${CMAKE_CURRENT_LIST_DIR}/debug/debugger_main.cpp;${CMAKE_CURRENT_LIST_DIR}/icoso_tile.cpp")
message("Removing binary files ${BINARY_FILES}")
FILE(GLOB SRCFILES *.cpp)
LIST(REMOVE_ITEM SRCFILES ${BINARY_FILES})
message("Project sources are: ${SRCFILES}")
add_executable(${PROJECT_NAME}_bin main.cpp ${SRCFILES} ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBTOURTRE_LIB} ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES} ${LIBTOURTRE_LIB})

add_executable(test_harmonic_bin test_harmonic.cpp 
  curvatureFlow.cpp 
  Helpers.cpp 
  ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(test_harmonic_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

add_executable(test_contour_tree contourtree.cpp TilingUtils.cpp ${LIBIGL_EXTRA_SOURCES} offsetSurface.cpp)
target_link_libraries(test_contour_tree ${LIBTOURTRE_LIB} ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

add_executable(test_icoso_tile icoso_tile.cpp  ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(test_icoso_tile ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

# Dynamic programming for correct topologies.
add_executable(tiler_main tiler_main.cpp Tiler.cpp ${SRCFILES})
target_link_libraries(tiler_main ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})

# Debugging runner for tiler_main (enter in tile stacks and check if valid).
# From the build directory -
#   $ make -j3 tiler_debugger && ./tiler_debugger < ../debug/tiler_debugger.in
add_executable(tiler_debugger debug/tiler_debugger.cpp ${SRCFILES})
target_link_libraries(tiler_debugger ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES} ${LIBTOURTRE_LIB})
