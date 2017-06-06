project(Remeshing)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

######################################################################
# OpenMesh
######################################################################
find_package(OpenMesh REQUIRED)
if (NOT OPENMESH_FOUND)
  message(FATAL_ERROR "OpenMesh not found --- You can download it from http://www.openmesh.org/download/ . Please not that only version 3.1 has been tested.")
endif(NOT OPENMESH_FOUND)
add_definitions(-DOM_STATIC_BUILD)
######################################################################

find_package(Eigen REQUIRED)
if (NOT EIGEN_FOUND)
  message(FATAL_ERROR "Eigen not found!")
endif(NOT EIGEN_FOUND)

SET(SRC_FILES
  IsotropicRemeshing.cpp
  NearestTriangleSearch.cpp
  OpenMeshConversion.cpp
  Remesh.cpp
  SquareDistancePointToTriangle.cpp
  TriangleSurfaceArea.cpp
  )
SET(INCLUDE_FILES
  IsotropicRemeshing.h
  KDTreeNode.h
  MeshTraits.h
  NearestTriangleSearch.h
  OpenMeshConversion.h
  OpenMeshGlobals.h
  Remesh.h
  )

INCLUDE_DIRECTORIES(${OPENMESH_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
add_library(Remeshing ${SRC_FILES})
target_link_libraries(Remeshing ${OPENMESH_LIBRARIES})
