# Try to find the meshfix library
# Once done this will define
#
#  MESHFIX_FOUND - system has meshfix
#  MESHFIX_INCLUDE_DIRS - **the** meshfix include directory
#  MESHFIX_LIBRARIES - The meshfix library

if(NOT MESHFIX_FOUND)

  FIND_PATH(MESHFIX_DIR meshfix.h
     ${PROJECT_SOURCE_DIR}/../../../meshfix/
     ${PROJECT_SOURCE_DIR}/../../meshfix/
     ${PROJECT_SOURCE_DIR}/../meshfix/
     ${PROJECT_SOURCE_DIR}/meshfix/
     $ENV{MESHFIX}/
     $ENV{MESHFIXROOT}/
     $ENV{MESHFIX_ROOT}/
     $ENV{MESHFIX_DIR}/
     /usr/include
     /usr/local/include
  )

  if(MESHFIX_DIR)
    set(MESHFIX_FOUND TRUE)
    set(JMESHEXT_DIR "${MESHFIX_DIR}/JMeshExt-1.0alpha_src")
    set(JMESH_DIR "${JMESHEXT_DIR}/JMeshLib-1.2/")
    set(MESHFIX_INCLUDE_DIRS
      ${JMESHEXT_DIR}/include
      ${JMESH_DIR}/include
      ${MESHFIX_DIR})
    set(MESHFIX_LIBRARIES 
      meshfix 
      jmesh 
      jmeshext 
      nl 
      superlu 
      ${LAPACK_LIBRARIES} 
      ${BLAS_LIBRARIES})
    #set(MESHFIX_LIBRARIES meshfix jmesh)
    add_definitions("-DMESHFIX_WITH_LIBIGL")
    add_subdirectory(${MESHFIX_DIR} "meshfix")
    link_directories(meshfix)
  endif(MESHFIX_DIR)

endif()
