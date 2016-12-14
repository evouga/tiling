# Try to find the LIBTOURTRE library
# Once done this will define
#
#  LIBTOURTRE_FOUND - system has LIBTOURTRE
#  LIBTOURTRE_INCLUDE_DIR - **the** LIBTOURTRE include directory
#  LIBTOURTRE_LIB - the LIBTOURTRE library

if(NOT LIBTOURTRE_FOUND)

  FIND_PATH(LIBTOURTRE_INCLUDE_DIR tourtre.h
     ${PROJECT_SOURCE_DIR}/../../../include
     ${PROJECT_SOURCE_DIR}/../../include
     ${PROJECT_SOURCE_DIR}/../include
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/../../../libtourtre/include
     ${PROJECT_SOURCE_DIR}/../../libtourtre/include
     ${PROJECT_SOURCE_DIR}/../libtourtre/include
     ${PROJECT_SOURCE_DIR}/libtourtre/include
     $ENV{LIBTOURTRE}/include
     $ENV{LIBTOURTREROOT}/include
     $ENV{LIBTOURTRE_ROOT}/include
     $ENV{LIBTOURTRE_DIR}/include
     /usr/include
     /usr/local/include
  )

  FIND_LIBRARY(LIBTOURTRE_LIB
    NAMES tourtre libtourtre
    HINTS ${LIBTOURTRE_INCLUDE_DIR}/../
  )

  message("include? ${LIBTOURTRE_INCLUDE_DIR} lib? ${LIBTOURTRE_LIB}")
  if(LIBTOURTRE_INCLUDE_DIR AND LIBTOURTRE_LIB)
    set(LIBTOURTRE_FOUND TRUE)
    set(LIBTOURTRE_INCLUDE_DIRS ${LIBTOURTRE_INCLUDE_DIR})
  endif(LIBTOURTRE_INCLUDE_DIR AND LIBTOURTRE_LIB)

endif()
