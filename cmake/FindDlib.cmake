# Try to find the Dlib library (http://dlib.net/)
# Once done this will define
#
#  DLIB_FOUND - system has Dlib
#  DLIB_INCLUDE_DIR - the Dlib include directory
#  DLIB_LIBRARY - The Dlib library

set(DLIB_INCLUDE_DIR "" CACHE PATH "Path for dlib include files")
if(NOT DLIB_FOUND)

  FIND_PATH(DLIB_INCLUDE_DIR optimization.h
    ${DLIB_INCLUDE_DIR}
    $ENV{DLIB}/
    $ENV{DLIBROOT}/
    $ENV{DLIB_ROOT}/
    $ENV{DLIB_DIR}/
    /usr/include
    /usr/local/include
    )

  if(DLIB_INCLUDE_DIR)
    set(DLIB_FOUND TRUE)
    set(DLIB_LIBRARY dlib)
    add_subdirectory(${DLIB_INCLUDE_DIR} "dlib")
  endif(DLIB_INCLUDE_DIR)

endif()
