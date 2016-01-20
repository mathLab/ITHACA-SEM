# - Try to find openPALM
# Once done this will define
#  OPENPALM_FOUND - System has openPALM
#  OPENPALM_INCLUDE_DIRS - The openPALM include directories
#  OPENPALM_LIBRARIES - The libraries needed to use openPALM
#  OPENPALM_DEFINITIONS - Compiler switches required for using openPALM

set(OPENPALM_DEFINITIONS "")

find_path(OPENPALM_INCLUDE_DIR palmlib.h
          HINTS ${OPENPALM_DIR}/include/
)

find_library(OPENPALM_LIBRARY
          NAMES "libpalm.a"
          HINTS "${OPENPALM_DIR}/lib64/" "${OPENPALM_DIR}/lib/"
)

set(OPENPALM_LIBRARIES
    ${OPENPALM_LIBRARY}
)
set(OPENPALM_INCLUDE_DIRS ${OPENPALM_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set OPENPALM_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(openPALM DEFAULT_MSG
                                  OPENPALM_LIBRARY
                                  OPENPALM_INCLUDE_DIR
                                 )

mark_as_advanced(OPENPALM_LIBRARY OPENPALM_INCLUDE_DIR)
