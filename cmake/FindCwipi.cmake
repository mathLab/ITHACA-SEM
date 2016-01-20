# - Try to find cwipi
# Once done this will define
#  CWIPI_FOUND - System has cwipi
#  CWIPI_INCLUDE_DIRS - The cwipi include directories
#  CWIPI_LIBRARIES - The libraries needed to use cwipi
#  CWIPI_DEFINITIONS - Compiler switches required for using cwipi

set(CWIPI_DEFINITIONS "")

find_path(CWIPI_INCLUDE_DIR cwipi.h
          HINTS ${CWIPI_DIR}/include/
)

find_library(CWIPI_LIBRARY
          NAMES "libcwipi.so"
          HINTS "${CWIPI_DIR}/lib64/" "${CWIPI_DIR}/lib/"
)

find_library(CWIPI_LIBRARY_FVMC
          NAMES "libfvmc.so"
          HINTS "${CWIPI_DIR}/lib64/" "${CWIPI_DIR}/lib/"
)

find_library(CWIPI_LIBRARY_BFTC
          NAMES "libbftc.so"
          HINTS "${CWIPI_DIR}/lib64/" "${CWIPI_DIR}/lib/"
)

set(CWIPI_LIBRARIES
    ${CWIPI_LIBRARY}
    ${CWIPI_LIBRARY_FVMC}
    ${CWIPI_LIBRARY_BFTC}
)
set(CWIPI_INCLUDE_DIRS ${CWIPI_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CWIPI_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(cwipi DEFAULT_MSG
                                  CWIPI_LIBRARY
                                  CWIPI_LIBRARY_FVMC
                                  CWIPI_LIBRARY_BFTC
                                  CWIPI_INCLUDE_DIR
                                 )

mark_as_advanced(CWIPI_LIBRARY CWIPI_INCLUDE_DIR)
