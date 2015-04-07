########################################################################
#
# ThirdParty configuration for Nektar++
#
# HDF5
#
########################################################################

OPTION(NEKTAR_USE_HDF5
    "Enable HDF5 I/O support." OFF)

IF (NEKTAR_USE_HDF5)
  INCLUDE(FindHDF5)

  INCLUDE_DIRECTORIES(SYSTEM ${HDF5_INCLUDE_DIRS})
ENDIF (NEKTAR_USE_HDF5)