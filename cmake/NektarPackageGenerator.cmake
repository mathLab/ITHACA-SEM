set(CPACK_PACKAGE_VENDOR "Imperial College London")
set(CPACK_PACKAGE_CONTACT 
    "Nektar++ users mailing list <nektar-users@imperial.ac.uk>")

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_PACKAGE_VERSION "${NEKTAR_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${NEKTAR_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${NEKTAR_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${NEKTAR_VERSION_PATCH}")

# Debian-specific options
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_PACKAGE_DEBUG ON)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER 
    "Chris Cantwell <c.cantwell@imperial.ac.uk>")
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://www.nektar.info")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# RPM-specific options
set(CPACK_RPM_PACKAGE_GROUP "Development/Libraries")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DEBUG 0)

# Set up generator-specific logic (i.e. that may change the above depending on
# generator) using additional CPack configuration file.
configure_file(${CMAKE_SOURCE_DIR}/cmake/NektarCPackConfig.cmake.in
               ${CMAKE_BINARY_DIR}/NektarCPackConfig.cmake @ONLY)
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_BINARY_DIR}/NektarCPackConfig.cmake)

# Finally, include the CPack module
include(CPack)
