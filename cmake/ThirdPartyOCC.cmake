########################################################################
#
# ThirdParty configuration for Nektar++
#
# OpenCascade
#
########################################################################

IF(NEKTAR_USE_OCC)
    SET(BUILD_OCC ON)

    OPTION(THIRDPARTY_DOWNLOAD_OCC
	"Get OpenCascade from ThirdParty." ${BUILD_OCC})

    IF (THIRDPARTY_DOWNLOAD_OCC)
	IF(WIN32)
	    MESSAGE(SEND_ERROR "Cannot use OpenCascade with Nektar++ on Windows")
	ELSEIF(APPLE)
	    INCLUDE(ExternalProject)
	    EXTERNALPROJECT_ADD(
		opencascade-6.8
		PREFIX ${TPSRC}
		URL ${TPURL}/OCC680osx64.tgz
		URL_MD5 626292523b0691304f0fa271989fbc44
		STAMP_DIR ${TPBUILD}/stamp
		BINARY_DIR ${TPBUILD}/opencascade-6.8
		DOWNLOAD_DIR ${TPSRC}
		SOURCE_DIR ${TPSRC}/opencascade-6.8
		INSTALL_DIR ${TPDIST}
		UPDATE_COMMAND ""
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND cp -a ${TPSRC}/opencascade-6.8/i686/lib/. ${TPDIST}/lib/ COMMAND cp -a ${TPSRC}/opencascade-6.8/i686/inc/. ${TPDIST}/include/
		)
	    LINK_DIRECTORIES(${TPDIST}/lib)
	    INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
	ELSE()
	    INCLUDE(ExternalProject)
	    EXTERNALPROJECT_ADD(
		opencascade-6.8
		PREFIX ${TPSRC}
		URL ${TPURL}/OCC680lin64.tgz
		URL_MD5 d655b6f50998bb9600e081907c247793
		STAMP_DIR ${TPBUILD}/stamp
		DOWNLOAD_DIR ${TPSRC}
		SOURCE_DIR ${TPSRC}/opencascade-6.8
		INSTALL_DIR ${TPDIST}
		UPDATE_COMMAND ""
		CONFIGURE_COMMAND ""
		BUILD_COMMAND ""
		INSTALL_COMMAND cp -a ${TPSRC}/opencascade-6.8/lib/. ${TPDIST}/lib/ COMMAND cp -a ${TPSRC}/opencascade-6.8/inc/. ${TPDIST}/include/
		)
	    LINK_DIRECTORIES(${TPDIST}/lib)
	    INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
	ENDIF()

	SET(OCC_LIBS PTKernel TKernel TKMath TKBRep TKIGES TKSTEP TKSTEPAttr
            TKSTEP209 TKSTEPBase TKShapeSchema TKGeomBase TKGeomAlgo TKG3d TKG2d
            TKXSBase TKPShape TKTopAlgo)
    ENDIF (THIRDPARTY_DOWNLOAD_OCC)
ENDIF(NEKTAR_USE_OCC)
