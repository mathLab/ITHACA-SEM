OPTION(NEKTAR_USE_OCC "Use opencascade for geometry interface." OFF)

IF(NEKTAR_USE_OCC)

	SET(BUILD_OCC ON)

	OPTION(THIRDPARTY_DOWNLOAD_OCC
	    "Get OpenCascade from ThirdParty." ${BUILD_OCC})

	IF (THIRDPARTY_DOWNLOAD_OCC)
	    INCLUDE(ExternalProject)
	    EXTERNALPROJECT_ADD(
	        opencascade-6.8
	        PREFIX ${TPSRC}
	        URL http://ae-nektar.ae.ic.ac.uk/~mt4313/OCC680osx64.tgz
	        URL_MD5 626292523b0691304f0fa271989fbc44
	        STAMP_DIR ${TPBUILD}/stamp
	        DOWNLOAD_DIR ${TPSRC}
	        SOURCE_DIR ${TPSRC}/opencascade-6.8
	        INSTALL_DIR ${TPDIST}
			UPDATE_COMMAND ""
			CONFIGURE_COMMAND ""
	        BUILD_COMMAND ""
			INSTALL_COMMAND ""
	    )
	SET(OCC_LIBS PTKernel TKernel TKMath TKBRep TKIGES TKSTEP TKSTEPAttr TKSTEP209 TKSTEPBase TKShapeSchema TKGeomBase TKGeomAlgo TKG3d TKG2d TKXSBase TKPShape TKTopAlgo)
	LINK_DIRECTORIES(${TPSRC}/opencascade-6.8/i686/lib)
	INCLUDE_DIRECTORIES(SYSTEM ${TPSRC}/opencascade-6.8/i686/inc)

	ENDIF (THIRDPARTY_DOWNLOAD_OCC)



ENDIF(NEKTAR_USE_OCC)
