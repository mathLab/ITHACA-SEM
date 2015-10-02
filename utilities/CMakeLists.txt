MACRO(ADD_UTILITIES_EXECUTABLE name component source)

    ADD_EXECUTABLE(${name} ${source} ${ARGN})

    SET_COMMON_PROPERTIES(${name})

    TARGET_LINK_LIBRARIES(${name} SolverUtils)
    
    IF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
        SET_TARGET_PROPERTIES(${name} 
            PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/${LIB_DIR} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
    ENDIF( ${CMAKE_SYSTEM} MATCHES "Darwin-*")
	SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${component})
    INSTALL(TARGETS ${name} 
		RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${component})

ENDMACRO(ADD_UTILITIES_EXECUTABLE name component source)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library ${CMAKE_SOURCE_DIR}/utilities)

ADD_UTILITIES_EXECUTABLE(FldAddFld         util FldAddFld.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddScalGrad    util FldAddScalGrad.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddScalGrad_elmt util FldAddScalGrad_elmt.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddWSS         util FldAddWSS.cpp)
ADD_UTILITIES_EXECUTABLE(FldAddMultiShear  util FldAddMultiShear.cpp)
ADD_UTILITIES_EXECUTABLE(ProbeFld          util-extra ProbeFld.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToTecplot      util XmlToTecplot.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToTecplotWireFrame  util XmlToTecplotWireFrame.cpp)
ADD_UTILITIES_EXECUTABLE(XmlToVtk          util XmlToVtk.cpp)

#SUBDIRS(PostProcessing FieldConvert MeshConvert SimpleMeshGen Extras)

# Generate list of available utilities
FILE(GLOB dir_list "*")
LIST(REMOVE_ITEM dir_list "${CMAKE_SOURCE_DIR}/.svn")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})

