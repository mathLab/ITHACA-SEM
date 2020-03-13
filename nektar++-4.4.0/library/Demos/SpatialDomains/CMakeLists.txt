SET(LinkLibraries 
    optimized SpatialDomains debug SpatialDomains-g
    optimized StdRegions debug StdRegions-g
    optimized LibUtilities debug LibUtilities-g
    optimized ${Boost_THREAD_LIBRARY_RELEASE} debug ${Boost_THREAD_LIBRARY_DEBUG}
    optimized ${Boost_IOSTREAMS_LIBRARY_RELEASE} debug ${Boost_IOSTREAMS_LIBRARY_DEBUG}
    optimized ${Boost_ZLIB_LIBRARY_RELEASE} debug ${Boost_ZLIB_LIBRARY_DEBUG}
    optimized ${TINYXML_LIBRARY} debug ${TINYXML_LIBRARY}
)

#SET(SpatialDomains1DSource  Graph1D.cpp)
#ADD_NEKTAR_EXECUTABLE(Graph1D SpatialDomains1DSource)
#TARGET_LINK_LIBRARIES(Graph1D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(Graph1D)

#SET(SpatialDomains2DSource  Graph2D.cpp)
#ADD_NEKTAR_EXECUTABLE(Graph2D SpatialDomains2DSource)
#TARGET_LINK_LIBRARIES(Graph2D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(Graph2D)

#SET(SpatialDomains3DSource  Graph3D.cpp)
#ADD_NEKTAR_EXECUTABLE(Graph3D SpatialDomains3DSource)
#TARGET_LINK_LIBRARIES(Graph3D ${LinkLibraries})
#SET_LAPACK_LINK_LIBRARIES(Graph3D)

# Generate list of available subdirectories
FILE(GLOB dir_list "*")
FOREACH(dir ${dir_list})
    IF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
        ADD_SUBDIRECTORY(${dir})
    ENDIF(IS_DIRECTORY ${dir} AND EXISTS ${dir}/CMakeLists.txt)
ENDFOREACH(dir ${dir_list})
