SET(COLLECTIONS_SOURCES
  CoalescedGeomData.cpp
  Collection.cpp
  CollectionOptimisation.cpp
  Operator.cpp
  BwdTrans.cpp
  IProductWRTBase.cpp
  PhysDeriv.cpp
  IProductWRTDerivBase.cpp
  IProduct.cpp
  )

SET(COLLECTIONS_HEADERS
  CoalescedGeomData.h
  Collection.h 
  CollectionOptimisation.h
  IProduct.h
  Operator.h
)

ADD_NEKTAR_LIBRARY(NAME Collections
    DESCRIPTION "Nektar++ Collections library"
    SOURCES ${COLLECTIONS_SOURCES}
    HEADERS ${COLLECTIONS_HEADERS}
    DEPENDS LocalRegions)

INSTALL(DIRECTORY ./ 
        DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections 
        COMPONENT dev 
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
