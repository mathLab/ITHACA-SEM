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

SET(SOLVER_UTILS_HEADERS
  CoalescedGeomData.h
  Collection.h 
  CollectionOptimisation.h
  IProduct.h
  Operator.h
)

ADD_NEKTAR_LIBRARY(Collections lib ${NEKTAR_LIBRARY_TYPE}
                    ${COLLECTIONS_SOURCES} ${COLLECTIONS_HEADERS})
TARGET_LINK_LIBRARIES(Collections LINK_PUBLIC LocalRegions)

INSTALL(DIRECTORY ./ 
        DESTINATION ${NEKTAR_INCLUDE_DIR}/Collections 
        COMPONENT dev 
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
