
SET(RectangularMeshSources 
    RectangularMesh.cpp)

SET(Regular2DMeshSources
    Regular2DMeshGenerator.cpp)


ADD_UTILITIES_EXECUTABLE(RectangularMesh util-extra ${RectangularMeshSources})
ADD_UTILITIES_EXECUTABLE(Regular2DMeshGenerator util-extra ${Regular2DMeshSources})

