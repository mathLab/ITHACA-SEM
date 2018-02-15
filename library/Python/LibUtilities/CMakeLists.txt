SET(LIBUTILITIES_SOURCES
    LibUtilities.cpp
    BasicUtils/SessionReader.cpp
    BasicUtils/SharedArray.cpp
    BasicUtils/ShapeType.cpp
    Foundations/Basis.cpp
    Foundations/Points.cpp
    LinearAlgebra/NekMatrix.cpp
    )

ADD_NEKPY_LIBRARY(LibUtilities LIBUTILITIES_SOURCES)
