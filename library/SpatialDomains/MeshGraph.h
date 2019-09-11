////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH_H

#include <unordered_map>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>

#include <SpatialDomains/MeshEntities.hpp>
#include <SpatialDomains/HexGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/TriGeom.h>

#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlDocument;

namespace Nektar
{
namespace SpatialDomains
{
typedef std::map<int, std::pair<LibUtilities::ShapeType, std::vector<int>>>
CompositeDescriptor;

enum ExpansionType
{
    eNoExpansionType,
    eModified,
    eModifiedQuadPlus1,
    eModifiedQuadPlus2,
    eModifiedGLLRadau10,
    eOrthogonal,
    eGLL_Lagrange,
    eGLL_Lagrange_SEM,
    eGauss_Lagrange,
    eGauss_Lagrange_SEM,
    eFourier,
    eFourierSingleMode,
    eFourierHalfModeRe,
    eFourierHalfModeIm,
    eChebyshev,
    eFourierChebyshev,
    eChebyshevFourier,
    eFourierModified,
    eExpansionTypeSize
};

// Keep this consistent with the enums in ExpansionType.
// This is used in the BC file to specify the expansion type.
const std::string kExpansionTypeStr[] = {"NOTYPE",
                                         "MODIFIED",
                                         "MODIFIEDQUADPLUS1",
                                         "MODIFIEDQUADPLUS2",
                                         "MODIFIEDGLLRADAU10",
                                         "ORTHOGONAL",
                                         "GLL_LAGRANGE",
                                         "GLL_LAGRANGE_SEM",
                                         "GAUSS_LAGRANGE",
                                         "GAUSS_LAGRANGE_SEM",
                                         "FOURIER",
                                         "FOURIERSINGLEMODE",
                                         "FOURIERHALFMODERE",
                                         "FOURIERHALFMODEIM",
                                         "CHEBYSHEV",
                                         "FOURIER-CHEBYSHEV",
                                         "CHEBYSHEV-FOURIER",
                                         "FOURIER-MODIFIED"};

typedef std::map<int, std::vector<unsigned int>> CompositeOrdering;
typedef std::map<int, std::vector<unsigned int>> BndRegionOrdering;

// set restriction on domain range for post-processing.
struct DomainRange
{
    bool m_doXrange;
    NekDouble m_xmin;
    NekDouble m_xmax;
    bool m_doYrange;
    NekDouble m_ymin;
    NekDouble m_ymax;
    bool m_doZrange;
    NekDouble m_zmin;
    NekDouble m_zmax;

    bool m_checkShape;
    LibUtilities::ShapeType m_shapeType;
};

typedef std::shared_ptr<DomainRange> DomainRangeShPtr;
static DomainRangeShPtr NullDomainRangeShPtr;

struct Composite
{
    std::vector<std::shared_ptr<Geometry>> m_geomVec;
};

typedef std::shared_ptr<Composite> CompositeSharedPtr;
typedef std::map<int, CompositeSharedPtr> CompositeMap;

struct Expansion
{
    Expansion(GeometrySharedPtr geomShPtr,
              const LibUtilities::BasisKeyVector basiskeyvec)
        : m_geomShPtr(geomShPtr), m_basisKeyVector(basiskeyvec)
    {
    }

    GeometrySharedPtr m_geomShPtr;
    LibUtilities::BasisKeyVector m_basisKeyVector;
};

typedef std::shared_ptr<Expansion> ExpansionShPtr;
typedef std::map<int, ExpansionShPtr> ExpansionMap;

typedef std::shared_ptr<ExpansionMap> ExpansionMapShPtr;
typedef std::map<std::string, ExpansionMapShPtr> ExpansionMapShPtrMap;

typedef std::map<std::string, std::string> GeomInfoMap;
typedef std::shared_ptr<std::vector<std::pair<GeometrySharedPtr, int>>>
    GeometryLinkSharedPtr;

typedef std::map<std::string, std::string> MeshMetaDataMap;

class MeshGraph;
typedef std::shared_ptr<MeshGraph> MeshGraphSharedPtr;

/// Base class for a spectral/hp element mesh.
class MeshGraph
{
public:
    SPATIAL_DOMAINS_EXPORT MeshGraph();
    SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph();

    SPATIAL_DOMAINS_EXPORT static MeshGraphSharedPtr Read(
        const LibUtilities::SessionReaderSharedPtr pSession,
        DomainRangeShPtr                           rng       = NullDomainRangeShPtr,
        bool                                       fillGraph = true);

    SPATIAL_DOMAINS_EXPORT virtual void WriteGeometry(
        std::string &outfilename,
        bool defaultExp = false,
        const LibUtilities::FieldMetaDataMap &metadata
                                     = LibUtilities::NullFieldMetaDataMap) = 0;

    void Empty(int dim, int space)
    {
        m_meshDimension  = dim;
        m_spaceDimension = space;
    }

    /*transfers the minial data structure to full meshgraph*/
    SPATIAL_DOMAINS_EXPORT void FillGraph();

    SPATIAL_DOMAINS_EXPORT void FillBoundingBoxTree();

    SPATIAL_DOMAINS_EXPORT std::vector<int> GetElementsContainingPoint(
        PointGeomSharedPtr p);

    ////////////////////
    ////////////////////

    SPATIAL_DOMAINS_EXPORT void ReadExpansions();

    /* ---- Helper functions ---- */
    /// Dimension of the mesh (can be a 1D curve in 3D space).
    int GetMeshDimension()
    {
        return m_meshDimension;
    }

    /// Dimension of the space (can be a 1D curve in 3D space).
    int GetSpaceDimension()
    {
        return m_spaceDimension;
    }

    /* Range definitions for postprorcessing */
    SPATIAL_DOMAINS_EXPORT void SetDomainRange(
        NekDouble xmin, NekDouble xmax,
        NekDouble ymin = NekConstants::kNekUnsetDouble,
        NekDouble ymax = NekConstants::kNekUnsetDouble,
        NekDouble zmin = NekConstants::kNekUnsetDouble,
        NekDouble zmax = NekConstants::kNekUnsetDouble);

    /// Check if goemetry is in range definition if activated
    SPATIAL_DOMAINS_EXPORT bool CheckRange(Geometry2D &geom);

    /// Check if goemetry is in range definition if activated
    SPATIAL_DOMAINS_EXPORT bool CheckRange(Geometry3D &geom);

    /* ---- Composites and Domain ---- */
    CompositeSharedPtr GetComposite(int whichComposite)
    {
        ASSERTL0(m_meshComposites.find(whichComposite) !=
                     m_meshComposites.end(),
                 "Composite not found.");
        return m_meshComposites.find(whichComposite)->second;
    }

    SPATIAL_DOMAINS_EXPORT GeometrySharedPtr GetCompositeItem(
        int whichComposite, int whichItem);

    SPATIAL_DOMAINS_EXPORT void GetCompositeList(
        const std::string &compositeStr,
        CompositeMap &compositeVector) const;

    std::map<int, CompositeSharedPtr> &GetComposites()
    {
        return m_meshComposites;
    }

    std::map<int, std::string> &GetCompositesLabels()
    {
        return m_compositesLabels;
    }

    std::vector<std::map<int, CompositeSharedPtr>> &GetDomain()
    {
        return m_domain;
    }

    std::map<int, CompositeSharedPtr> &GetDomain(int domain)
    {
        ASSERTL1(domain < m_domain.size(),
                 "Request for domain which does not exist");
        return m_domain[domain];
    }

    SPATIAL_DOMAINS_EXPORT const ExpansionMap &GetExpansions(
        const std::string variable = "DefaultVar");

    SPATIAL_DOMAINS_EXPORT ExpansionShPtr GetExpansion(
        GeometrySharedPtr geom, const std::string variable = "DefaultVar");

    /// Sets expansions given field definitions
    SPATIAL_DOMAINS_EXPORT void SetExpansions(
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef);

    /// Sets expansions given field definition, quadrature points.
    SPATIAL_DOMAINS_EXPORT void SetExpansions(
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef,
        std::vector<std::vector<LibUtilities::PointsType>> &pointstype);

    /// Sets expansions to have equispaced points
    SPATIAL_DOMAINS_EXPORT void SetExpansionsToEvenlySpacedPoints(
        int npoints = 0);

    /// Reset expansion to have specified polynomial order \a nmodes
    SPATIAL_DOMAINS_EXPORT void SetExpansionsToPolyOrder(int nmodes);

    /// Reset expansion to have specified point order \a
    /// npts
    SPATIAL_DOMAINS_EXPORT void SetExpansionsToPointOrder(int npts);
    /// This function sets the expansion #exp in map with
    /// entry #variable

    inline void SetExpansions(const std::string variable,
                              ExpansionMapShPtr &exp);

    inline void SetSession(LibUtilities::SessionReaderSharedPtr pSession);

    /// Sets the basis key for all expansions of the given shape.
    SPATIAL_DOMAINS_EXPORT void SetBasisKey(LibUtilities::ShapeType shape,
                                            LibUtilities::BasisKeyVector &keys,
                                            std::string var = "DefaultVar");

    inline bool SameExpansions(const std::string var1, const std::string var2);

    inline bool CheckForGeomInfo(std::string parameter);

    inline const std::string GetGeomInfo(std::string parameter);

    SPATIAL_DOMAINS_EXPORT static LibUtilities::BasisKeyVector
    DefineBasisKeyFromExpansionType(GeometrySharedPtr in, ExpansionType type,
                                    const int order);

    SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKeyVector
    DefineBasisKeyFromExpansionTypeHomo(
        GeometrySharedPtr in, ExpansionType type_x, ExpansionType type_y,
        ExpansionType type_z, const int nummodes_x, const int nummodes_y,
        const int nummodes_z);

    /* ---- Manipulation of mesh ---- */
    int GetNvertices()
    {
        return m_vertSet.size();
    }

    PointGeomSharedPtr GetVertex(int id)
    {
        return m_vertSet[id];
    }

    SPATIAL_DOMAINS_EXPORT SegGeomSharedPtr GetSegGeom(int id)
    {
        return m_segGeoms[id];
    }

    SPATIAL_DOMAINS_EXPORT CurveMap &GetCurvedEdges()
    {
        return m_curvedEdges;
    }
    SPATIAL_DOMAINS_EXPORT CurveMap &GetCurvedFaces()
    {
        return m_curvedFaces;
    }

    SPATIAL_DOMAINS_EXPORT std::map<int, PointGeomSharedPtr> &GetAllPointGeoms()
    {
        return m_vertSet;
    }
    SPATIAL_DOMAINS_EXPORT std::map<int, SegGeomSharedPtr> &GetAllSegGeoms()
    {
        return m_segGeoms;
    }
    SPATIAL_DOMAINS_EXPORT TriGeomMap &GetAllTriGeoms()
    {
        return m_triGeoms;
    }
    SPATIAL_DOMAINS_EXPORT QuadGeomMap &GetAllQuadGeoms()
    {
        return m_quadGeoms;
    }
    SPATIAL_DOMAINS_EXPORT TetGeomMap &GetAllTetGeoms()
    {
        return m_tetGeoms;
    }
    SPATIAL_DOMAINS_EXPORT PyrGeomMap &GetAllPyrGeoms()
    {
        return m_pyrGeoms;
    }
    SPATIAL_DOMAINS_EXPORT PrismGeomMap &GetAllPrismGeoms()
    {
        return m_prismGeoms;
    }
    SPATIAL_DOMAINS_EXPORT HexGeomMap &GetAllHexGeoms()
    {
        return m_hexGeoms;
    }

    SPATIAL_DOMAINS_EXPORT int GetNumElements();

    Geometry2DSharedPtr GetGeometry2D(int gID)
    {
        auto it1 = m_triGeoms.find(gID);
        if (it1 != m_triGeoms.end())
            return it1->second;

        auto it2 = m_quadGeoms.find(gID);
        if (it2 != m_quadGeoms.end())
            return it2->second;

        return Geometry2DSharedPtr();
    };

    SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKey GetEdgeBasisKey(
        SegGeomSharedPtr edge, const std::string variable = "DefaultVar");

    SPATIAL_DOMAINS_EXPORT GeometryLinkSharedPtr GetElementsFromEdge(
        Geometry1DSharedPtr edge);

    SPATIAL_DOMAINS_EXPORT GeometryLinkSharedPtr GetElementsFromFace(
        Geometry2DSharedPtr face);

    SPATIAL_DOMAINS_EXPORT LibUtilities::BasisKey GetFaceBasisKey(
        Geometry2DSharedPtr face, const int facedir,
        const std::string variable = "DefaultVar");

    CompositeOrdering &GetCompositeOrdering()
    {
        return m_compOrder;
    }

    BndRegionOrdering &GetBndRegionOrdering()
    {
        return m_bndRegOrder;
    }

    /*an inital read which loads a very light weight data structure*/
    SPATIAL_DOMAINS_EXPORT virtual void ReadGeometry(
        DomainRangeShPtr rng,
        bool             fillGraph) = 0;
    SPATIAL_DOMAINS_EXPORT virtual void PartitionMesh(
        LibUtilities::SessionReaderSharedPtr session) = 0;

    SPATIAL_DOMAINS_EXPORT std::map<int, MeshEntity>
        CreateMeshEntities();
    SPATIAL_DOMAINS_EXPORT CompositeDescriptor CreateCompositeDescriptor();

protected:

    void PopulateFaceToElMap(Geometry3DSharedPtr element, int kNfaces);
    ExpansionMapShPtr SetUpExpansionMap();
    std::string GetCompositeString(CompositeSharedPtr comp);

    LibUtilities::SessionReaderSharedPtr m_session;
    PointGeomMap m_vertSet;

    CurveMap m_curvedEdges;
    CurveMap m_curvedFaces;

    SegGeomMap m_segGeoms;

    TriGeomMap m_triGeoms;
    QuadGeomMap m_quadGeoms;
    TetGeomMap m_tetGeoms;
    PyrGeomMap m_pyrGeoms;
    PrismGeomMap m_prismGeoms;
    HexGeomMap m_hexGeoms;

    int m_meshDimension;
    int m_spaceDimension;
    int m_partition;
    bool m_meshPartitioned;

    CompositeMap m_meshComposites;
    std::map<int, std::string> m_compositesLabels;
    std::vector<CompositeMap> m_domain;
    DomainRangeShPtr m_domainRange;

    ExpansionMapShPtrMap m_expansionMapShPtrMap;

    GeomInfoMap m_geomInfo;

    std::unordered_map<int, GeometryLinkSharedPtr> m_faceToElMap;

    TiXmlElement *m_xmlGeom;

    CompositeOrdering m_compOrder;
    BndRegionOrdering m_bndRegOrder;

    struct GeomRTree;
    std::unique_ptr<GeomRTree> m_boundingBoxTree;
};
typedef std::shared_ptr<MeshGraph> MeshGraphSharedPtr;
typedef LibUtilities::NekFactory<std::string, MeshGraph> MeshGraphFactory;

SPATIAL_DOMAINS_EXPORT MeshGraphFactory &GetMeshGraphFactory();

/**
 *
 */
void MeshGraph::SetExpansions(const std::string variable,
                              ExpansionMapShPtr &exp)
{
    if (m_expansionMapShPtrMap.count(variable) != 0)
    {
        ASSERTL0(false,
                 (std::string("Expansion field is already set for variable ") +
                  variable)
                     .c_str());
    }
    else
    {
        m_expansionMapShPtrMap[variable] = exp;
    }
}

/**
 *
 */
void MeshGraph::SetSession(LibUtilities::SessionReaderSharedPtr pSession)
{
    m_session = pSession;
}

/**
 *
 */
inline bool MeshGraph::SameExpansions(const std::string var1,
                                      const std::string var2)
{
    ExpansionMapShPtr expVec1 = m_expansionMapShPtrMap.find(var1)->second;
    ExpansionMapShPtr expVec2 = m_expansionMapShPtrMap.find(var2)->second;

    if (expVec1.get() == expVec2.get())
    {
        return true;
    }

    return false;
}

/**
 *
 */
inline bool MeshGraph::CheckForGeomInfo(std::string parameter)
{
    return m_geomInfo.find(parameter) != m_geomInfo.end();
}

/**
 *
 */
inline const std::string MeshGraph::GetGeomInfo(std::string parameter)
{
    ASSERTL1(m_geomInfo.find(parameter) != m_geomInfo.end(),
             "Parameter " + parameter + " does not exist.");
    return m_geomInfo[parameter];
}
}
}

#endif
