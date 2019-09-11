////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraph.cpp
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

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/BasicUtils/CompressData.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdTriExp.h>

#include <cstring>
#include <iomanip>
#include <sstream>

#include <SpatialDomains/MeshGraph.h>

// These are required for the Write(...) and Import(...) functions.
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
namespace bg = boost::geometry;

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

/**
 * Returns an instance of the MeshGraph factory, held as a singleton.
 */
MeshGraphFactory &GetMeshGraphFactory()
{
    static MeshGraphFactory instance;
    return instance;
}

struct MeshGraph::GeomRTree
{
    typedef bg::model::point<NekDouble, 3, bg::cs::cartesian> BgPoint;
    typedef bg::model::box<BgPoint> BgBox;
    typedef std::pair<BgBox, int> BgRtreeValue;

    bg::index::rtree< BgRtreeValue, bg::index::rstar<16, 4> > m_bgTree;

    void InsertGeom(GeometrySharedPtr const &geom)
    {
        std::array<NekDouble, 6> minMax = geom->GetBoundingBox();
        BgPoint ptMin(minMax[0], minMax[1], minMax[2]);
        BgPoint ptMax(minMax[3], minMax[4], minMax[5]);
        m_bgTree.insert(
            std::make_pair(BgBox(ptMin, ptMax), geom->GetGlobalID()));
    }
};

MeshGraph::MeshGraph()
{
    m_boundingBoxTree = std::unique_ptr<MeshGraph::GeomRTree>(
        new MeshGraph::GeomRTree());
}

/**
 *
 */
MeshGraph::~MeshGraph()
{
}

MeshGraphSharedPtr MeshGraph::Read(
    const LibUtilities::SessionReaderSharedPtr session,
    DomainRangeShPtr                           rng,
    bool                                       fillGraph)
{
    LibUtilities::CommSharedPtr comm = session->GetComm();
    ASSERTL0(comm.get(), "Communication not initialised.");

    // Populate SessionReader. This should be done only on the root process so
    // that we can partition appropriately without all processes having to read
    // in the input file.
    const bool isRoot = comm->TreatAsRankZero();
    std::string geomType;

    if (isRoot)
    {
        // Parse the XML document.
        session->InitSession();

        // Get geometry type, i.e. XML (compressed/uncompressed) or HDF5.
        geomType = session->GetGeometryType();

        // Convert to a vector of chars so that we can broadcast.
        std::vector<char> v(geomType.c_str(),
                            geomType.c_str() + geomType.length());

        size_t length = v.size();
        comm->Bcast(length, 0);
        comm->Bcast(v, 0);
    }
    else
    {
        size_t length;
        comm->Bcast(length, 0);

        std::vector<char> v(length);
        comm->Bcast(v, 0);

        geomType = std::string(v.begin(),v.end());
    }

    // Every process then creates a mesh. Partitioning logic takes place inside
    // the PartitionMesh function so that we can support different options for
    // XML and HDF5.
    MeshGraphSharedPtr mesh = GetMeshGraphFactory().CreateInstance(geomType);
    mesh->PartitionMesh(session);

    // Finally, read the geometry information.
    mesh->ReadGeometry(rng, fillGraph);

    return mesh;
}

void MeshGraph::FillGraph()
{
    ReadExpansions();

    switch (m_meshDimension)
    {
        case 3:
        {
            for (auto &x : m_pyrGeoms)
            {
                x.second->Setup();
            }
            for (auto &x : m_prismGeoms)
            {
                x.second->Setup();
            }
            for (auto &x : m_tetGeoms)
            {
                x.second->Setup();
            }
            for (auto &x : m_hexGeoms)
            {
                x.second->Setup();
            }
        }
        break;
        case 2:
        {
            for (auto &x : m_triGeoms)
            {
                x.second->Setup();
            }
            for (auto &x : m_quadGeoms)
            {
                x.second->Setup();
            }
        }
        break;
        case 1:
        {
            for (auto x : m_segGeoms)
            {
                x.second->Setup();
            }
        }
        break;
    }
}

void MeshGraph::FillBoundingBoxTree()
{

    m_boundingBoxTree->m_bgTree.clear();
    switch (m_meshDimension)
    {
        case 1:
            for (auto &x : m_segGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            break;
        case 2:
            for (auto &x : m_triGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            for (auto &x : m_quadGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            break;
        case 3:
            for (auto &x : m_tetGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            for (auto &x : m_prismGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            for (auto &x : m_pyrGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            for (auto &x : m_hexGeoms)
            {
                m_boundingBoxTree->InsertGeom(x.second);
            }
            break;
        default:
            ASSERTL0(false, "Unknown dim");
    }
}

std::vector<int> MeshGraph::GetElementsContainingPoint(
    PointGeomSharedPtr p)
{
    if (m_boundingBoxTree->m_bgTree.empty())
    {
        FillBoundingBoxTree();
    }

    NekDouble x = 0.0;
    NekDouble y = 0.0;
    NekDouble z = 0.0;
    std::vector<GeomRTree::BgRtreeValue> matches;

    p->GetCoords(x, y, z);

    GeomRTree::BgBox b(GeomRTree::BgPoint(x, y, z),
                       GeomRTree::BgPoint(x, y, z));

    m_boundingBoxTree->m_bgTree.query(bg::index::intersects(b),
                                      std::back_inserter(matches));

    std::vector<int> vals(matches.size());

    for (int i = 0; i < matches.size(); ++i)
    {
        vals[i] = matches[i].second;
    }

    return vals;
}

void MeshGraph::SetDomainRange(NekDouble xmin, NekDouble xmax, NekDouble ymin,
                               NekDouble ymax, NekDouble zmin, NekDouble zmax)
{
    m_domainRange->m_checkShape = false;

    if (m_domainRange == NullDomainRangeShPtr)
    {
        m_domainRange = MemoryManager<DomainRange>::AllocateSharedPtr();
        m_domainRange->m_doXrange = true;
    }

    m_domainRange->m_xmin = xmin;
    m_domainRange->m_xmax = xmax;

    if (ymin == NekConstants::kNekUnsetDouble)
    {
        m_domainRange->m_doYrange = false;
    }
    else
    {
        m_domainRange->m_doYrange = true;
        m_domainRange->m_ymin     = ymin;
        m_domainRange->m_ymax     = ymax;
    }

    if (zmin == NekConstants::kNekUnsetDouble)
    {
        m_domainRange->m_doZrange = false;
    }
    else
    {
        m_domainRange->m_doZrange = true;
        m_domainRange->m_zmin     = zmin;
        m_domainRange->m_zmax     = zmax;
    }
}

int MeshGraph::GetNumElements()
{
    switch (m_meshDimension)
    {
        case 1:
        {
            return m_segGeoms.size();
        }
        break;
        case 2:
        {
            return m_triGeoms.size() + m_quadGeoms.size();
        }
        break;
        case 3:
        {
            return m_tetGeoms.size() + m_pyrGeoms.size() + m_prismGeoms.size() +
                   m_hexGeoms.size();
        }
    }

    return 0;
}

bool MeshGraph::CheckRange(Geometry2D &geom)
{
    bool returnval = true;

    if (m_domainRange != NullDomainRangeShPtr)
    {
        int nverts  = geom.GetNumVerts();
        int coordim = geom.GetCoordim();

        // exclude elements outside x range if all vertices not in region
        if (m_domainRange->m_doXrange)
        {
            int ncnt_low = 0;
            int ncnt_up  = 0;
            for (int i = 0; i < nverts; ++i)
            {
                NekDouble xval = (*geom.GetVertex(i))[0];
                if (xval < m_domainRange->m_xmin)
                {
                    ncnt_low++;
                }

                if (xval > m_domainRange->m_xmax)
                {
                    ncnt_up++;
                }
            }

            // check for all verts to be less or greater than
            // range so that if element spans thin range then
            // it is still included
            if ((ncnt_up == nverts) || (ncnt_low == nverts))
            {
                returnval = false;
            }
        }

        // exclude elements outside y range if all vertices not in region
        if (m_domainRange->m_doYrange)
        {
            int ncnt_low = 0;
            int ncnt_up  = 0;
            for (int i = 0; i < nverts; ++i)
            {
                NekDouble yval = (*geom.GetVertex(i))[1];
                if (yval < m_domainRange->m_ymin)
                {
                    ncnt_low++;
                }

                if (yval > m_domainRange->m_ymax)
                {
                    ncnt_up++;
                }
            }

            // check for all verts to be less or greater than
            // range so that if element spans thin range then
            // it is still included
            if ((ncnt_up == nverts) || (ncnt_low == nverts))
            {
                returnval = false;
            }
        }

        if (coordim > 2)
        {
            // exclude elements outside z range if all vertices not in region
            if (m_domainRange->m_doZrange)
            {
                int ncnt_low = 0;
                int ncnt_up  = 0;

                for (int i = 0; i < nverts; ++i)
                {
                    NekDouble zval = (*geom.GetVertex(i))[2];

                    if (zval < m_domainRange->m_zmin)
                    {
                        ncnt_low++;
                    }

                    if (zval > m_domainRange->m_zmax)
                    {
                        ncnt_up++;
                    }
                }

                // check for all verts to be less or greater than
                // range so that if element spans thin range then
                // it is still included
                if ((ncnt_up == nverts) || (ncnt_low == nverts))
                {
                    returnval = false;
                }
            }
        }
    }
    return returnval;
}

/* Domain checker for 3D geometries */
bool MeshGraph::CheckRange(Geometry3D &geom)
{
    bool returnval = true;

    if (m_domainRange != NullDomainRangeShPtr)
    {
        int nverts = geom.GetNumVerts();

        if (m_domainRange->m_doXrange)
        {
            int ncnt_low = 0;
            int ncnt_up  = 0;

            for (int i = 0; i < nverts; ++i)
            {
                NekDouble xval = (*geom.GetVertex(i))[0];
                if (xval < m_domainRange->m_xmin)
                {
                    ncnt_low++;
                }

                if (xval > m_domainRange->m_xmax)
                {
                    ncnt_up++;
                }
            }

            // check for all verts to be less or greater than
            // range so that if element spans thin range then
            // it is still included
            if ((ncnt_up == nverts) || (ncnt_low == nverts))
            {
                returnval = false;
            }
        }

        if (m_domainRange->m_doYrange)
        {
            int ncnt_low = 0;
            int ncnt_up  = 0;
            for (int i = 0; i < nverts; ++i)
            {
                NekDouble yval = (*geom.GetVertex(i))[1];
                if (yval < m_domainRange->m_ymin)
                {
                    ncnt_low++;
                }

                if (yval > m_domainRange->m_ymax)
                {
                    ncnt_up++;
                }
            }

            // check for all verts to be less or greater than
            // range so that if element spans thin range then
            // it is still included
            if ((ncnt_up == nverts) || (ncnt_low == nverts))
            {
                returnval = false;
            }
        }

        if (m_domainRange->m_doZrange)
        {
            int ncnt_low = 0;
            int ncnt_up  = 0;
            for (int i = 0; i < nverts; ++i)
            {
                NekDouble zval = (*geom.GetVertex(i))[2];

                if (zval < m_domainRange->m_zmin)
                {
                    ncnt_low++;
                }

                if (zval > m_domainRange->m_zmax)
                {
                    ncnt_up++;
                }
            }

            // check for all verts to be less or greater than
            // range so that if element spans thin range then
            // it is still included
            if ((ncnt_up == nverts) || (ncnt_low == nverts))
            {
                returnval = false;
            }
        }

        if (m_domainRange->m_checkShape)
        {
            if (geom.GetShapeType() != m_domainRange->m_shapeType)
            {
                returnval = false;
            }
        }
    }

    return returnval;
}

/**
 *
 */
GeometrySharedPtr MeshGraph::GetCompositeItem(int whichComposite, int whichItem)
{
    GeometrySharedPtr returnval;
    bool error = false;

    if (whichComposite >= 0 && whichComposite < int(m_meshComposites.size()))
    {
        if (whichItem >= 0 &&
            whichItem < int(m_meshComposites[whichComposite]->m_geomVec.size()))
        {
            returnval = m_meshComposites[whichComposite]->m_geomVec[whichItem];
        }
        else
        {
            error = true;
        }
    }
    else
    {
        error = true;
    }

    if (error)
    {
        std::ostringstream errStream;
        errStream << "Unable to access composite item [" << whichComposite
                  << "][" << whichItem << "].";

        std::string testStr = errStream.str();

        NEKERROR(ErrorUtil::efatal, testStr.c_str());
    }

    return returnval;
}

/**
 *
 */
void MeshGraph::GetCompositeList(const std::string &compositeStr,
                                 CompositeMap &compositeVector) const
{
    // Parse the composites into a list.
    vector<unsigned int> seqVector;
    bool parseGood =
        ParseUtils::GenerateSeqVector(compositeStr.c_str(), seqVector);

    ASSERTL0(
        parseGood && !seqVector.empty(),
        (std::string("Unable to read composite index range: ") + compositeStr)
            .c_str());

    vector<unsigned int>
        addedVector; // Vector of those composites already added to
                     // compositeVector;
    for (auto iter = seqVector.begin(); iter != seqVector.end(); ++iter)
    {
        // Only add a new one if it does not already exist in vector.
        // Can't go back and delete with a vector, so prevent it from
        // being added in the first place.
        if (std::find(addedVector.begin(), addedVector.end(), *iter) ==
            addedVector.end())
        {

            // If the composite listed is not found and we are working
            // on a partitioned mesh, silently ignore it.
            if (m_meshComposites.find(*iter) == m_meshComposites.end() &&
                m_meshPartitioned)
            {
                continue;
            }

            addedVector.push_back(*iter);
            ASSERTL0(m_meshComposites.find(*iter) != m_meshComposites.end(),
                     "Composite not found.");
            CompositeSharedPtr composite = m_meshComposites.find(*iter)->second;

            if (composite)
            {
                compositeVector[*iter] = composite;
            }
            else
            {
                char str[64];
                ::sprintf(str, "%d", *iter);
                NEKERROR(ErrorUtil::ewarning,
                         (std::string("Undefined composite: ") + str).c_str());
            }
        }
    }
}

/**
 *
 */
const ExpansionMap &MeshGraph::GetExpansions(const std::string variable)
{
    ExpansionMapShPtr returnval;

    if (m_expansionMapShPtrMap.count(variable))
    {
        returnval = m_expansionMapShPtrMap.find(variable)->second;
    }
    else
    {
        if (m_expansionMapShPtrMap.count("DefaultVar") == 0)
        {
            NEKERROR(
                ErrorUtil::efatal,
                (std::string(
                     "Unable to find expansion vector definition for field: ") +
                 variable)
                    .c_str());
        }
        returnval = m_expansionMapShPtrMap.find("DefaultVar")->second;
        m_expansionMapShPtrMap[variable] = returnval;

        NEKERROR(
            ErrorUtil::ewarning,
            (std::string(
                 "Using Default variable expansion definition for field: ") +
             variable)
                .c_str());
    }

    return *returnval;
}

/**
 *
 */
ExpansionShPtr MeshGraph::GetExpansion(GeometrySharedPtr geom,
                                       const std::string variable)
{
    ExpansionMapShPtr expansionMap =
        m_expansionMapShPtrMap.find(variable)->second;

    auto iter = expansionMap->find(geom->GetGlobalID());
    ASSERTL1(iter != expansionMap->end(),
             "Could not find expansion " +
                 boost::lexical_cast<string>(geom->GetGlobalID()) +
                 " in expansion for variable " + variable);
    return iter->second;
}

/**
 *
 */
void MeshGraph::SetExpansions(
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
{
    int i, j, k, cnt, id;
    GeometrySharedPtr geom;

    ExpansionMapShPtr expansionMap;

    // Loop over fields and determine unique fields string and
    // declare whole expansion list
    for (i = 0; i < fielddef.size(); ++i)
    {
        for (j = 0; j < fielddef[i]->m_fields.size(); ++j)
        {
            std::string field = fielddef[i]->m_fields[j];
            if (m_expansionMapShPtrMap.count(field) == 0)
            {
                expansionMap                  = SetUpExpansionMap();
                m_expansionMapShPtrMap[field] = expansionMap;

                // check to see if DefaultVar also not set and
                // if so assign it to this expansion
                if (m_expansionMapShPtrMap.count("DefaultVar") == 0)
                {
                    m_expansionMapShPtrMap["DefaultVar"] = expansionMap;
                }
            }
        }
    }

    // loop over all elements find the geometry shared ptr and
    // set up basiskey vector
    for (i = 0; i < fielddef.size(); ++i)
    {
        cnt                                        = 0;
        std::vector<std::string> fields            = fielddef[i]->m_fields;
        std::vector<unsigned int> nmodes           = fielddef[i]->m_numModes;
        std::vector<LibUtilities::BasisType> basis = fielddef[i]->m_basis;
        bool pointDef                              = fielddef[i]->m_pointsDef;
        bool numPointDef = fielddef[i]->m_numPointsDef;

        // Check points and numpoints
        std::vector<unsigned int> npoints            = fielddef[i]->m_numPoints;
        std::vector<LibUtilities::PointsType> points = fielddef[i]->m_points;

        bool UniOrder = fielddef[i]->m_uniOrder;

        for (j = 0; j < fielddef[i]->m_elementIDs.size(); ++j)
        {

            LibUtilities::BasisKeyVector bkeyvec;
            id = fielddef[i]->m_elementIDs[j];

            switch (fielddef[i]->m_shapeType)
            {
                case LibUtilities::eSegment:
                {
                    if (m_segGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                    {
                        // skip element likely from parallel read
                        if (!UniOrder)
                        {
                            cnt++;
                            cnt += fielddef[i]->m_numHomogeneousDir;
                        }
                        continue;
                    }
                    geom = m_segGeoms[fielddef[i]->m_elementIDs[j]];

                    LibUtilities::PointsKey pkey(
                        nmodes[cnt] + 1, LibUtilities::eGaussLobattoLegendre);

                    if (numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey1(npoints[cnt],
                                                            points[0]);
                        pkey = pkey1;
                    }
                    else if (!numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey1(nmodes[cnt] + 1,
                                                            points[0]);
                        pkey = pkey1;
                    }
                    else if (numPointDef && !pointDef)
                    {
                        const LibUtilities::PointsKey pkey1(
                            npoints[cnt], LibUtilities::eGaussLobattoLegendre);
                        pkey = pkey1;
                    }

                    LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);

                    if (!UniOrder)
                    {
                        cnt++;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                    bkeyvec.push_back(bkey);
                }
                break;
                case LibUtilities::eTriangle:
                {
                    if (m_triGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                    {
                        // skip element likely from parallel read
                        if (!UniOrder)
                        {
                            cnt += 2;
                            cnt += fielddef[i]->m_numHomogeneousDir;
                        }
                        continue;
                    }
                    geom = m_triGeoms[fielddef[i]->m_elementIDs[j]];

                    LibUtilities::PointsKey pkey(
                        nmodes[cnt] + 1, LibUtilities::eGaussLobattoLegendre);
                    if (numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(npoints[cnt],
                                                            points[0]);
                        pkey = pkey2;
                    }
                    else if (!numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(nmodes[cnt] + 1,
                                                            points[0]);
                        pkey = pkey2;
                    }
                    else if (numPointDef && !pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(
                            npoints[cnt], LibUtilities::eGaussLobattoLegendre);
                        pkey = pkey2;
                    }
                    LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);

                    bkeyvec.push_back(bkey);

                    LibUtilities::PointsKey pkey1(
                        nmodes[cnt + 1], LibUtilities::eGaussRadauMAlpha1Beta0);
                    if (numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(npoints[cnt + 1],
                                                            points[1]);
                        pkey1 = pkey2;
                    }
                    else if (!numPointDef && pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(nmodes[cnt + 1],
                                                            points[1]);
                        pkey1 = pkey2;
                    }
                    else if (numPointDef && !pointDef)
                    {
                        const LibUtilities::PointsKey pkey2(
                            npoints[cnt + 1],
                            LibUtilities::eGaussRadauMAlpha1Beta0);
                        pkey1 = pkey2;
                    }
                    LibUtilities::BasisKey bkey1(basis[1], nmodes[cnt + 1],
                                                 pkey1);
                    bkeyvec.push_back(bkey1);

                    if (!UniOrder)
                    {
                        cnt += 2;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    if (m_quadGeoms.count(fielddef[i]->m_elementIDs[j]) == 0)
                    {
                        // skip element likely from parallel read
                        if (!UniOrder)
                        {
                            cnt += 2;
                            cnt += fielddef[i]->m_numHomogeneousDir;
                        }
                        continue;
                    }

                    geom = m_quadGeoms[fielddef[i]->m_elementIDs[j]];

                    for (int b = 0; b < 2; ++b)
                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + b] + 1,
                            LibUtilities::eGaussLobattoLegendre);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b], points[b]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + b] + 1, points[b]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 2;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                }
                break;

                case LibUtilities::eTetrahedron:
                {
                    k = fielddef[i]->m_elementIDs[j];

                    // allow for possibility that fielddef is
                    // larger than m_graph which can happen in
                    // parallel runs
                    if (m_tetGeoms.count(k) == 0)
                    {
                        if (!UniOrder)
                        {
                            cnt += 3;
                        }
                        continue;
                    }
                    geom = m_tetGeoms[k];

                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt] + 1,
                            LibUtilities::eGaussLobattoLegendre);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(npoints[cnt],
                                                                points[0]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(nmodes[cnt] + 1,
                                                                points[0]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[0], nmodes[cnt],
                                                    pkey);

                        bkeyvec.push_back(bkey);
                    }
                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + 1],
                            LibUtilities::eGaussRadauMAlpha1Beta0);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 1], points[1]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + 1] + 1, points[1]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 1],
                                LibUtilities::eGaussRadauMAlpha1Beta0);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[1], nmodes[cnt + 1],
                                                    pkey);

                        bkeyvec.push_back(bkey);
                    }

                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + 2],
                            LibUtilities::eGaussRadauMAlpha2Beta0);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2], points[2]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + 2] + 1, points[2]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2],
                                LibUtilities::eGaussRadauMAlpha1Beta0);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[2], nmodes[cnt + 2],
                                                    pkey);

                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::ePrism:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    if (m_prismGeoms.count(k) == 0)
                    {
                        if (!UniOrder)
                        {
                            cnt += 3;
                        }
                        continue;
                    }
                    geom = m_prismGeoms[k];

                    for (int b = 0; b < 2; ++b)
                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + b] + 1,
                            LibUtilities::eGaussLobattoLegendre);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b], points[b]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + b] + 1, points[b]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + 2],
                            LibUtilities::eGaussRadauMAlpha1Beta0);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2], points[2]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + 2] + 1, points[2]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[2], nmodes[cnt + 2],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::ePyramid:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_pyrGeoms.find(k) != m_pyrGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_pyrGeoms[k];

                    for (int b = 0; b < 2; ++b)
                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + b] + 1,
                            LibUtilities::eGaussLobattoLegendre);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b], points[b]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + b] + 1, points[b]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + 2],
                            LibUtilities::eGaussRadauMAlpha2Beta0);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2], points[2]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + 2] + 1, points[2]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + 2],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[2], nmodes[cnt + 2],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    if (m_hexGeoms.count(k) == 0)
                    {
                        if (!UniOrder)
                        {
                            cnt += 3;
                        }
                        continue;
                    }

                    geom = m_hexGeoms[k];

                    for (int b = 0; b < 3; ++b)
                    {
                        LibUtilities::PointsKey pkey(
                            nmodes[cnt + b],
                            LibUtilities::eGaussLobattoLegendre);

                        if (numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b], points[b]);
                            pkey = pkey2;
                        }
                        else if (!numPointDef && pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                nmodes[cnt + b] + 1, points[b]);
                            pkey = pkey2;
                        }
                        else if (numPointDef && !pointDef)
                        {
                            const LibUtilities::PointsKey pkey2(
                                npoints[cnt + b],
                                LibUtilities::eGaussLobattoLegendre);
                            pkey = pkey2;
                        }

                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                default:
                    ASSERTL0(
                        false,
                        "Need to set up for pyramid and prism 3D Expansions");
                    break;
            }

            for (k = 0; k < fields.size(); ++k)
            {
                expansionMap = m_expansionMapShPtrMap.find(fields[k])->second;
                if ((*expansionMap).find(id) != (*expansionMap).end())
                {
                    (*expansionMap)[id]->m_geomShPtr      = geom;
                    (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                }
            }
        }
    }
}

/**
 *
 */
void MeshGraph::SetExpansions(
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef,
    std::vector<std::vector<LibUtilities::PointsType>> &pointstype)
{
    int i, j, k, cnt, id;
    GeometrySharedPtr geom;

    ExpansionMapShPtr expansionMap;

    // Loop over fields and determine unique fields string and
    // declare whole expansion list
    for (i = 0; i < fielddef.size(); ++i)
    {
        for (j = 0; j < fielddef[i]->m_fields.size(); ++j)
        {
            std::string field = fielddef[i]->m_fields[j];
            if (m_expansionMapShPtrMap.count(field) == 0)
            {
                expansionMap                  = SetUpExpansionMap();
                m_expansionMapShPtrMap[field] = expansionMap;

                // check to see if DefaultVar also not set and
                // if so assign it to this expansion
                if (m_expansionMapShPtrMap.count("DefaultVar") == 0)
                {
                    m_expansionMapShPtrMap["DefaultVar"] = expansionMap;
                }
            }
        }
    }

    // loop over all elements find the geometry shared ptr and
    // set up basiskey vector
    for (i = 0; i < fielddef.size(); ++i)
    {
        cnt                                        = 0;
        std::vector<std::string> fields            = fielddef[i]->m_fields;
        std::vector<unsigned int> nmodes           = fielddef[i]->m_numModes;
        std::vector<LibUtilities::BasisType> basis = fielddef[i]->m_basis;
        bool UniOrder                              = fielddef[i]->m_uniOrder;

        for (j = 0; j < fielddef[i]->m_elementIDs.size(); ++j)
        {
            LibUtilities::BasisKeyVector bkeyvec;
            id = fielddef[i]->m_elementIDs[j];

            switch (fielddef[i]->m_shapeType)
            {
                case LibUtilities::eSegment:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_segGeoms.find(k) != m_segGeoms.end(),
                             "Failed to find geometry with same global id.");
                    geom = m_segGeoms[k];

                    const LibUtilities::PointsKey pkey(nmodes[cnt],
                                                       pointstype[i][0]);
                    LibUtilities::BasisKey bkey(basis[0], nmodes[cnt], pkey);
                    if (!UniOrder)
                    {
                        cnt++;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                    bkeyvec.push_back(bkey);
                }
                break;
                case LibUtilities::eTriangle:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_triGeoms.find(k) != m_triGeoms.end(),
                             "Failed to find geometry with same global id.");
                    geom = m_triGeoms[k];
                    for (int b = 0; b < 2; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 2;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_quadGeoms.find(k) != m_quadGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_quadGeoms[k];

                    for (int b = 0; b < 2; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 2;
                        cnt += fielddef[i]->m_numHomogeneousDir;
                    }
                }
                break;
                case LibUtilities::eTetrahedron:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_tetGeoms.find(k) != m_tetGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_tetGeoms[k];

                    for (int b = 0; b < 3; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::ePyramid:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_pyrGeoms.find(k) != m_pyrGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_pyrGeoms[k];

                    for (int b = 0; b < 3; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::ePrism:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_prismGeoms.find(k) != m_prismGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_prismGeoms[k];

                    for (int b = 0; b < 3; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    k = fielddef[i]->m_elementIDs[j];
                    ASSERTL0(m_hexGeoms.find(k) != m_hexGeoms.end(),
                             "Failed to find geometry with same global id");
                    geom = m_hexGeoms[k];

                    for (int b = 0; b < 3; ++b)
                    {
                        const LibUtilities::PointsKey pkey(nmodes[cnt + b],
                                                           pointstype[i][b]);
                        LibUtilities::BasisKey bkey(basis[b], nmodes[cnt + b],
                                                    pkey);
                        bkeyvec.push_back(bkey);
                    }

                    if (!UniOrder)
                    {
                        cnt += 3;
                    }
                }
                break;
                default:
                    ASSERTL0(
                        false,
                        "Need to set up for pyramid and prism 3D Expansions");
                    break;
            }

            for (k = 0; k < fields.size(); ++k)
            {
                expansionMap = m_expansionMapShPtrMap.find(fields[k])->second;
                if ((*expansionMap).find(id) != (*expansionMap).end())
                {
                    (*expansionMap)[id]->m_geomShPtr      = geom;
                    (*expansionMap)[id]->m_basisKeyVector = bkeyvec;
                }
            }
        }
    }
}

/**
 * \brief Reset all points keys to have equispaced points with
 * optional arguemt of \a npoints which redefines how many
 * points are to be used.
 */
void MeshGraph::SetExpansionsToEvenlySpacedPoints(int npoints)
{
    // iterate over all defined expansions
    for (auto it = m_expansionMapShPtrMap.begin();
         it != m_expansionMapShPtrMap.end(); ++it)
    {
        for (auto expIt = it->second->begin(); expIt != it->second->end();
             ++expIt)
        {
            for (int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey bkeyold =
                    expIt->second->m_basisKeyVector[i];

                int npts;

                if (npoints) // use input
                {
                    npts = npoints;
                }
                else
                {
                    npts = bkeyold.GetNumModes();
                }
                npts = max(npts, 2);

                const LibUtilities::PointsKey pkey(
                    npts, LibUtilities::ePolyEvenlySpaced);
                LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),
                                               bkeyold.GetNumModes(), pkey);
                expIt->second->m_basisKeyVector[i] = bkeynew;
            }
        }
    }
}

/**
 * \brief Reset all points keys to have expansion order of \a
 *  nmodes.  we keep the point distribution the same and make
 *  the number of points the same difference from the number
 *  of modes as the original expansion definition
 */
void MeshGraph::SetExpansionsToPolyOrder(int nmodes)
{
    // iterate over all defined expansions
    for (auto it = m_expansionMapShPtrMap.begin();
         it != m_expansionMapShPtrMap.end(); ++it)
    {
        for (auto expIt = it->second->begin(); expIt != it->second->end();
             ++expIt)
        {
            for (int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey bkeyold =
                    expIt->second->m_basisKeyVector[i];

                int npts =
                    nmodes + (bkeyold.GetNumPoints() - bkeyold.GetNumModes());

                const LibUtilities::PointsKey pkey(npts,
                                                   bkeyold.GetPointsType());
                LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(), nmodes,
                                               pkey);
                expIt->second->m_basisKeyVector[i] = bkeynew;
            }
        }
    }
}

/**
 * \brief Reset all points keys to have expansion order of \a
 *  nmodes.  we keep the point distribution the same and make
 *  the number of points the same difference from the number
 *  of modes as the original expansion definition
 */
void MeshGraph::SetExpansionsToPointOrder(int npts)
{
    // iterate over all defined expansions
    for (auto it = m_expansionMapShPtrMap.begin();
         it != m_expansionMapShPtrMap.end(); ++it)
    {
        for (auto expIt = it->second->begin(); expIt != it->second->end();
             ++expIt)
        {
            for (int i = 0; i < expIt->second->m_basisKeyVector.size(); ++i)
            {
                LibUtilities::BasisKey bkeyold =
                    expIt->second->m_basisKeyVector[i];

                const LibUtilities::PointsKey pkey(npts,
                                                   bkeyold.GetPointsType());

                LibUtilities::BasisKey bkeynew(bkeyold.GetBasisType(),
                                               bkeyold.GetNumModes(), pkey);
                expIt->second->m_basisKeyVector[i] = bkeynew;
            }
        }
    }
}

/**
 * For each element of shape given by \a shape in field \a
 * var, replace the current BasisKeyVector describing the
 * expansion in each dimension, with the one provided by \a
 * keys.
 *
 * @TODO: Allow selection of elements through a CompositeVector,
 * as well as by type.
 *
 * @param   shape     The shape of elements to be changed.
 * @param   keys      The new basis vector to apply to those elements.
 */
void MeshGraph::SetBasisKey(LibUtilities::ShapeType shape,
                            LibUtilities::BasisKeyVector &keys, std::string var)
{
    ExpansionMapShPtr expansionMap = m_expansionMapShPtrMap.find(var)->second;

    for (auto elemIter = expansionMap->begin(); elemIter != expansionMap->end();
         ++elemIter)
    {
        if ((elemIter->second)->m_geomShPtr->GetShapeType() == shape)
        {
            (elemIter->second)->m_basisKeyVector = keys;
        }
    }
}

/**
 *
 */
LibUtilities::BasisKeyVector MeshGraph::DefineBasisKeyFromExpansionType(
    GeometrySharedPtr in, ExpansionType type, const int nummodes)
{
    LibUtilities::BasisKeyVector returnval;

    LibUtilities::ShapeType shape = in->GetShapeType();

    int quadoffset = 1;
    switch (type)
    {
        case eModified:
        case eModifiedGLLRadau10:
            quadoffset = 1;
            break;
        case eModifiedQuadPlus1:
            quadoffset = 2;
            break;
        case eModifiedQuadPlus2:
            quadoffset = 3;
            break;
        default:
            break;
    }

    switch (type)
    {
        case eModified:
        case eModifiedQuadPlus1:
        case eModifiedQuadPlus2:
        case eModifiedGLLRadau10:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eTriangle:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes + quadoffset - 1,
                        LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModified_B,
                                                 nummodes, pkey1);

                    returnval.push_back(bkey1);
                }
                break;
                case LibUtilities::eTetrahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes + quadoffset - 1,
                        LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModified_B,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);

                    if (type == eModifiedGLLRadau10)
                    {
                        const LibUtilities::PointsKey pkey2(
                            nummodes + quadoffset - 1,
                            LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey bkey2(LibUtilities::eModified_C,
                                                     nummodes, pkey2);
                        returnval.push_back(bkey2);
                    }
                    else
                    {
                        const LibUtilities::PointsKey pkey2(
                            nummodes + quadoffset - 1,
                            LibUtilities::eGaussRadauMAlpha2Beta0);
                        LibUtilities::BasisKey bkey2(LibUtilities::eModified_C,
                                                     nummodes, pkey2);
                        returnval.push_back(bkey2);
                    }
                }
                break;
                case LibUtilities::ePyramid:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes + quadoffset,
                        LibUtilities::eGaussRadauMAlpha2Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModifiedPyr_C,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                case LibUtilities::ePrism:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + quadoffset,
                        LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eModified_A,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes + quadoffset - 1,
                        LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModified_B,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch for this shape");
                }
                break;
            }
        }
        break;

        case eGLL_Lagrange:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eTriangle: // define with corrects points key
                    // and change to Ortho on construction
                    {
                        const LibUtilities::PointsKey pkey(
                            nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                        LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                    nummodes, pkey);
                        returnval.push_back(bkey);

                        const LibUtilities::PointsKey pkey1(
                            nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                        LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,
                                                     nummodes, pkey1);
                        returnval.push_back(bkey1);
                    }
                    break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eGauss_Lagrange:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGauss_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eOrthogonal:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eTriangle:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,
                                                nummodes, pkey);

                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,
                                                 nummodes, pkey1);

                    returnval.push_back(bkey1);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eTetrahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eOrtho_A,
                                                nummodes, pkey);

                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes, LibUtilities::eGaussRadauMAlpha1Beta0);
                    LibUtilities::BasisKey bkey1(LibUtilities::eOrtho_B,
                                                 nummodes, pkey1);

                    returnval.push_back(bkey1);

                    const LibUtilities::PointsKey pkey2(
                        nummodes, LibUtilities::eGaussRadauMAlpha2Beta0);
                    LibUtilities::BasisKey bkey2(LibUtilities::eOrtho_C,
                                                 nummodes, pkey2);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eGLL_Lagrange_SEM:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey(LibUtilities::eGLL_Lagrange,
                                                nummodes, pkey);

                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourier:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourierSingleMode:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierSingleMode, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierSingleMode, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierSingleMode, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourierHalfModeRe:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeRe, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourierHalfModeIm:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey(
                        LibUtilities::eFourierHalfModeIm, nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eChebyshev:
        {
            switch (shape)
            {
                case LibUtilities::eSegment:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                case LibUtilities::eHexahedron:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey(LibUtilities::eChebyshev,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourierChebyshev:
        {
            switch (shape)
            {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eChebyshevFourier:
        {
            switch (shape)
            {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);

                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        case eFourierModified:
        {
            switch (shape)
            {
                case LibUtilities::eQuadrilateral:
                {
                    const LibUtilities::PointsKey pkey(
                        nummodes, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey(LibUtilities::eFourier,
                                                nummodes, pkey);
                    returnval.push_back(bkey);

                    const LibUtilities::PointsKey pkey1(
                        nummodes + 1, LibUtilities::eGaussLobattoLegendre);
                    LibUtilities::BasisKey bkey1(LibUtilities::eModified_A,
                                                 nummodes, pkey1);
                    returnval.push_back(bkey1);
                }
                break;
                default:
                {
                    ASSERTL0(false,
                             "Expansion not defined in switch  for this shape");
                }
                break;
            }
        }
        break;

        default:
        {
            ASSERTL0(false, "Expansion type not defined");
        }
        break;
    }

    return returnval;
}

/**
 *
 */
LibUtilities::BasisKeyVector MeshGraph::DefineBasisKeyFromExpansionTypeHomo(
    GeometrySharedPtr in, ExpansionType type_x, ExpansionType type_y,
    ExpansionType type_z, const int nummodes_x, const int nummodes_y,
    const int nummodes_z)
{
    LibUtilities::BasisKeyVector returnval;

    LibUtilities::ShapeType shape = in->GetShapeType();

    switch (shape)
    {
        case LibUtilities::eSegment:
        {
            ASSERTL0(false, "Homogeneous expansion not defined for this shape");
        }
        break;

        case LibUtilities::eQuadrilateral:
        {
            ASSERTL0(false, "Homogeneous expansion not defined for this shape");
        }
        break;

        case LibUtilities::eHexahedron:
        {
            switch (type_x)
            {
                case eFourier:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes_x, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey1(LibUtilities::eFourier,
                                                 nummodes_x, pkey1);
                    returnval.push_back(bkey1);
                }
                break;

                case eFourierSingleMode:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes_x, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey1(
                        LibUtilities::eFourierSingleMode, nummodes_x, pkey1);
                    returnval.push_back(bkey1);
                }
                break;

                case eFourierHalfModeRe:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes_x, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey1(
                        LibUtilities::eFourierHalfModeRe, nummodes_x, pkey1);
                    returnval.push_back(bkey1);
                }
                break;

                case eFourierHalfModeIm:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes_x, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey1(
                        LibUtilities::eFourierHalfModeIm, nummodes_x, pkey1);
                    returnval.push_back(bkey1);
                }
                break;

                case eChebyshev:
                {
                    const LibUtilities::PointsKey pkey1(
                        nummodes_x, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey1(LibUtilities::eChebyshev,
                                                 nummodes_x, pkey1);
                    returnval.push_back(bkey1);
                }
                break;

                default:
                {
                    ASSERTL0(false, "Homogeneous expansion can be of Fourier "
                                    "or Chebyshev type only");
                }
                break;
            }

            switch (type_y)
            {
                case eFourier:
                {
                    const LibUtilities::PointsKey pkey2(
                        nummodes_y, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey2(LibUtilities::eFourier,
                                                 nummodes_y, pkey2);
                    returnval.push_back(bkey2);
                }
                break;

                case eFourierSingleMode:
                {
                    const LibUtilities::PointsKey pkey2(
                        nummodes_y, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey2(
                        LibUtilities::eFourierSingleMode, nummodes_y, pkey2);
                    returnval.push_back(bkey2);
                }
                break;

                case eFourierHalfModeRe:
                {
                    const LibUtilities::PointsKey pkey2(
                        nummodes_y, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey2(
                        LibUtilities::eFourierHalfModeRe, nummodes_y, pkey2);
                    returnval.push_back(bkey2);
                }
                break;

                case eFourierHalfModeIm:
                {
                    const LibUtilities::PointsKey pkey2(
                        nummodes_y, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey2(
                        LibUtilities::eFourierHalfModeIm, nummodes_y, pkey2);
                    returnval.push_back(bkey2);
                }
                break;

                case eChebyshev:
                {
                    const LibUtilities::PointsKey pkey2(
                        nummodes_y, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey2(LibUtilities::eChebyshev,
                                                 nummodes_y, pkey2);
                    returnval.push_back(bkey2);
                }
                break;

                default:
                {
                    ASSERTL0(false, "Homogeneous expansion can be of Fourier "
                                    "or Chebyshev type only");
                }
                break;
            }

            switch (type_z)
            {
                case eFourier:
                {
                    const LibUtilities::PointsKey pkey3(
                        nummodes_z, LibUtilities::eFourierEvenlySpaced);
                    LibUtilities::BasisKey bkey3(LibUtilities::eFourier,
                                                 nummodes_z, pkey3);
                    returnval.push_back(bkey3);
                }
                break;

                case eFourierSingleMode:
                {
                    const LibUtilities::PointsKey pkey3(
                        nummodes_z, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey3(
                        LibUtilities::eFourierSingleMode, nummodes_z, pkey3);
                    returnval.push_back(bkey3);
                }
                break;

                case eFourierHalfModeRe:
                {
                    const LibUtilities::PointsKey pkey3(
                        nummodes_z, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey3(
                        LibUtilities::eFourierHalfModeRe, nummodes_z, pkey3);
                    returnval.push_back(bkey3);
                }
                break;

                case eFourierHalfModeIm:
                {
                    const LibUtilities::PointsKey pkey3(
                        nummodes_z, LibUtilities::eFourierSingleModeSpaced);
                    LibUtilities::BasisKey bkey3(
                        LibUtilities::eFourierHalfModeIm, nummodes_z, pkey3);
                    returnval.push_back(bkey3);
                }
                break;

                case eChebyshev:
                {
                    const LibUtilities::PointsKey pkey3(
                        nummodes_z, LibUtilities::eGaussGaussChebyshev);
                    LibUtilities::BasisKey bkey3(LibUtilities::eChebyshev,
                                                 nummodes_z, pkey3);
                    returnval.push_back(bkey3);
                }
                break;

                default:
                {
                    ASSERTL0(false, "Homogeneous expansion can be of Fourier "
                                    "or Chebyshev type only");
                }
                break;
            }
        }
        break;

        case LibUtilities::eTriangle:
        {
            ASSERTL0(false, "Homogeneous expansion not defined for this shape");
        }
        break;

        case LibUtilities::eTetrahedron:
        {
            ASSERTL0(false, "Homogeneous expansion not defined for this shape");
        }
        break;

        default:
            ASSERTL0(false, "Expansion not defined in switch  for this shape");
            break;
    }

    return returnval;
}

/**
 * Generate a single vector of Expansion structs mapping global element
 * ID to a corresponding Geometry shared pointer and basis key.
 *
 * Expansion map ensures elements which appear in multiple composites
 * within the domain are only listed once.
 */
ExpansionMapShPtr MeshGraph::SetUpExpansionMap(void)
{
    ExpansionMapShPtr returnval;
    returnval = MemoryManager<ExpansionMap>::AllocateSharedPtr();

    for (int d = 0; d < m_domain.size(); ++d)
    {
        for (auto compIter = m_domain[d].begin(); compIter != m_domain[d].end();
             ++compIter)
        {
            for (auto x = compIter->second->m_geomVec.begin();
                 x != compIter->second->m_geomVec.end(); ++x)
            {
                LibUtilities::BasisKeyVector def;
                ExpansionShPtr expansionElementShPtr =
                    MemoryManager<Expansion>::AllocateSharedPtr(*x, def);
                int id           = (*x)->GetGlobalID();
                (*returnval)[id] = expansionElementShPtr;
            }
        }
    }

    return returnval;
}

/**
 * @brief Returns a string representation of a composite.
 */
std::string MeshGraph::GetCompositeString(CompositeSharedPtr comp)
{
    if (comp->m_geomVec.size() == 0)
    {
        return "";
    }

    // Create a map that gets around the issue of mapping faces -> F and edges
    // -> E inside the tag.
    map<LibUtilities::ShapeType, pair<string, string>> compMap;
    compMap[LibUtilities::ePoint]         = make_pair("V", "V");
    compMap[LibUtilities::eSegment]       = make_pair("S", "E");
    compMap[LibUtilities::eQuadrilateral] = make_pair("Q", "F");
    compMap[LibUtilities::eTriangle]      = make_pair("T", "F");
    compMap[LibUtilities::eTetrahedron]   = make_pair("A", "A");
    compMap[LibUtilities::ePyramid]       = make_pair("P", "P");
    compMap[LibUtilities::ePrism]         = make_pair("R", "R");
    compMap[LibUtilities::eHexahedron]    = make_pair("H", "H");

    stringstream s;

    GeometrySharedPtr firstGeom = comp->m_geomVec[0];
    int shapeDim                = firstGeom->GetShapeDim();
    string tag                  = (shapeDim < m_meshDimension)
        ? compMap[firstGeom->GetShapeType()].second
        : compMap[firstGeom->GetShapeType()].first;

    std::vector<unsigned int> idxList;
    std::transform(
        comp->m_geomVec.begin(), comp->m_geomVec.end(),
        std::back_inserter(idxList),
        [] (GeometrySharedPtr geom) { return geom->GetGlobalID(); });

    s << " " << tag << "[" << ParseUtils::GenerateSeqString(idxList) << "] ";
    return s.str();
}

void MeshGraph::ReadExpansions()
{
    // Find the Expansions tag
    TiXmlElement *expansionTypes = m_session->GetElement("NEKTAR/EXPANSIONS");
    ASSERTL0(expansionTypes, "Unable to find EXPANSIONS tag in file.");

    if (expansionTypes)
    {
        // Find the Expansion type
        TiXmlElement *expansion = expansionTypes->FirstChildElement();
        std::string expType     = expansion->Value();

        if (expType == "E")
        {
            int i;
            ExpansionMapShPtr expansionMap;

            /// Expansiontypes will contain composite,
            /// nummodes, and expansiontype (eModified, or
            /// eOrthogonal) Or a full list of data of
            /// basistype, nummodes, pointstype, numpoints;

            /// Expansiontypes may also contain a list of
            /// fields that this expansion relates to. If this
            /// does not exist the variable is only set to
            /// "DefaultVar".

            while (expansion)
            {

                const char *fStr = expansion->Attribute("FIELDS");
                std::vector<std::string> fieldStrings;

                if (fStr) // extract other fields.
                {
                    std::string fieldStr = fStr;
                    bool valid = ParseUtils::GenerateVector(
                        fieldStr.c_str(), fieldStrings);
                    ASSERTL0(valid, "Unable to correctly parse the field "
                                    "string in ExpansionTypes.");
                }

                // check to see if m_expasionVectorShPtrMap has
                // already been intiailised and if not intiailse
                // vector.
                if (m_expansionMapShPtrMap.count("DefaultVar") ==
                    0) // no previous definitions
                {
                    expansionMap = SetUpExpansionMap();

                    m_expansionMapShPtrMap["DefaultVar"] = expansionMap;

                    // make sure all fields in this search point
                    // to same expansion vector;
                    for (i = 0; i < fieldStrings.size(); ++i)
                    {
                        m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                    }
                }
                else // default variable is defined
                {

                    if (fieldStrings.size()) // fields are defined
                    {
                        // see if field exists
                        if (m_expansionMapShPtrMap.count(fieldStrings[0]))
                        {
                            expansionMap =
                                m_expansionMapShPtrMap.find(fieldStrings[0])
                                    ->second;
                        }
                        else
                        {
                            expansionMap = SetUpExpansionMap();
                            // make sure all fields in this search point
                            // to same expansion vector;
                            for (i = 0; i < fieldStrings.size(); ++i)
                            {
                                if (m_expansionMapShPtrMap.count(
                                        fieldStrings[i]) == 0)
                                {
                                    m_expansionMapShPtrMap[fieldStrings[i]] =
                                        expansionMap;
                                }
                                else
                                {
                                    ASSERTL0(false, "Expansion vector for this "
                                                    "field is already  setup");
                                }
                            }
                        }
                    }
                    else // use default variable list
                    {
                        expansionMap =
                            m_expansionMapShPtrMap.find("DefaultVar")->second;
                    }
                }

                /// Mandatory components...optional are to follow later.
                std::string compositeStr = expansion->Attribute("COMPOSITE");
                ASSERTL0(compositeStr.length() > 3,
                         "COMPOSITE must be specified in expansion definition");
                int beg = compositeStr.find_first_of("[");
                int end = compositeStr.find_first_of("]");
                std::string compositeListStr =
                    compositeStr.substr(beg + 1, end - beg - 1);

                map<int, CompositeSharedPtr> compositeVector;
                GetCompositeList(compositeListStr, compositeVector);

                bool useExpansionType = false;
                ExpansionType expansion_type;
                int num_modes;

                LibUtilities::BasisKeyVector basiskeyvec;
                const char *tStr = expansion->Attribute("TYPE");

                if (tStr) // use type string to define expansion
                {
                    std::string typeStr       = tStr;
                    const std::string *begStr = kExpansionTypeStr;
                    const std::string *endStr =
                        kExpansionTypeStr + eExpansionTypeSize;
                    const std::string *expStr =
                        std::find(begStr, endStr, typeStr);

                    ASSERTL0(expStr != endStr, "Invalid expansion type.");
                    expansion_type = (ExpansionType)(expStr - begStr);

                    /// \todo solvers break the pattern 'instantiate Session ->
                    /// instantiate MeshGraph'
                    /// and parse command line arguments by themselves; one
                    /// needs to unify command
                    /// line arguments handling.
                    /// Solvers tend to call MeshGraph::Read statically ->
                    /// m_session
                    /// is not defined -> no info about command line arguments
                    /// presented
                    /// ASSERTL0(m_session != 0, "One needs to instantiate
                    /// SessionReader first");

                    const char *nStr = expansion->Attribute("NUMMODES");
                    ASSERTL0(nStr, "NUMMODES was not defined in EXPANSION "
                                   "section of input");
                    std::string nummodesStr = nStr;

                    // ASSERTL0(m_session,"Session should be defined to evaluate
                    // nummodes ");
                    if (m_session)
                    {
                        LibUtilities::Equation nummodesEqn(m_session->GetInterpreter(), nummodesStr);
                        num_modes = (int)nummodesEqn.Evaluate();
                    }
                    else
                    {
                        num_modes = boost::lexical_cast<int>(nummodesStr);
                    }

                    useExpansionType = true;
                }
                else // assume expansion is defined individually
                {
                    // Extract the attributes.
                    const char *bTypeStr = expansion->Attribute("BASISTYPE");
                    ASSERTL0(bTypeStr, "TYPE or BASISTYPE was not defined in "
                                       "EXPANSION section of input");
                    std::string basisTypeStr = bTypeStr;

                    // interpret the basis type string.
                    std::vector<std::string> basisStrings;
                    std::vector<LibUtilities::BasisType> basis;
                    bool valid = ParseUtils::GenerateVector(
                        basisTypeStr.c_str(), basisStrings);
                    ASSERTL0(valid,
                             "Unable to correctly parse the basis types.");
                    for (vector<std::string>::size_type i = 0;
                         i < basisStrings.size(); i++)
                    {
                        valid = false;
                        for (unsigned int j = 0;
                             j < LibUtilities::SIZE_BasisType; j++)
                        {
                            if (LibUtilities::BasisTypeMap[j] ==
                                basisStrings[i])
                            {
                                basis.push_back((LibUtilities::BasisType)j);
                                valid = true;
                                break;
                            }
                        }
                        ASSERTL0(
                            valid,
                            std::string(
                                "Unable to correctly parse the basis type: ")
                                .append(basisStrings[i])
                                .c_str());
                    }
                    const char *nModesStr = expansion->Attribute("NUMMODES");
                    ASSERTL0(nModesStr, "NUMMODES was not defined in EXPANSION "
                                        "section of input");

                    std::string numModesStr = nModesStr;
                    std::vector<unsigned int> numModes;
                    valid = ParseUtils::GenerateVector(
                        numModesStr.c_str(), numModes);
                    ASSERTL0(valid,
                             "Unable to correctly parse the number of modes.");
                    ASSERTL0(numModes.size() == basis.size(),
                             "information for num modes does not match the "
                             "number of basis");

                    const char *pTypeStr = expansion->Attribute("POINTSTYPE");
                    ASSERTL0(pTypeStr, "POINTSTYPE was not defined in "
                                       "EXPANSION section of input");
                    std::string pointsTypeStr = pTypeStr;
                    // interpret the points type string.
                    std::vector<std::string> pointsStrings;
                    std::vector<LibUtilities::PointsType> points;
                    valid = ParseUtils::GenerateVector(
                        pointsTypeStr.c_str(), pointsStrings);
                    ASSERTL0(valid,
                             "Unable to correctly parse the points types.");
                    for (vector<std::string>::size_type i = 0;
                         i < pointsStrings.size(); i++)
                    {
                        valid = false;
                        for (unsigned int j = 0;
                             j < LibUtilities::SIZE_PointsType; j++)
                        {
                            if (LibUtilities::kPointsTypeStr[j] ==
                                pointsStrings[i])
                            {
                                points.push_back((LibUtilities::PointsType)j);
                                valid = true;
                                break;
                            }
                        }
                        ASSERTL0(
                            valid,
                            std::string(
                                "Unable to correctly parse the points type: ")
                                .append(pointsStrings[i])
                                .c_str());
                    }

                    const char *nPointsStr = expansion->Attribute("NUMPOINTS");
                    ASSERTL0(nPointsStr, "NUMPOINTS was not defined in "
                                         "EXPANSION section of input");
                    std::string numPointsStr = nPointsStr;
                    std::vector<unsigned int> numPoints;
                    valid = ParseUtils::GenerateVector(
                        numPointsStr.c_str(), numPoints);
                    ASSERTL0(valid,
                             "Unable to correctly parse the number of points.");
                    ASSERTL0(numPoints.size() == numPoints.size(),
                             "information for num points does not match the "
                             "number of basis");

                    for (int i = 0; i < basis.size(); ++i)
                    {
                        // Generate Basis key  using information
                        const LibUtilities::PointsKey pkey(numPoints[i],
                                                           points[i]);
                        basiskeyvec.push_back(LibUtilities::BasisKey(
                            basis[i], numModes[i], pkey));
                    }
                }

                // Now have composite and basiskeys.  Cycle through
                // all composites for the geomShPtrs and set the modes
                // and types for the elements contained in the element
                // list.
                for (auto compVecIter = compositeVector.begin();
                     compVecIter != compositeVector.end(); ++compVecIter)
                {
                    for (auto geomVecIter =
                             compVecIter->second->m_geomVec.begin();
                         geomVecIter != compVecIter->second->m_geomVec.end();
                         ++geomVecIter)
                    {
                        auto x =
                            expansionMap->find((*geomVecIter)->GetGlobalID());
                        ASSERTL0(x != expansionMap->end(),
                                 "Expansion not found!!");
                        if (useExpansionType)
                        {
                            (x->second)->m_basisKeyVector =
                                MeshGraph::DefineBasisKeyFromExpansionType(
                                    *geomVecIter, expansion_type, num_modes);
                        }
                        else
                        {
                            ASSERTL0((*geomVecIter)->GetShapeDim() ==
                                         basiskeyvec.size(),
                                     " There is an incompatible expansion "
                                     "dimension with geometry dimension");
                            (x->second)->m_basisKeyVector = basiskeyvec;
                        }
                    }
                }

                expansion = expansion->NextSiblingElement("E");
            }
        }
        else if (expType == "H")
        {
            int i;
            ExpansionMapShPtr expansionMap;

            while (expansion)
            {

                const char *fStr = expansion->Attribute("FIELDS");
                std::vector<std::string> fieldStrings;

                if (fStr) // extract other fields.
                {
                    std::string fieldStr = fStr;
                    bool valid = ParseUtils::GenerateVector(
                        fieldStr.c_str(), fieldStrings);
                    ASSERTL0(valid, "Unable to correctly parse the field "
                                    "string in ExpansionTypes.");
                }

                // check to see if m_expasionVectorShPtrMap has
                // already been intiailised and if not intiailse
                // vector.
                if (m_expansionMapShPtrMap.count("DefaultVar") ==
                    0) // no previous definitions
                {
                    expansionMap = SetUpExpansionMap();

                    m_expansionMapShPtrMap["DefaultVar"] = expansionMap;

                    // make sure all fields in this search point
                    // to same expansion vector;
                    for (i = 0; i < fieldStrings.size(); ++i)
                    {
                        m_expansionMapShPtrMap[fieldStrings[i]] = expansionMap;
                    }
                }
                else // default variable is defined
                {

                    if (fieldStrings.size()) // fields are defined
                    {
                        // see if field exists
                        if (m_expansionMapShPtrMap.count(fieldStrings[0]))
                        {
                            expansionMap =
                                m_expansionMapShPtrMap.find(fieldStrings[0])
                                    ->second;
                        }
                        else
                        {
                            expansionMap = SetUpExpansionMap();
                            // make sure all fields in this search point
                            // to same expansion vector;
                            for (i = 0; i < fieldStrings.size(); ++i)
                            {
                                if (m_expansionMapShPtrMap.count(
                                        fieldStrings[i]) == 0)
                                {
                                    m_expansionMapShPtrMap[fieldStrings[i]] =
                                        expansionMap;
                                }
                                else
                                {
                                    ASSERTL0(false, "Expansion vector for this "
                                                    "field is already  setup");
                                }
                            }
                        }
                    }
                    else // use default variable list
                    {
                        expansionMap =
                            m_expansionMapShPtrMap.find("DefaultVar")->second;
                    }
                }

                /// Mandatory components...optional are to follow later.
                std::string compositeStr = expansion->Attribute("COMPOSITE");
                ASSERTL0(compositeStr.length() > 3,
                         "COMPOSITE must be specified in expansion definition");
                int beg = compositeStr.find_first_of("[");
                int end = compositeStr.find_first_of("]");
                std::string compositeListStr =
                    compositeStr.substr(beg + 1, end - beg - 1);

                map<int, CompositeSharedPtr> compositeVector;
                GetCompositeList(compositeListStr, compositeVector);

                ExpansionType expansion_type_x = eNoExpansionType;
                ExpansionType expansion_type_y = eNoExpansionType;
                ExpansionType expansion_type_z = eNoExpansionType;
                int num_modes_x                = 0;
                int num_modes_y                = 0;
                int num_modes_z                = 0;

                LibUtilities::BasisKeyVector basiskeyvec;

                const char *tStr_x = expansion->Attribute("TYPE-X");

                if (tStr_x) // use type string to define expansion
                {
                    std::string typeStr       = tStr_x;
                    const std::string *begStr = kExpansionTypeStr;
                    const std::string *endStr =
                        kExpansionTypeStr + eExpansionTypeSize;
                    const std::string *expStr =
                        std::find(begStr, endStr, typeStr);

                    ASSERTL0(expStr != endStr, "Invalid expansion type.");
                    expansion_type_x = (ExpansionType)(expStr - begStr);

                    const char *nStr = expansion->Attribute("NUMMODES-X");
                    ASSERTL0(nStr, "NUMMODES-X was not defined in EXPANSION "
                                   "section of input");
                    std::string nummodesStr = nStr;

                    // ASSERTL0(m_session,"Session should be defined to evaluate
                    // nummodes ");

                    if (m_session)
                    {
                        LibUtilities::Equation nummodesEqn(m_session->GetInterpreter(), nummodesStr);
                        num_modes_x = (int)nummodesEqn.Evaluate();
                    }
                    else
                    {
                        num_modes_x = boost::lexical_cast<int>(nummodesStr);
                    }
                }

                const char *tStr_y = expansion->Attribute("TYPE-Y");

                if (tStr_y) // use type string to define expansion
                {
                    std::string typeStr       = tStr_y;
                    const std::string *begStr = kExpansionTypeStr;
                    const std::string *endStr =
                        kExpansionTypeStr + eExpansionTypeSize;
                    const std::string *expStr =
                        std::find(begStr, endStr, typeStr);

                    ASSERTL0(expStr != endStr, "Invalid expansion type.");
                    expansion_type_y = (ExpansionType)(expStr - begStr);

                    const char *nStr = expansion->Attribute("NUMMODES-Y");
                    ASSERTL0(nStr, "NUMMODES-Y was not defined in EXPANSION "
                                   "section of input");
                    std::string nummodesStr = nStr;

                    // ASSERTL0(m_session,"Session should be defined to evaluate
                    // nummodes ");
                    if (m_session)
                    {
                        LibUtilities::Equation nummodesEqn(m_session->GetInterpreter(), nummodesStr);
                        num_modes_y = (int)nummodesEqn.Evaluate();
                    }
                    else
                    {
                        num_modes_y = boost::lexical_cast<int>(nummodesStr);
                    }
                }

                const char *tStr_z = expansion->Attribute("TYPE-Z");

                if (tStr_z) // use type string to define expansion
                {
                    std::string typeStr       = tStr_z;
                    const std::string *begStr = kExpansionTypeStr;
                    const std::string *endStr =
                        kExpansionTypeStr + eExpansionTypeSize;
                    const std::string *expStr =
                        std::find(begStr, endStr, typeStr);

                    ASSERTL0(expStr != endStr, "Invalid expansion type.");
                    expansion_type_z = (ExpansionType)(expStr - begStr);

                    const char *nStr = expansion->Attribute("NUMMODES-Z");
                    ASSERTL0(nStr, "NUMMODES-Z was not defined in EXPANSION "
                                   "section of input");
                    std::string nummodesStr = nStr;

                    // ASSERTL0(m_session,"Session should be defined to evaluate
                    // nummodes ");
                    if (m_session)
                    {
                        LibUtilities::Equation nummodesEqn(m_session->GetInterpreter(), nummodesStr);
                        num_modes_z = (int)nummodesEqn.Evaluate();
                    }
                    else
                    {
                        num_modes_z = boost::lexical_cast<int>(nummodesStr);
                    }
                }

                for (auto compVecIter = compositeVector.begin();
                     compVecIter != compositeVector.end(); ++compVecIter)
                {
                    for (auto geomVecIter =
                             compVecIter->second->m_geomVec.begin();
                         geomVecIter != compVecIter->second->m_geomVec.end();
                         ++geomVecIter)
                    {
                        for (auto expVecIter = expansionMap->begin();
                             expVecIter != expansionMap->end(); ++expVecIter)
                        {

                            (expVecIter->second)->m_basisKeyVector =
                                DefineBasisKeyFromExpansionTypeHomo(
                                    *geomVecIter, expansion_type_x,
                                    expansion_type_y, expansion_type_z,
                                    num_modes_x, num_modes_y, num_modes_z);
                        }
                    }
                }

                expansion = expansion->NextSiblingElement("H");
            }
        }
        else if (expType ==
                 "ELEMENTS") // Reading a file with the expansion definition
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;

            // This has to use the XML reader since we are treating the already
            // parsed XML as a standard FLD file.
            std::shared_ptr<LibUtilities::FieldIOXml> f =
                make_shared<LibUtilities::FieldIOXml>(m_session->GetComm(), false);
            f->ImportFieldDefs(LibUtilities::XmlDataSource::create(m_session->GetDocument()),
                               fielddefs, true);
            cout << "    Number of elements: " << fielddefs.size() << endl;
            SetExpansions(fielddefs);
        }
        else if (expType == "F")
        {
            ASSERTL0(expansion->Attribute("FILE"),
                     "Attribute FILE expected for type F expansion");
            std::string filenameStr = expansion->Attribute("FILE");
            ASSERTL0(!filenameStr.empty(),
                     "A filename must be specified for the FILE "
                     "attribute of expansion");

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;
            LibUtilities::FieldIOSharedPtr f =
                LibUtilities::FieldIO::CreateForFile(m_session, filenameStr);
            f->Import(filenameStr, fielddefs);
            SetExpansions(fielddefs);
        }
        else
        {
            ASSERTL0(false, "Expansion type not defined");
        }
    }
}

GeometryLinkSharedPtr MeshGraph::GetElementsFromEdge(Geometry1DSharedPtr edge)
{
    // Search tris and quads
    // Need to iterate through vectors because there may be multiple
    // occurrences.

    GeometryLinkSharedPtr ret =
        GeometryLinkSharedPtr(new vector<pair<GeometrySharedPtr, int>>);

    TriGeomSharedPtr triGeomShPtr;
    QuadGeomSharedPtr quadGeomShPtr;

    for (int d = 0; d < m_domain.size(); ++d)
    {
        for (auto compIter = m_domain[d].begin(); compIter != m_domain[d].end();
             ++compIter)
        {
            for (auto &geomIter : compIter->second->m_geomVec)
            {
                triGeomShPtr  = std::dynamic_pointer_cast<TriGeom>(geomIter);
                quadGeomShPtr = std::dynamic_pointer_cast<QuadGeom>(geomIter);

                if (triGeomShPtr || quadGeomShPtr)
                {
                    if (triGeomShPtr)
                    {
                        for (int i = 0; i < triGeomShPtr->GetNumEdges(); i++)
                        {
                            if (triGeomShPtr->GetEdge(i)->GetGlobalID() ==
                                edge->GetGlobalID())
                            {
                                ret->push_back(make_pair(triGeomShPtr, i));
                                break;
                            }
                        }
                    }
                    else if (quadGeomShPtr)
                    {
                        for (int i = 0; i < quadGeomShPtr->GetNumEdges(); i++)
                        {
                            if (quadGeomShPtr->GetEdge(i)->GetGlobalID() ==
                                edge->GetGlobalID())
                            {
                                ret->push_back(make_pair(quadGeomShPtr, i));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    return ret;
}

LibUtilities::BasisKey MeshGraph::GetEdgeBasisKey(SegGeomSharedPtr edge,
                                                  const std::string variable)
{
    GeometryLinkSharedPtr elmts = GetElementsFromEdge(edge);
    // Perhaps, a check should be done here to ensure that
    // in case elements->size!=1, all elements to which
    // the edge belongs have the same type and order of
    // expansion such that no confusion can arise.
    GeometrySharedPtr geom   = elmts->at(0).first;
    ExpansionShPtr expansion = GetExpansion(geom, variable);
    int edge_id              = elmts->at(0).second;
    if (geom->GetShapeType() == LibUtilities::eTriangle)
    {
        edge_id = (edge_id) ? 1 : 0;
    }
    else
    {
        edge_id = edge_id % 2;
    }
    int nummodes  = expansion->m_basisKeyVector[edge_id].GetNumModes();
    int numpoints = expansion->m_basisKeyVector[edge_id].GetNumPoints();
    if (geom->GetShapeType() == LibUtilities::eTriangle)
    {
        // Use edge 0 to define basis of order relevant to edge
        switch (expansion->m_basisKeyVector[edge_id].GetBasisType())
        {
            case LibUtilities::eGLL_Lagrange:
            {
                switch (expansion->m_basisKeyVector[edge_id].GetPointsType())
                {
                    case LibUtilities::eGaussLobattoLegendre:
                    {
                        const LibUtilities::PointsKey pkey(
                            numpoints, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Unexpected points distribution");
                        // It doesn't matter what we return
                        // here since the ASSERT will stop
                        // execution.  Just return something
                        // to prevent warnings messages.
                        const LibUtilities::PointsKey pkey(
                            numpoints, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                        break;
                }
            }
            break;
            case LibUtilities::eOrtho_B: // Assume this is called from nodal
                                         // triangular basis
                {
                    switch (
                        expansion->m_basisKeyVector[edge_id].GetPointsType())
                    {
                        case LibUtilities::eGaussRadauMAlpha1Beta0:
                        {
                            const LibUtilities::PointsKey pkey(
                                numpoints + 1,
                                LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(
                                LibUtilities::eGLL_Lagrange, nummodes, pkey);
                        }
                        break;
                        default:
                            ASSERTL0(false, "Unexpected points distribution");
                            // It doesn't matter what we return
                            // here since the ASSERT will stop
                            // execution.  Just return something
                            // to prevent warnings messages.
                            const LibUtilities::PointsKey pkey(
                                numpoints + 1,
                                LibUtilities::eGaussLobattoLegendre);
                            return LibUtilities::BasisKey(
                                expansion->m_basisKeyVector[0].GetBasisType(),
                                nummodes, pkey);
                            break;
                    }
                }
                break;
            case LibUtilities::eModified_B:
            {
                switch (expansion->m_basisKeyVector[edge_id].GetPointsType())
                {
                    case LibUtilities::eGaussRadauMAlpha1Beta0:
                    {
                        const LibUtilities::PointsKey pkey(
                            numpoints + 1, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                    }
                    break;
                    case LibUtilities::eGaussLobattoLegendre:
                    {
                        const LibUtilities::PointsKey pkey(
                            numpoints, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Unexpected points distribution");
                        // It doesn't matter what we return
                        // here since the ASSERT will stop
                        // execution.  Just return something
                        // to prevent warnings messages.
                        const LibUtilities::PointsKey pkey(
                            numpoints + 1, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                        break;
                }
            }
            break;
            case LibUtilities::eModified_A:
            {
                switch (expansion->m_basisKeyVector[edge_id].GetPointsType())
                {
                    case LibUtilities::eGaussLobattoLegendre:
                    {
                        const LibUtilities::PointsKey pkey(
                            numpoints, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                    }
                    break;
                    default:
                        ASSERTL0(false, "Unexpected points distribution");
                        // It doesn't matter what we return here
                        // since the ASSERT will stop execution.
                        // Just return something to prevent
                        // warnings messages.
                        const LibUtilities::PointsKey pkey(
                            numpoints, LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                            expansion->m_basisKeyVector[0].GetBasisType(),
                            nummodes, pkey);
                        break;
                }
            }
            break;
            default:
                ASSERTL0(false, "Unexpected basis distribution");
                // It doesn't matter what we return here since the
                // ASSERT will stop execution.  Just return
                // something to prevent warnings messages.
                const LibUtilities::PointsKey pkey(
                    numpoints + 1, LibUtilities::eGaussLobattoLegendre);
                return LibUtilities::BasisKey(
                    expansion->m_basisKeyVector[0].GetBasisType(), nummodes,
                    pkey);
        }
    }
    else
    {
        // Quadrilateral
        const LibUtilities::PointsKey pkey(
            numpoints, expansion->m_basisKeyVector[edge_id].GetPointsType());
        return LibUtilities::BasisKey(
            expansion->m_basisKeyVector[edge_id].GetBasisType(), nummodes,
            pkey);
    }

    ASSERTL0(false, "Unable to determine edge points type.");
    return LibUtilities::NullBasisKey;
}

/// 3D functions

LibUtilities::BasisKey MeshGraph::GetFaceBasisKey(Geometry2DSharedPtr face,
                                                  const int facedir,
                                                  const std::string variable)
{
    // Retrieve the list of elements and the associated face index
    // to which the face geometry belongs.
    GeometryLinkSharedPtr elements = GetElementsFromFace(face);
    ASSERTL0(elements->size() > 0,
             "No elements for the given face."
             " Check all elements belong to the domain composite.");
    // Perhaps, a check should be done here to ensure that in case
    // elements->size!=1, all elements to which the edge belongs have
    // the same type and order of expansion such that no confusion can
    // arise.
    // Get the Expansion structure detailing the basis keys used for
    // this element.
    GeometrySharedPtr geom   = elements->at(0).first;
    ExpansionShPtr expansion = GetExpansion(geom, variable);
    ASSERTL0(expansion, "Could not find expansion connected to face " +
                            boost::lexical_cast<string>(face->GetGlobalID()));
    // Retrieve the geometry object of the element as a Geometry3D.
    Geometry3DSharedPtr geom3d =
        std::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
            expansion->m_geomShPtr);
    // Use the geometry of the element to calculate the coordinate
    // direction of the element which corresponds to the requested
    // coordinate direction of the given face.
    int dir = geom3d->GetDir(elements->at(0).second, facedir);
    if (face->GetNumVerts() == 3)
    {
        return StdRegions::EvaluateTriFaceBasisKey(
            facedir, expansion->m_basisKeyVector[dir].GetBasisType(),
            expansion->m_basisKeyVector[dir].GetNumPoints(),
            expansion->m_basisKeyVector[dir].GetNumModes());
    }
    else
    {
        return StdRegions::EvaluateQuadFaceBasisKey(
            facedir, expansion->m_basisKeyVector[dir].GetBasisType(),
            expansion->m_basisKeyVector[dir].GetNumPoints(),
            expansion->m_basisKeyVector[dir].GetNumModes());
    }

    // Keep things happy by returning a value.
    return LibUtilities::NullBasisKey;
}

GeometryLinkSharedPtr MeshGraph::GetElementsFromFace(Geometry2DSharedPtr face)
{
    auto it = m_faceToElMap.find(face->GetGlobalID());

    ASSERTL0(it != m_faceToElMap.end(), "Unable to find corresponding face!");

    return it->second;
}

/**
 * @brief Given a 3D geometry object #element, populate the face to
 * element map #m_faceToElMap which maps faces to their corresponding
 * element(s).
 *
 * @param element  Element to process.
 * @param kNfaces  Number of faces of #element. Should be removed and
 * put into Geometry3D as a virtual member function.
 */
void MeshGraph::PopulateFaceToElMap(Geometry3DSharedPtr element, int kNfaces)
{
    // Set up face -> element map
    for (int i = 0; i < kNfaces; ++i)
    {
        int faceId = element->GetFace(i)->GetGlobalID();

        // Search map to see if face already exists.
        auto it = m_faceToElMap.find(faceId);

        if (it == m_faceToElMap.end())
        {
            GeometryLinkSharedPtr tmp =
                GeometryLinkSharedPtr(new vector<pair<GeometrySharedPtr, int>>);
            tmp->push_back(make_pair(element, i));
            m_faceToElMap[faceId] = tmp;
        }
        else
        {
            it->second->push_back(make_pair(element, i));
        }
    }
}

/**
 * @brief Create mesh entities for this graph.
 *
 * This function will create a map of all mesh entities of the current graph,
 * which can then be used within the mesh partitioner to construct an
 * appropriate partitioning.
 */
std::map<int, MeshEntity> MeshGraph::CreateMeshEntities()
{
    std::map<int, MeshEntity> elements;
    switch (m_meshDimension)
    {
        case 1:
        {
            for (auto &i : m_segGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetVertex(0)->GetGlobalID());
                e.list.push_back(i.second->GetVertex(1)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
        }
        break;
        case 2:
        {
            for (auto &i : m_triGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetEdge(0)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(1)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(2)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
            for (auto &i : m_quadGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetEdge(0)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(1)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(2)->GetGlobalID());
                e.list.push_back(i.second->GetEdge(3)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
        }
        break;
        case 3:
        {
            for (auto &i : m_tetGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
            for (auto &i : m_pyrGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
            for (auto &i : m_prismGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
            for (auto &i : m_hexGeoms)
            {
                MeshEntity e;
                e.id = e.origId = i.first;
                e.list.push_back(i.second->GetFace(0)->GetGlobalID());
                e.list.push_back(i.second->GetFace(1)->GetGlobalID());
                e.list.push_back(i.second->GetFace(2)->GetGlobalID());
                e.list.push_back(i.second->GetFace(3)->GetGlobalID());
                e.list.push_back(i.second->GetFace(4)->GetGlobalID());
                e.list.push_back(i.second->GetFace(5)->GetGlobalID());
                e.ghost = false;
                elements[e.id] = e;
            }
        }
        break;
    }

    return elements;
}

CompositeDescriptor MeshGraph::CreateCompositeDescriptor()
{
    CompositeDescriptor ret;

    for (auto &comp : m_meshComposites)
    {
        std::pair<LibUtilities::ShapeType, vector<int>> tmp;
        tmp.first = comp.second->m_geomVec[0]->GetShapeType();

        tmp.second.resize(comp.second->m_geomVec.size());
        for (size_t i = 0; i < tmp.second.size(); ++i)
        {
            tmp.second[i] = comp.second->m_geomVec[i]->GetGlobalID();
        }

        ret[comp.first] = tmp;
    }

    return ret;
}



} // end of namespace
} // end of namespace
