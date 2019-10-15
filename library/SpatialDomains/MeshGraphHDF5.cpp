////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshGraphHDF5.cpp
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
//  Description: HDF5-based mesh format for Nektar++.
//
////////////////////////////////////////////////////////////////////////////////

#include <type_traits>
#include <tinyxml.h>
#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshPartition.h>
#include <SpatialDomains/MeshGraphHDF5.h>

#define TIME_RESULT(verb, msg, timer)                            \
    if (verb)                                                    \
    {                                                            \
        std::cout << "  - " << msg << ": "                       \
                  << timer.TimePerTest(1) << std::endl;          \
    }

using namespace std;
using namespace Nektar::LibUtilities;

namespace Nektar
{
namespace SpatialDomains
{

/// Version of the Nektar++ HDF5 geometry format, which is embedded into the
/// main NEKTAR/GEOMETRY group as an attribute.
const unsigned int MeshGraphHDF5::FORMAT_VERSION = 1;

std::string MeshGraphHDF5::className =
    GetMeshGraphFactory().RegisterCreatorFunction(
        "HDF5", MeshGraphHDF5::create, "IO with HDF5 geometry");

void MeshGraphHDF5::ReadGeometry(
    DomainRangeShPtr rng,
    bool fillGraph)
{
    boost::ignore_unused(rng);

    ReadComposites();
    ReadDomain();
    ReadExpansions();

    // Close up shop.
    m_mesh->Close();
    m_maps->Close();
    m_file->Close();

    if (fillGraph)
    {
        MeshGraph::FillGraph();
    }
}

/**
 * @brief Utility function to split a vector equally amongst a number of
 * processors.
 *
 * @param vecsize  Size of the total amount of work
 * @param rank     Rank of this process
 * @param nprocs   Number of processors in the group
 *
 * @return A pair with the offset this process should occupy, along with the
 *         count for the amount of work.
 */
std::pair<size_t, size_t> SplitWork(size_t vecsize, int rank, int nprocs)
{
    size_t div = vecsize / nprocs;
    size_t rem = vecsize % nprocs;
    if (rank < rem)
    {
        return std::make_pair(rank * (div + 1), div + 1);
    }
    else
    {
        return std::make_pair((rank - rem) * div + rem * (div + 1), div);
    }
}

template<class T, typename std::enable_if<T::kDim == 0, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap)
{
    boost::ignore_unused(geomMap);
    return 3;
}

template<class T, typename std::enable_if<T::kDim == 1, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap)
{
    boost::ignore_unused(geomMap);
    return T::kNverts;
}

template<class T, typename std::enable_if<T::kDim == 2, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap)
{
    boost::ignore_unused(geomMap);
    return T::kNedges;
}

template<class T, typename std::enable_if<T::kDim == 3, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap)
{
    boost::ignore_unused(geomMap);
    return T::kNfaces;
}

template<class ...T>
inline void UniqueValues(std::unordered_set<int> &unique)
{
    boost::ignore_unused(unique);
}

template<class ...T>
inline void UniqueValues(std::unordered_set<int> &unique,
                         const std::vector<int> &input,
                         T&... args)
{
    for (auto i : input)
    {
        unique.insert(i);
    }

    UniqueValues(unique, args...);
}

/**
 * @brief Partition the mesh
 */
void MeshGraphHDF5::PartitionMesh(LibUtilities::SessionReaderSharedPtr session)
{
    LibUtilities::Timer all;
    all.Start();
    int err;
    LibUtilities::CommSharedPtr comm = session->GetComm();
    int rank = comm->GetRank(), nproc = comm->GetSize();

    // By default, only the root process will have read the session file, which
    // is done to avoid every process needing to read the XML file. For HDF5, we
    // don't care about this, so just have every process parse the session file.
    if (rank > 0)
    {
        session->InitSession();
    }

    // We use the XML geometry to find information about the HDF5 file.
    m_session            = session;
    m_xmlGeom            = m_session->GetElement("NEKTAR/GEOMETRY");
    TiXmlAttribute *attr = m_xmlGeom->FirstAttribute();
    m_meshPartitioned    = true;
    m_meshDimension      = 3;
    m_spaceDimension     = 3;

    while (attr)
    {
        std::string attrName(attr->Name());
        if (attrName == "DIM")
        {
            err = attr->QueryIntValue(&m_meshDimension);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read mesh dimension.");
        }
        else if (attrName == "SPACE")
        {
            err = attr->QueryIntValue(&m_spaceDimension);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read space dimension.");
        }
        else if (attrName == "PARTITION")
        {
            ASSERTL0(false,
                     "PARTITION parameter should only be used in XML meshes");
        }
        else if(attrName == "HDF5FILE")
        {
            m_hdf5Name = attr->Value();
            ASSERTL1(err == TIXML_SUCCESS, "Unable to read hdf5 name.");
        }
        else if(attrName == "PARTITIONED")
        {
            ASSERTL0(false,
                     "PARTITIONED parameter should only be used in XML meshes");
        }
        else
        {
            std::string errstr("Unknown attribute: ");
            errstr += attrName;
            ASSERTL1(false, errstr.c_str());
        }
        // Get the next attribute.
        attr = attr->Next();
    }

    ASSERTL0(m_hdf5Name.size() > 0, "Unable to obtain mesh file name.");
    ASSERTL0(m_meshDimension <= m_spaceDimension,
             "Mesh dimension greater than space dimension.");

    // Open handle to the HDF5 mesh
    LibUtilities::H5::PListSharedPtr parallelProps = H5::PList::Default();
    m_readPL = H5::PList::Default();

    if (nproc > 1)
    {
        // Use MPI/O to access the file
        parallelProps = H5::PList::FileAccess();
        parallelProps->SetMpio(comm);
        // Use collective IO
        m_readPL = H5::PList::DatasetXfer();
        m_readPL->SetDxMpioCollective();
    }

    m_file = H5::File::Open(m_hdf5Name, H5F_ACC_RDONLY, parallelProps);

    auto root = m_file->OpenGroup("NEKTAR");
    ASSERTL0(root, "Cannot find NEKTAR group in HDF5 file.");

    auto root2 = root->OpenGroup("GEOMETRY");
    ASSERTL0(root2, "Cannot find NEKTAR/GEOMETRY group in HDF5 file.");

    // Check format version
    unsigned int formatVersion;
    H5::Group::AttrIterator attrIt  = root2->attr_begin();
    H5::Group::AttrIterator attrEnd = root2->attr_end();
    for (; attrIt != attrEnd; ++attrIt)
    {
        if (*attrIt == "FORMAT_VERSION")
        {
            break;
        }
    }
    ASSERTL0(attrIt != attrEnd,
             "Unable to determine Nektar++ geometry HDF5 file version.");
    root2->GetAttribute("FORMAT_VERSION", formatVersion);

    ASSERTL0(formatVersion <= FORMAT_VERSION,
             "File format in " + m_hdf5Name + " is higher than supported in "
             "this version of Nektar++");

    m_mesh = root2->OpenGroup("MESH");
    ASSERTL0(m_mesh, "Cannot find NEKTAR/GEOMETRY/MESH group in HDF5 file.");
    m_maps = root2->OpenGroup("MAPS");
    ASSERTL0(m_mesh, "Cannot find NEKTAR/GEOMETRY/MAPS group in HDF5 file.");

    // Depending on dimension, read element IDs.
    std::map<int, std::vector<std::tuple<
        std::string, int, LibUtilities::ShapeType>>> dataSets;

    dataSets[1] = { make_tuple("SEG", 2, LibUtilities::eSegment) };
    dataSets[2] = { make_tuple("TRI", 3, LibUtilities::eTriangle),
                    make_tuple("QUAD", 4, LibUtilities::eQuadrilateral) };
    dataSets[3] = { make_tuple("TET", 4, LibUtilities::eTetrahedron),
                    make_tuple("PYR", 5, LibUtilities::ePyramid),
                    make_tuple("PRISM", 5, LibUtilities::ePrism),
                    make_tuple("HEX", 6, LibUtilities::eHexahedron) };

    // Read IDs for partitioning purposes
    std::vector<MeshEntity> elmts;
    std::vector<int> ids;

    // Map from element ID to 'row' which is a contiguous ordering required for
    // parallel partitioning.
    std::unordered_map<int, int> row2id, id2row;

    LibUtilities::Timer t;
    t.Start();
    int rowCount = 0;
    for (auto &it : dataSets[m_meshDimension])
    {
        std::string ds = std::get<0>(it);

        if (!m_mesh->ContainsDataSet(ds))
        {
            continue;
        }

        // Open metadata dataset
        H5::DataSetSharedPtr data = m_mesh->OpenDataSet(ds);
        H5::DataSpaceSharedPtr space = data->GetSpace();
        vector<hsize_t> dims = space->GetDims();

        H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(ds);
        H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
        vector<hsize_t> mdims = mspace->GetDims();

        // TODO: This could be done more intelligently; reads all IDs so that we
        // can construct the dual graph of the mesh.
        vector<int> tmpElmts, tmpIds;
        mdata->Read(tmpIds, mspace, m_readPL);
        data->Read(tmpElmts, space, m_readPL);

        const int nGeomData = std::get<1>(it);

        for (int i = 0, cnt = 0; i < tmpIds.size(); ++i, ++rowCount)
        {
            MeshEntity e;
            row2id[rowCount] = tmpIds[i];
            id2row[tmpIds[i]] = row2id[rowCount];
            e.id = rowCount;
            e.origId = tmpIds[i];
            e.list = std::vector<unsigned int>(
                &tmpElmts[cnt], &tmpElmts[cnt+nGeomData]);
            elmts.push_back(e);
            cnt += nGeomData;
        }
    }

    comm->Block();
    t.Stop();

    bool verbRoot = rank == 0 &&
        m_session->DefinesCmdLineArgument("verbose");

    if (verbRoot)
    {
        std::cout << "Reading HDF5 geometry..." << std::endl;
    }

    TIME_RESULT(verbRoot, "initial read", t);

    // Check to see we have at least as many processors as elements.
    size_t numElmt = elmts.size();
    ASSERTL0(nproc <= numElmt,
             "This mesh has more processors than elements!");

    // Calculate reasonably even distribution of processors for calling ptScotch
    auto elRange = SplitWork(numElmt, rank, nproc);

    t.Start();

    // Construct map of element entities for partitioner.
    std::map<int, MeshEntity> partElmts;
    std::unordered_set<int> facetIDs;

    int vcnt = 0;

    for (int el = elRange.first; el < elRange.first + elRange.second;
         ++el, ++vcnt)
    {
        MeshEntity elmt = elmts[el];
        elmt.ghost = false;
        partElmts[el] = elmt;

        for (auto &facet : elmt.list)
        {
            facetIDs.insert(facet);
        }
    }

    // Now identify ghost vertices for the graph. This could probably be
    // improved.
    int nLocal = vcnt;
    for (int i = 0; i < numElmt; ++i)
    {
        // Ignore anything we already read.
        if (i >= elRange.first && i < elRange.first + elRange.second)
        {
            //i += elRange.second - elRange.first - 1;
            continue;
        }

        MeshEntity elmt = elmts[i];
        bool insert = false;

        // Check for connections to local elements.
        for (auto &eId : elmt.list)
        {
            if (facetIDs.find(eId) != facetIDs.end())
            {
                insert = true;
                break;
            }
        }

        if (insert)
        {
            elmt.ghost = true;
            partElmts[elmt.id] = elmt;
        }
    }

    // Create partitioner. Default partitioner to use is PtScotch. Use ParMetis
    // as default if it is installed. Override default with command-line flags
    // if they are set.
    string partitionerName = nproc > 1 ? "PtScotch" : "Scotch";
    if (GetMeshPartitionFactory().ModuleExists("ParMetis"))
    {
        partitionerName = "ParMetis";
    }
    if (session->DefinesCmdLineArgument("use-parmetis"))
    {
        partitionerName = "ParMetis";
    }
    if (session->DefinesCmdLineArgument("use-ptscotch"))
    {
        partitionerName = "PtScotch";
    }

    MeshPartitionSharedPtr partitioner =
        GetMeshPartitionFactory().CreateInstance(
            partitionerName, session, m_meshDimension,
            partElmts, CreateCompositeDescriptor(id2row));

    partitioner->PartitionMesh(nproc, true, false, nLocal);

    // Now we need to distribute vertex IDs to all the different processors.
    t.Stop();
    TIME_RESULT(verbRoot, "partitioning", t);

    // Each process now knows which rows of the dataset it needs to read for the
    // elements of dimension m_meshDimension.
    std::vector<unsigned int> tmp;
    std::unordered_set<int> toRead;
    partitioner->GetElementIDs(comm->GetRank(), tmp);

    for (auto &tmpId : tmp)
    {
        toRead.insert(row2id[tmpId]);
    }

    // Since objects are going to be constructed starting from vertices, we now
    // need to recurse down the geometry facet dimensions to figure out which
    // rows to read from each dataset.
    std::vector<int> vertIDs, segIDs, triIDs, quadIDs;
    std::vector<int> tetIDs, prismIDs, pyrIDs, hexIDs;
    std::vector<int> segData, triData, quadData, tetData;
    std::vector<int> prismData, pyrData, hexData;
    std::vector<NekDouble> vertData;

    if (m_meshDimension == 3)
    {
        t.Start();
        // Read 3D data
        ReadGeometryData(m_hexGeoms, "HEX", toRead, hexIDs, hexData);
        ReadGeometryData(m_pyrGeoms, "PYR", toRead, pyrIDs, pyrData);
        ReadGeometryData(m_prismGeoms, "PRISM", toRead, prismIDs, prismData);
        ReadGeometryData(m_tetGeoms, "TET", toRead, tetIDs, tetData);

        toRead.clear();
        UniqueValues(toRead, hexData, pyrData, prismData, tetData);
        t.Stop();
        TIME_RESULT(verbRoot, "read 3D elements", t);
    }

    if (m_meshDimension >= 2)
    {
        t.Start();
        // Read 2D data
        ReadGeometryData(m_triGeoms, "TRI", toRead, triIDs, triData);
        ReadGeometryData(m_quadGeoms, "QUAD", toRead, quadIDs, quadData);

        toRead.clear();
        UniqueValues(toRead, triData, quadData);
        t.Stop();
        TIME_RESULT(verbRoot, "read 2D elements", t);
    }

    if (m_meshDimension >= 1)
    {
        t.Start();
        // Read 2D data
        ReadGeometryData(m_segGeoms, "SEG", toRead, segIDs, segData);

        toRead.clear();
        UniqueValues(toRead, segData);
        t.Stop();
        TIME_RESULT(verbRoot, "read 1D elements", t);
    }

    t.Start();
    ReadGeometryData(m_vertSet, "VERT", toRead, vertIDs, vertData);
    t.Stop();
    TIME_RESULT(verbRoot, "read 0D elements", t);

    // Now start to construct geometry objects, starting from vertices upwards.
    t.Start();
    FillGeomMap(m_vertSet, CurveMap(), vertIDs, vertData);
    t.Stop();
    TIME_RESULT(verbRoot, "construct 0D elements", t);

    if (m_meshDimension >= 1)
    {
        // Read curves
        toRead.clear();
        for (auto &edge : segIDs)
        {
            toRead.insert(edge);
        }
        ReadCurveMap(m_curvedEdges, "CURVE_EDGE", toRead);

        t.Start();
        FillGeomMap(m_segGeoms, m_curvedEdges, segIDs, segData);
        t.Stop();
        TIME_RESULT(verbRoot, "construct 1D elements", t);
    }

    if (m_meshDimension >= 2)
    {
        // Read curves
        toRead.clear();
        for (auto &face : triIDs)
        {
            toRead.insert(face);
        }
        for (auto &face : quadIDs)
        {
            toRead.insert(face);
        }
        ReadCurveMap(m_curvedFaces, "CURVE_FACE", toRead);

        t.Start();
        FillGeomMap(m_triGeoms, m_curvedFaces, triIDs, triData);
        FillGeomMap(m_quadGeoms, m_curvedFaces, quadIDs, quadData);
        t.Stop();
        TIME_RESULT(verbRoot, "construct 2D elements", t);
    }

    if (m_meshDimension >= 3)
    {
        t.Start();
        FillGeomMap(m_hexGeoms, CurveMap(), hexIDs, hexData);
        FillGeomMap(m_prismGeoms, CurveMap(), prismIDs, prismData);
        FillGeomMap(m_pyrGeoms, CurveMap(), pyrIDs, pyrData);
        FillGeomMap(m_tetGeoms, CurveMap(), tetIDs, tetData);
        t.Stop();
        TIME_RESULT(verbRoot, "construct 3D elements", t);
    }
    all.Stop();
    TIME_RESULT(verbRoot, "total time", all);
}

template<class T, typename DataType> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<T>> &geomMap, int id,
    DataType *data, CurveSharedPtr curve)
{
    boost::ignore_unused(geomMap, id, data, curve);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PointGeom>> &geomMap, int id,
    NekDouble *data, CurveSharedPtr curve)
{
    boost::ignore_unused(curve);
    geomMap[id] = MemoryManager<PointGeom>::AllocateSharedPtr(
        m_spaceDimension, id, data[0], data[1], data[2]);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<SegGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    PointGeomSharedPtr pts[2] = { GetVertex(data[0]), GetVertex(data[1]) };
    geomMap[id] = MemoryManager<SegGeom>::AllocateSharedPtr(
        id, m_spaceDimension, pts, curve);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<TriGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    SegGeomSharedPtr segs[3] = {
        GetSegGeom(data[0]), GetSegGeom(data[1]), GetSegGeom(data[2]) };
    geomMap[id] = MemoryManager<TriGeom>::AllocateSharedPtr(id, segs, curve);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<QuadGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    SegGeomSharedPtr segs[4] = {
        GetSegGeom(data[0]), GetSegGeom(data[1]), GetSegGeom(data[2]),
        GetSegGeom(data[3])
    };
    geomMap[id] = MemoryManager<QuadGeom>::AllocateSharedPtr(id, segs, curve);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<TetGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    boost::ignore_unused(curve);
    TriGeomSharedPtr faces[4] = {
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[0])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[1])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[2])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[3]))
    };

    auto tetGeom = MemoryManager<TetGeom>::AllocateSharedPtr(id, faces);
    PopulateFaceToElMap(tetGeom, TetGeom::kNfaces);
    geomMap[id] = tetGeom;
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PyrGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    boost::ignore_unused(curve);
    Geometry2DSharedPtr faces[5] = {
        GetGeometry2D(data[0]), GetGeometry2D(data[1]), GetGeometry2D(data[2]),
        GetGeometry2D(data[3]), GetGeometry2D(data[4])
    };

    auto pyrGeom = MemoryManager<PyrGeom>::AllocateSharedPtr(id, faces);
    PopulateFaceToElMap(pyrGeom, PyrGeom::kNfaces);
    geomMap[id] = pyrGeom;
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PrismGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    boost::ignore_unused(curve);
    Geometry2DSharedPtr faces[5] = {
        GetGeometry2D(data[0]), GetGeometry2D(data[1]), GetGeometry2D(data[2]),
        GetGeometry2D(data[3]), GetGeometry2D(data[4])
    };

    auto prismGeom = MemoryManager<PrismGeom>::AllocateSharedPtr(id, faces);
    PopulateFaceToElMap(prismGeom, PrismGeom::kNfaces);
    geomMap[id] = prismGeom;
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<HexGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    boost::ignore_unused(curve);
    QuadGeomSharedPtr faces[6] = {
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[0])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[1])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[2])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[3])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[4])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[5]))
    };

    auto hexGeom = MemoryManager<HexGeom>::AllocateSharedPtr(id, faces);
    PopulateFaceToElMap(hexGeom, HexGeom::kNfaces);
    geomMap[id] = hexGeom;
}

template<class T, typename DataType>
void MeshGraphHDF5::FillGeomMap(
    std::map<int, std::shared_ptr<T>> &geomMap,
    const CurveMap                    &curveMap,
    std::vector<int>                  &ids,
    std::vector<DataType>             &geomData)
{
    const int nGeomData = GetGeomDataDim(geomMap);
    const int nRows = geomData.size() / nGeomData;
    CurveSharedPtr empty;

    // Construct geometry object.
    if (curveMap.size() > 0)
    {
        for(int i = 0, cnt = 0; i < nRows; i++, cnt += nGeomData)
        {
            auto cIt = curveMap.find(ids[i]);
            ConstructGeomObject(
                geomMap, ids[i], &geomData[cnt],
                cIt == curveMap.end() ? empty : cIt->second);
        }
    }
    else
    {
        for(int i = 0, cnt = 0; i < nRows; i++, cnt += nGeomData)
        {
            ConstructGeomObject(geomMap, ids[i], &geomData[cnt], empty);
        }
    }
}

template<class T, typename DataType>
void MeshGraphHDF5::ReadGeometryData(
    std::map<int, std::shared_ptr<T>>      &geomMap,
    std::string                             dataSet,
    const std::unordered_set<int>          &readIds,
    std::vector<int>                       &ids,
    std::vector<DataType>                  &geomData)
{
    if (!m_mesh->ContainsDataSet(dataSet))
    {
        return;
    }

    // Open mesh dataset
    H5::DataSetSharedPtr data = m_mesh->OpenDataSet(dataSet);
    H5::DataSpaceSharedPtr space = data->GetSpace();
    vector<hsize_t> dims = space->GetDims();

    // Open metadata dataset
    H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(dataSet);
    H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
    vector<hsize_t> mdims = mspace->GetDims();

    ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

    const int nGeomData = GetGeomDataDim(geomMap);

    // Read all IDs
    vector<int> allIds;
    mdata->Read(allIds, mspace);

    // Selective reading; clear data space range so that we can select certain
    // rows from the datasets.
    space->ClearRange();

    int i = 0;
    std::vector<hsize_t> coords;
    for (auto &id : allIds)
    {
        if (readIds.find(id) != readIds.end())
        {
            for (int j = 0; j < nGeomData; ++j)
            {
                coords.push_back(i);
                coords.push_back(j);
            }
            ids.push_back(id);
        }
        ++i;
    }

    space->SetSelection(coords.size() / 2, coords);

    // Read selected data.
    data->Read(geomData, space, m_readPL);
}

void MeshGraphHDF5::ReadCurveMap(
    CurveMap                      &curveMap,
    std::string                    dsName,
    const std::unordered_set<int> &readIds)
{
    // If dataset does not exist, exit.
    if (!m_mesh->ContainsDataSet(dsName))
    {
        return;
    }

    // Open up curve map data.
    H5::DataSetSharedPtr curveData = m_mesh->OpenDataSet(dsName);
    H5::DataSpaceSharedPtr curveSpace = curveData->GetSpace();

    // Open up ID data set.
    H5::DataSetSharedPtr idData = m_maps->OpenDataSet(dsName);
    H5::DataSpaceSharedPtr idSpace = idData->GetSpace();

    // Read all IDs and clear data space.
    vector<int> ids, newIds;
    idData->Read(ids, idSpace);
    curveSpace->ClearRange();

    // Search IDs to figure out which curves to read.
    vector<hsize_t> curveSel;

    int cnt = 0;
    for (auto &id : ids)
    {
        if (readIds.find(id) != readIds.end())
        {
            curveSel.push_back(cnt);
            curveSel.push_back(0);
            curveSel.push_back(cnt);
            curveSel.push_back(1);
            curveSel.push_back(cnt);
            curveSel.push_back(2);
            newIds.push_back(id);
        }

        ++cnt;
    }

    // Check to see whether any processor will read anything
    auto toRead = newIds.size();
    m_session->GetComm()->AllReduce(toRead, LibUtilities::ReduceSum);

    if (toRead == 0)
    {
        return;
    }

    // Now read curve map and read data.
    vector<int> curveInfo;
    curveSpace->SetSelection(curveSel.size() / 2, curveSel);
    curveData->Read(curveInfo, curveSpace, m_readPL);

    curveSel.clear();

    std::unordered_map<int, int> curvePtOffset;

    // Construct curves. We'll populate nodes in a minute!
    for (int i = 0, cnt = 0, cnt2 = 0; i < curveInfo.size() / 3; ++i, cnt += 3)
    {
        CurveSharedPtr curve = MemoryManager<Curve>::AllocateSharedPtr(
            newIds[i], (LibUtilities::PointsType)curveInfo[cnt + 1]);

        curve->m_points.resize(curveInfo[cnt]);

        const int ptOffset = curveInfo[cnt + 2];

        for (int j = 0; j < curveInfo[cnt]; ++j)
        {
            // ptoffset gives us the row, multiply by 3 for number of
            // coordinates.
            curveSel.push_back(ptOffset + j);
            curveSel.push_back(0);
            curveSel.push_back(ptOffset + j);
            curveSel.push_back(1);
            curveSel.push_back(ptOffset + j);
            curveSel.push_back(2);
        }

        // Store the offset so we know to come back later on to fill in these
        // points.
        curvePtOffset[newIds[i]] = 3 * cnt2;
        cnt2 += curveInfo[cnt];

        curveMap[newIds[i]] = curve;
    }

    curveInfo.clear();

    // Open node data spacee.
    H5::DataSetSharedPtr nodeData = m_mesh->OpenDataSet("CURVE_NODES");
    H5::DataSpaceSharedPtr nodeSpace = nodeData->GetSpace();

    nodeSpace->ClearRange();
    nodeSpace->SetSelection(curveSel.size() / 2, curveSel);

    vector<NekDouble> nodeRawData;
    nodeData->Read(nodeRawData, nodeSpace, m_readPL);

    // Go back and populate data from nodes.
    for (auto &cIt : curvePtOffset)
    {
        CurveSharedPtr curve = curveMap[cIt.first];

        // Create nodes.
        int cnt = cIt.second;
        for (int i = 0; i < curve->m_points.size(); ++i, cnt += 3)
        {
            curve->m_points[i] = MemoryManager<PointGeom>::AllocateSharedPtr(
                0, m_spaceDimension, nodeRawData[cnt], nodeRawData[cnt+1],
                nodeRawData[cnt+2]);
        }
    }
}

void MeshGraphHDF5::ReadDomain()
{
    map<int, CompositeSharedPtr> fullDomain;
    H5::DataSetSharedPtr dst = m_mesh->OpenDataSet("DOMAIN");
    H5::DataSpaceSharedPtr space = dst->GetSpace();

    vector<string> data;
    dst->ReadVectorString(data, space, m_readPL);
    GetCompositeList(data[0], fullDomain);
    m_domain.push_back(fullDomain);
}

void MeshGraphHDF5::ReadComposites()
{
    string nm = "COMPOSITE";

    H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
    H5::DataSpaceSharedPtr space = data->GetSpace();
    vector<hsize_t> dims = space->GetDims();

    vector<string> comps;
    data->ReadVectorString(comps, space);

    H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
    H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
    vector<hsize_t> mdims = mspace->GetDims();

    vector<int> ids;
    mdata->Read(ids, mspace);

    for(int i = 0; i < dims[0]; i++)
    {
        string compStr = comps[i];

        char type;
        istringstream strm(compStr);

        strm >> type;

        CompositeSharedPtr comp =
            MemoryManager<Composite>::AllocateSharedPtr();

        string::size_type indxBeg = compStr.find_first_of('[') + 1;
        string::size_type indxEnd = compStr.find_last_of(']') - 1;

        string indxStr = compStr.substr(indxBeg, indxEnd - indxBeg + 1);
        vector<unsigned int> seqVector;

        ParseUtils::GenerateSeqVector(indxStr, seqVector);

        switch (type)
        {
            case 'V':
                for (auto &i : seqVector)
                {
                    auto it = m_vertSet.find(i);
                    if (it != m_vertSet.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'S':
            case 'E':
                for (auto &i : seqVector)
                {
                    auto it = m_segGeoms.find(i);
                    if (it != m_segGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'Q':
                for (auto &i : seqVector)
                {
                    auto it = m_quadGeoms.find(i);
                    if (it != m_quadGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'T':
                for (auto &i : seqVector)
                {
                    auto it = m_triGeoms.find(i);
                    if (it != m_triGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'F':
                for(auto &i : seqVector)
                {
                    auto it1 = m_quadGeoms.find(i);
                    if (it1 != m_quadGeoms.end())
                    {
                        comp->m_geomVec.push_back(it1->second);
                        continue;
                    }
                    auto it2 = m_triGeoms.find(i);
                    if (it2 != m_triGeoms.end())
                    {
                        comp->m_geomVec.push_back(it2->second);
                    }
                }
                break;
            case 'A':
                for(auto &i : seqVector)
                {
                    auto it = m_tetGeoms.find(i);
                    if (it != m_tetGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'P':
                for(auto &i : seqVector)
                {
                    auto it = m_pyrGeoms.find(i);
                    if (it != m_pyrGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'R':
                for(auto &i : seqVector)
                {
                    auto it = m_prismGeoms.find(i);
                    if (it != m_prismGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
            case 'H':
                for(auto &i : seqVector)
                {
                    auto it = m_hexGeoms.find(i);
                    if (it != m_hexGeoms.end())
                    {
                        comp->m_geomVec.push_back(it->second);
                    }
                }
                break;
        }

        if (comp->m_geomVec.size() > 0)
        {
            m_meshComposites[ids[i]] = comp;
        }
    }
}

CompositeDescriptor MeshGraphHDF5::CreateCompositeDescriptor(
    std::unordered_map<int, int> &id2row)
{
    CompositeDescriptor ret;

    string nm = "COMPOSITE";

    H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
    H5::DataSpaceSharedPtr space = data->GetSpace();
    vector<hsize_t> dims = space->GetDims();

    vector<string> comps;
    data->ReadVectorString(comps, space);

    H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
    H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
    vector<hsize_t> mdims = mspace->GetDims();

    vector<int> ids;
    mdata->Read(ids, mspace);

    for (int i = 0; i < dims[0]; i++)
    {
        string compStr = comps[i];

        char type;
        istringstream strm(compStr);

        strm >> type;

        string::size_type indxBeg = compStr.find_first_of('[') + 1;
        string::size_type indxEnd = compStr.find_last_of(']') - 1;

        string indxStr = compStr.substr(indxBeg, indxEnd - indxBeg + 1);
        vector<unsigned int> seqVector;
        ParseUtils::GenerateSeqVector(indxStr, seqVector);

        LibUtilities::ShapeType shapeType;

        switch (type)
        {
            case 'V':
                shapeType = LibUtilities::ePoint;
                break;
            case 'S':
            case 'E':
                shapeType = LibUtilities::eSegment;
                break;
            case 'Q':
            case 'F':
                // Note that for HDF5, the composite descriptor is only used for
                // partitioning purposes so 'F' tag is not really going to be
                // critical in this context.
                shapeType = LibUtilities::eQuadrilateral;
                break;
            case 'T':
                shapeType = LibUtilities::eTriangle;
                break;
            case 'A':
                shapeType = LibUtilities::eTetrahedron;
                break;
            case 'P':
                shapeType = LibUtilities::ePyramid;
                break;
            case 'R':
                shapeType = LibUtilities::ePrism;
                break;
            case 'H':
                shapeType = LibUtilities::eHexahedron;
                break;
        }

        std::vector<int> filteredVector;
        for (auto &compElmt : seqVector)
        {
            if (id2row.find(compElmt) == id2row.end())
            {
                continue;
            }

            filteredVector.push_back(compElmt);
        }

        if (filteredVector.size() == 0)
        {
            continue;
        }

        ret[ids[i]] = std::make_pair(shapeType, filteredVector);
    }

    return ret;
}


template<class T, typename std::enable_if<T::kDim == 0, int>::type = 0>
inline NekDouble GetGeomData(std::shared_ptr<T> &geom, int i)
{
    return (*geom)(i);
}

template<class T, typename std::enable_if<T::kDim == 1, int>::type = 0>
inline int GetGeomData(std::shared_ptr<T> &geom, int i)
{
    return geom->GetVid(i);
}

template<class T, typename std::enable_if<T::kDim == 2, int>::type = 0>
inline int GetGeomData(std::shared_ptr<T> &geom, int i)
{
    return geom->GetEid(i);
}

template<class T, typename std::enable_if<T::kDim == 3, int>::type = 0>
inline int GetGeomData(std::shared_ptr<T> &geom, int i)
{
    return geom->GetFid(i);
}

template<class T>
void MeshGraphHDF5::WriteGeometryMap(std::map<int, std::shared_ptr<T>> &geomMap,
                                     std::string datasetName)
{
    typedef typename std::conditional<
        std::is_same<T, PointGeom>::value, NekDouble, int>::type DataType;

    const int nGeomData = GetGeomDataDim(geomMap);
    const size_t nGeom = geomMap.size();

    if (nGeom == 0)
    {
        return;
    }

    // Construct a map storing IDs
    vector<int> idMap(nGeom);
    vector<DataType> data(nGeom * nGeomData);

    int cnt1 = 0, cnt2 = 0;
    for (auto &it : geomMap)
    {
        idMap[cnt1++] = it.first;

        for (int j = 0; j < nGeomData; ++j)
        {
            data[cnt2 + j] = GetGeomData(it.second, j);
        }

        cnt2 += nGeomData;
    }

    vector<hsize_t> dims = { static_cast<hsize_t>(nGeom),
                             static_cast<hsize_t>(nGeomData) };
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(data[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet(datasetName, tp, ds);
    dst->Write(data, ds);

    tp = H5::DataType::OfObject(idMap[0]);
    dims = { nGeom };
    ds = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet(datasetName, tp, ds);
    dst->Write(idMap, ds);
}

void MeshGraphHDF5::WriteCurveMap(CurveMap &curves,
                                  std::string dsName,
                                  MeshCurvedPts &curvedPts,
                                  int &ptOffset,
                                  int &newIdx)
{
    vector<int> data, map;

    // Compile curve data.
    for (auto &c : curves)
    {
        map.push_back(c.first);
        data.push_back(c.second->m_points.size());
        data.push_back(c.second->m_ptype);
        data.push_back(ptOffset);

        ptOffset += c.second->m_points.size();

        for (auto &pt : c.second->m_points)
        {
            MeshVertex v;
            v.id = newIdx;
            pt->GetCoords(v.x, v.y, v.z);
            curvedPts.pts.push_back(v);
            curvedPts.index.push_back(newIdx++);
        }
    }

    // Write data.
    vector<hsize_t> dims = { data.size() / 3, 3 };
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(data[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet(dsName, tp, ds);
    dst->Write(data, ds);

    tp = H5::DataType::OfObject(map[0]);
    dims = { map.size() };
    ds = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet(dsName, tp, ds);
    dst->Write(map, ds);
}

void MeshGraphHDF5::WriteCurvePoints(MeshCurvedPts &curvedPts)
{
    vector<double> vertData(curvedPts.pts.size() * 3);

    int cnt = 0;
    for (auto &pt : curvedPts.pts)
    {
        vertData[cnt++] = pt.x;
        vertData[cnt++] = pt.y;
        vertData[cnt++] = pt.z;
    }

    vector<hsize_t> dims = { curvedPts.pts.size(), 3 };
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(vertData[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("CURVE_NODES", tp, ds);
    dst->Write(vertData, ds);
}

void MeshGraphHDF5::WriteComposites(CompositeMap &composites)
{
    vector<string> comps;

    //dont need location map only a id map
    //will filter the composites per parition on read, its easier
    //composites do not need to be written in paralell.
    vector<int> c_map;

    for (auto &cIt : composites)
    {
        if (cIt.second->m_geomVec.size() == 0)
        {
            continue;
        }

        comps.push_back(GetCompositeString(cIt.second));
        c_map.push_back(cIt.first);
    }

    H5::DataTypeSharedPtr tp  = H5::DataType::String();
    H5::DataSpaceSharedPtr ds = H5::DataSpace::OneD(comps.size());
    H5::DataSetSharedPtr dst  = m_mesh->CreateDataSet("COMPOSITE", tp, ds);
    dst->WriteVectorString(comps, ds, tp);

    tp  = H5::DataType::OfObject(c_map[0]);
    ds  = H5::DataSpace::OneD(c_map.size());
    dst = m_maps->CreateDataSet("COMPOSITE", tp, ds);
    dst->Write(c_map, ds);
}

void MeshGraphHDF5::WriteDomain(vector<CompositeMap> &domain)
{
    vector<unsigned int> idxList;
    for (auto cIt = domain[0].begin(); cIt != domain[0].end(); ++cIt)
    {
        idxList.push_back(cIt->first);
    }
    stringstream domString;
    vector<string> doms;
    doms.push_back(ParseUtils::GenerateSeqString(idxList));

    H5::DataTypeSharedPtr tp  = H5::DataType::String();
    H5::DataSpaceSharedPtr ds = H5::DataSpace::OneD(doms.size());
    H5::DataSetSharedPtr dst  = m_mesh->CreateDataSet("DOMAIN", tp, ds);
    dst->WriteVectorString(doms, ds, tp);
}

void MeshGraphHDF5::WriteGeometry(
    std::string                          &outfilename,
    bool                                  defaultExp,
    const LibUtilities::FieldMetaDataMap &metadata)
{
    boost::ignore_unused(metadata);

    vector<string> tmp;
    boost::split(tmp, outfilename, boost::is_any_of("."));
    string filenameXml  = tmp[0] + ".xml";
    string filenameHdf5 = tmp[0] + ".nekg";

    //////////////////
    // XML part
    //////////////////

    //Check to see if a xml of the same name exists
    //if might have boundary conditions etc, we will just alter the geometry
    //tag if needed
    TiXmlDocument *doc = new TiXmlDocument;
    TiXmlElement *root;
    TiXmlElement *geomTag;

    if(boost::filesystem::exists(filenameXml.c_str()))
    {
        ifstream file(filenameXml.c_str());
        file >> (*doc);
        TiXmlHandle docHandle(doc);
        root = docHandle.FirstChildElement("NEKTAR").Element();
        ASSERTL0(root, "Unable to find NEKTAR tag in file.");
        geomTag = root->FirstChildElement("GEOMETRY");
        defaultExp = false;
    }
    else
    {
        TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
        doc->LinkEndChild(decl);
        root = new TiXmlElement("NEKTAR");
        doc->LinkEndChild(root);

        geomTag = new TiXmlElement("GEOMETRY");
        root->LinkEndChild(geomTag);
    }

    // Update attributes with dimensions.
    geomTag->SetAttribute("DIM", m_meshDimension);
    geomTag->SetAttribute("SPACE", m_spaceDimension);
    geomTag->SetAttribute("HDF5FILE", filenameHdf5);

    geomTag->Clear();

    if (defaultExp)
    {
        TiXmlElement *expTag = new TiXmlElement("EXPANSIONS");

        for (auto it = m_meshComposites.begin(); it != m_meshComposites.end();
             it++)
        {
            if (it->second->m_geomVec[0]->GetShapeDim() == m_meshDimension)
            {
                TiXmlElement *exp = new TiXmlElement("E");
                exp->SetAttribute(
                    "COMPOSITE",
                    "C[" + boost::lexical_cast<string>(it->first) + "]");
                exp->SetAttribute("NUMMODES", 4);
                exp->SetAttribute("TYPE", "MODIFIED");
                exp->SetAttribute("FIELDS", "u");

                expTag->LinkEndChild(exp);
            }
        }
        root->LinkEndChild(expTag);
    }

    doc->SaveFile(filenameXml);

    //////////////////
    // HDF5 part
    //////////////////

    // This is serial IO so we will just override any existing file.
    m_file = H5::File::Create(filenameHdf5, H5F_ACC_TRUNC);
    auto hdfRoot = m_file->CreateGroup("NEKTAR");
    auto hdfRoot2 = hdfRoot->CreateGroup("GEOMETRY");

    // Write format version.
    hdfRoot2->SetAttribute("FORMAT_VERSION", FORMAT_VERSION);

    // Create main groups.
    m_mesh = hdfRoot2->CreateGroup("MESH");
    m_maps = hdfRoot2->CreateGroup("MAPS");

    WriteGeometryMap(m_vertSet, "VERT");
    WriteGeometryMap(m_segGeoms, "SEG");
    if (m_meshDimension > 1)
    {
        WriteGeometryMap(m_triGeoms, "TRI");
        WriteGeometryMap(m_quadGeoms, "QUAD");
    }
    if (m_meshDimension > 2)
    {
        WriteGeometryMap(m_tetGeoms, "TET");
        WriteGeometryMap(m_pyrGeoms, "PYR");
        WriteGeometryMap(m_prismGeoms, "PRISM");
        WriteGeometryMap(m_hexGeoms, "HEX");
    }

    // Write curves
    int ptOffset = 0, newIdx = 0;
    MeshCurvedPts curvePts;
    WriteCurveMap(m_curvedEdges, "CURVE_EDGE", curvePts, ptOffset, newIdx);
    WriteCurveMap(m_curvedFaces, "CURVE_FACE", curvePts, ptOffset, newIdx);
    WriteCurvePoints(curvePts);

    // Write composites and domain.
    WriteComposites(m_meshComposites);
    WriteDomain(m_domain);
}

}
}
