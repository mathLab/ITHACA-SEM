////////////////////////////////////////////////////////////////////////////////
//
//  File: BoundaryConditions.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//
////////////////////////////////////////////////////////////////////////////////

#include "MeshGraphHDF5.h"

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SpatialDomains/MeshEntities.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <type_traits>


#include <tinyxml.h>
using namespace std;
using namespace Nektar::LibUtilities;

namespace Nektar
{
namespace SpatialDomains
{

std::string MeshGraphHDF5::className =
    GetMeshGraphFactory().RegisterCreatorFunction(
        "HDF5", MeshGraphHDF5::create, "IO with HDF5 geometry");

void MeshGraphHDF5::ReadGeometry(
    DomainRangeShPtr rng,
    bool fillGraph)
{
    int err;

    //we use the xml geom to find information about the HDF5 file
    m_xmlGeom            = m_session->GetElement("NEKTAR/GEOMETRY");
    TiXmlAttribute *attr = m_xmlGeom->FirstAttribute();
    m_meshPartitioned    = false;
    m_meshDimension      = 3;
    m_spaceDimension     = 3;
    while (attr)
    {
        std::string attrName(attr->Name());
        if (attrName == "DIM")
        {
            err = attr->QueryIntValue(&m_meshDimension);
            ASSERTL1(err == TIXML_SUCCESS, "Unable to read mesh dimension.");
        }
        else if (attrName == "SPACE")
        {
            err = attr->QueryIntValue(&m_spaceDimension);
            ASSERTL1(err == TIXML_SUCCESS, "Unable to read space dimension.");
        }
        else if (attrName == "PARTITION")
        {
            ASSERTL0(false,"PARTITION parameter should only be use in xml meshes");
        }
        else if(attrName == "HDF5FILE")
        {
            m_hdf5Name = attr->Value();
            ASSERTL1(err == TIXML_SUCCESS, "Unable to read hdf5 name.");
        }
        else if(attrName == "PARTITIONED")
        {
            m_meshPartitioned = true;
            if(m_session->GetComm()->GetSize() == 1)
            {
                //mesh may be paritioned but we want the whole mesh
                m_meshPartitioned = false;
            }
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

    ASSERTL0(m_hdf5Name.size() > 0, "unable to obtain mesh file name");

    ASSERTL1(m_meshDimension <= m_spaceDimension,
             "Mesh dimension greater than space dimension");

    //load the HDF5 mesh
    m_file = H5::File::Open(m_hdf5Name, H5F_ACC_RDWR);
    m_mesh = m_file->OpenGroup("mesh");
    m_maps = m_file->OpenGroup("maps");

    ReadGeometryMap(m_vertSet, "vert");
    //ReadCurves();
    if (m_meshDimension >= 2)
    {
        ReadGeometryMap(m_segGeoms, "seg", m_curvedEdges);
        if (m_meshDimension == 3)
        {
            ReadFaces();
        }
    }
    ReadElements();
    ReadComposites();
    ReadDomain();
    ReadExpansions();
}

void MeshGraphHDF5::PartitionMesh(LibUtilities::SessionReaderSharedPtr session)
{
    m_session = session;
    // Don't do anything yet!
}

template<class T, typename std::enable_if<T::kDim == 0, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap, int s)
{
    return s;
}

template<class T, typename std::enable_if<T::kDim == 1, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap, int s)
{
    return T::kNverts;
}

template<class T, typename std::enable_if<T::kDim == 2, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap, int s)
{
    return T::kNedges;
}

template<class T, typename std::enable_if<T::kDim == 3, int>::type = 0>
inline int GetGeomDataDim(std::map<int, std::shared_ptr<T>> &geomMap, int s)
{
    return T::kNfaces;
}

template<class T, typename DataType> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<T>> &geomMap, int id,
    DataType *data, CurveSharedPtr curve)
{
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PointGeom>> &geomMap, int id,
    NekDouble *data, CurveSharedPtr curve)
{
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
    TriGeomSharedPtr faces[4] = {
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[0])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[1])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[2])),
        std::static_pointer_cast<TriGeom>(GetGeometry2D(data[3]))
    };
    geomMap[id] = MemoryManager<TetGeom>::AllocateSharedPtr(id, faces);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PyrGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    Geometry2DSharedPtr faces[5] = {
        GetGeometry2D(data[0]), GetGeometry2D(data[1]), GetGeometry2D(data[2]),
        GetGeometry2D(data[3]), GetGeometry2D(data[4])
    };
    geomMap[id] = MemoryManager<PyrGeom>::AllocateSharedPtr(id, faces);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<PrismGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    Geometry2DSharedPtr faces[5] = {
        GetGeometry2D(data[0]), GetGeometry2D(data[1]), GetGeometry2D(data[2]),
        GetGeometry2D(data[3]), GetGeometry2D(data[4])
    };
    geomMap[id] = MemoryManager<PrismGeom>::AllocateSharedPtr(id, faces);
}

template<> void MeshGraphHDF5::ConstructGeomObject(
    std::map<int, std::shared_ptr<HexGeom>> &geomMap, int id, int *data,
    CurveSharedPtr curve)
{
    QuadGeomSharedPtr faces[6] = {
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[0])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[1])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[2])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[3])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[4])),
        std::static_pointer_cast<QuadGeom>(GetGeometry2D(data[5]))
    };
    geomMap[id] = MemoryManager<HexGeom>::AllocateSharedPtr(id, faces);
}

template<class T>
void MeshGraphHDF5::ReadGeometryMap(
    std::map<int, std::shared_ptr<T>> &geomMap,
    std::string dataSet,
    const CurveMap &curveMap,
    const std::unordered_set<unsigned int> &readIds)
{
    typedef typename std::conditional<
        std::is_same<T, PointGeom>::value, NekDouble, int>::type DataType;

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

    // Read IDs: TODO: could be done in chunks
    vector<int> ids, newids;
    mdata->Read(ids, mspace);

    vector<DataType> geomData;

    const int nGeomData = GetGeomDataDim(geomMap, m_spaceDimension);

    if (readIds.size() > 0)
    {
        // Selective reading
        int i = 0;
        for (auto &id : ids)
        {
            if (readIds.find(id) != readIds.end())
            {
                space->AppendRange(i, nGeomData);
                newids.push_back(id);
            }
            i += nGeomData;
        }
    }
    else
    {
        newids = ids;
    }

    // Read data (collectively)
    data->Read(geomData, space);

    const int nRows = geomData.size() / nGeomData;

    CurveSharedPtr empty;

    // Construct geometry object.
    if (curveMap.size() > 0)
    {
        for(int i = 0, cnt = 0; i < nRows; i++, cnt += nGeomData)
        {
            auto cIt = curveMap.find(newids[i]);
            ConstructGeomObject(
                geomMap, newids[i], &geomData[cnt],
                cIt == curveMap.end() ? cIt->second : empty);
        }
    }
    else
    {
        for(int i = 0, cnt = 0; i < nRows; i++, cnt += nGeomData)
        {
            ConstructGeomObject(geomMap, newids[i], &geomData[cnt], empty);
        }
    }
}

void MeshGraphHDF5::ReadCurves()
{
    string nmEdge = "curve_edge";
    string nmFace = "curve_face";
    string nmNodes = "curve_nodes";

    /*ta[i * 4 + 0] = edgeInfo[i].entityid;
    data[i * 4 + 1] = edgeInfo[i].npoints;
    data[i * 4 + 2] = edgeInfo[i].ptype;
    data[i * 4 + 3] = edgeInfo[i].ptoffset;  */

    H5::DataSetSharedPtr dataEdge = m_mesh->OpenDataSet(nmEdge);
    H5::DataSpaceSharedPtr spaceEdge = dataEdge->GetSpace();
    vector<hsize_t> dimsEdge = spaceEdge->GetDims();

    vector<int> rawDataEdge;

    dataEdge->Read(rawDataEdge, spaceEdge);

    H5::DataSetSharedPtr dataFace = m_mesh->OpenDataSet(nmEdge);
    H5::DataSpaceSharedPtr spaceFace = dataFace->GetSpace();
    vector<hsize_t> dimsFace = spaceFace->GetDims();

    vector<int> rawDataFace;

    dataEdge->Read(rawDataFace, spaceFace);

    H5::DataSetSharedPtr dataNodes = m_mesh->OpenDataSet(nmNodes);
    H5::DataSpaceSharedPtr spaceNodes = dataNodes->GetSpace();

    for(int i = 0; i < dimsEdge[0]; i++)
    {
        int numPoints = rawDataEdge[i*4+1];
        int ptoffset = rawDataEdge[i*4+3];
        vector<vector<int> > required;
        for(int j = 0; j < numPoints; j++)
        {
            vector<int> x, y, z;
            x.push_back(ptoffset+j);
            x.push_back(0);
            y.push_back(ptoffset+j);
            y.push_back(1);
            z.push_back(ptoffset+j);
            z.push_back(2);
            required.push_back(x);
            required.push_back(y);
            required.push_back(z);
        }

        int edgeId = rawDataEdge[i*4+0];
        LibUtilities::PointsType type;
        type = (LibUtilities::PointsType)rawDataEdge[i*4+2];

        CurveSharedPtr curve =
            MemoryManager<Curve>::AllocateSharedPtr(edgeId, type);

        vector<NekDouble> points;
        dataNodes->Read(points, spaceNodes, required);

        for(int j = 0; j < numPoints; j++)
        {
            PointGeomSharedPtr vert = MemoryManager<PointGeom>::AllocateSharedPtr
                (m_meshDimension, i, points[j*3+0], points[j*3+1], points[j*3+2]);

            curve->m_points.push_back(vert);
        }

        m_curvedEdges[edgeId] = curve;
    }

    for(int i = 0; i < dimsFace[0]; i++)
    {
        int numPoints = rawDataFace[i*4+1];
        int ptoffset = rawDataFace[i*4+3];
        vector<vector<int> > required;
        for(int j = 0; j < numPoints; j++)
        {
            vector<int> x, y, z;
            x.push_back(ptoffset+j);
            x.push_back(0);
            y.push_back(ptoffset+j);
            y.push_back(1);
            z.push_back(ptoffset+j);
            z.push_back(2);
            required.push_back(x);
            required.push_back(y);
            required.push_back(z);
        }

        int faceId = rawDataFace[i*4+0];
        LibUtilities::PointsType type;
        type = (LibUtilities::PointsType)rawDataFace[i*4+2];

        CurveSharedPtr curve =
            MemoryManager<Curve>::AllocateSharedPtr(faceId, type);

        vector<NekDouble> points;
        dataNodes->Read(points, spaceNodes, required);

        for(int j = 0; j < numPoints; j++)
        {
            PointGeomSharedPtr vert = MemoryManager<PointGeom>::AllocateSharedPtr
                (m_meshDimension, i, points[j*3+0], points[j*3+1], points[j*3+2]);

            curve->m_points.push_back(vert);
        }

        m_curvedFaces[faceId] = curve;
    }
}

void MeshGraphHDF5::ReadDomain()
{
    map<int, CompositeSharedPtr> fullDomain;
    GetCompositeList("0", fullDomain);
    m_domain.push_back(fullDomain);
}

void MeshGraphHDF5::ReadFaces()
{
    ReadGeometryMap(m_triGeoms, "tri", m_curvedFaces);
    ReadGeometryMap(m_quadGeoms, "quad", m_curvedFaces);
}

void MeshGraphHDF5::ReadElements()
{
    if(m_meshDimension == 1)
    {
        ReadGeometryMap(m_segGeoms, "seg", m_curvedEdges);
        return;
    }
    else if (m_meshDimension == 2)
    {
        ReadFaces();
        return;
    }

    ReadGeometryMap(m_tetGeoms, "tet");
    ReadGeometryMap(m_pyrGeoms, "pyr");
    ReadGeometryMap(m_prismGeoms, "prism");
    ReadGeometryMap(m_hexGeoms, "hex");
}

void MeshGraphHDF5::ReadComposites()
{
    string nm = "composite";

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
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_vertSet[i]);
                }
                break;
            case 'S':
            case 'E':
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_segGeoms[i]);
                }
                break;
            case 'T':
            case 'Q':
            case 'F':
                for(auto & i : seqVector)
                {
                    comp->m_geomVec.push_back(GetGeometry2D(i));
                }
                break;
            case 'A':
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_tetGeoms[i]);
                }
                break;
            case 'P':
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_pyrGeoms[i]);
                }
                break;
            case 'R':
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_prismGeoms[i]);
                }
                break;
            case 'H':
                for(auto &i : seqVector)
                {
                    comp->m_geomVec.push_back(m_hexGeoms[i]);
                }
                break;
        }

        m_meshComposites[ids[i]] = comp;

    }
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

    const int nGeomData = GetGeomDataDim(geomMap, m_spaceDimension);
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

void MeshGraphHDF5::WriteCurves(CurveMap &edges, CurveMap &faces)
{
    if(edges.size() == 0 && faces.size() == 0)
    {
        return;
    }

    vector<MeshCurvedInfo> edgeInfo;
    vector<MeshCurvedInfo> faceInfo;
    MeshCurvedPts curvedPts;
    curvedPts.id = 0;
    int ptOffset = 0;
    int newIdx = 0;

    for (auto& i : edges)
    {
        MeshCurvedInfo cinfo;
        cinfo.entityid = i.first;
        cinfo.npoints = i.second->m_points.size();
        cinfo.ptype = i.second->m_ptype;
        cinfo.ptid = 0;
        cinfo.ptoffset = ptOffset;

        edgeInfo.push_back(cinfo);

        for(int j = 0; j < i.second->m_points.size(); j++)
        {
            MeshVertex v;
            v.id = newIdx;
            i.second->m_points[j]->GetCoords(v.x,v.y,v.z);
            curvedPts.pts.push_back(v);
            curvedPts.index.push_back(newIdx);
            newIdx++;
        }
        ptOffset += cinfo.npoints;
    }

    for (auto& i : faces)
    {
        MeshCurvedInfo cinfo;
        cinfo.entityid = i.first;
        cinfo.npoints = i.second->m_points.size();
        cinfo.ptype = i.second->m_ptype;
        cinfo.ptid = 0;
        cinfo.ptoffset = ptOffset;

        faceInfo.push_back(cinfo);

        for(int j = 0; j < i.second->m_points.size(); j++)
        {
            MeshVertex v;
            v.id = newIdx;
            i.second->m_points[j]->GetCoords(v.x,v.y,v.z);
            curvedPts.pts.push_back(v);
            curvedPts.index.push_back(newIdx);
            newIdx++;
        }
        ptOffset += cinfo.npoints;
    }

    if(edgeInfo.size())
    {
        vector<int> data(edgeInfo.size() * 3);
        vector<int> map(edgeInfo.size() * 2);
        for (int i = 0; i < edgeInfo.size(); i++)
        {
            map[i * 2 + 0] = edgeInfo[i].entityid;
            map[i * 2 + 1] = i;
            data[i * 3 + 0] = edgeInfo[i].npoints;
            data[i * 3 + 1] = edgeInfo[i].ptype;
            data[i * 3 + 2] = edgeInfo[i].ptoffset;
        }
        vector<hsize_t> dims;
        dims.push_back(edgeInfo.size());
        dims.push_back(3);
        H5::DataTypeSharedPtr tp = H5::DataType::OfObject(data[0]);
        H5::DataSpaceSharedPtr ds =
            std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
        H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("curve_edge", tp, ds);
        dst->Write(data, ds);

        tp = H5::DataType::OfObject(map[0]);
        dims.clear();
        dims.push_back(edgeInfo.size());
        dims.push_back(2);
        ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
        dst = m_maps->CreateDataSet("curve_edge", tp, ds);
        dst->Write(map, ds);
    }

    if(faceInfo.size())
    {
        vector<int> data(faceInfo.size() * 3);
        vector<int> map(faceInfo.size() * 2);
        for (int i = 0; i < faceInfo.size(); i++)
        {
            map[i * 2 + 0] = faceInfo[i].entityid;
            map[i * 2 + 1] = i;
            data[i * 3 + 0] = faceInfo[i].npoints;
            data[i * 3 + 1] = faceInfo[i].ptype;
            data[i * 3 + 2] = faceInfo[i].ptoffset;
        }
        vector<hsize_t> dims;
        dims.push_back(faceInfo.size());
        dims.push_back(3);
        H5::DataTypeSharedPtr tp = H5::DataType::OfObject(data[0]);
        H5::DataSpaceSharedPtr ds =
            std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
        H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("curve_face", tp, ds);
        dst->Write(data, ds);

        tp = H5::DataType::OfObject(map[0]);
        dims.clear();
        dims.push_back(faceInfo.size());
        dims.push_back(2);
        ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
        dst = m_maps->CreateDataSet("curve_face", tp, ds);
        dst->Write(map, ds);
    }

    if(edgeInfo.size() || faceInfo.size())
    {
        vector<double> vertData(curvedPts.pts.size() * m_spaceDimension);
        for (int i = 0; i < curvedPts.pts.size(); i++)
        {
            for (int j = 0; j < m_spaceDimension; j++)
            {
                vertData[i * 3 + 0] = curvedPts.pts[i].x;
                if(m_meshDimension > 1)
                {
                    vertData[i * 3 + 1] = curvedPts.pts[i].y;
                }
                if(m_meshDimension > 2)
                {
                    vertData[i * 3 + 2] = curvedPts.pts[i].z;
                }
            }
        }
        vector<hsize_t> dims;
        dims.push_back(curvedPts.pts.size());
        dims.push_back(m_spaceDimension);

        H5::DataTypeSharedPtr tp = H5::DataType::OfObject(vertData[0]);
        H5::DataSpaceSharedPtr ds =
            std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
        H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("curve_nodes", tp, ds);
        dst->Write(vertData, ds);
    }
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
    H5::DataSetSharedPtr dst  = m_mesh->CreateDataSet("composite", tp, ds);
    dst->WriteVectorString(comps, tp);

    tp  = H5::DataType::OfObject(c_map[0]);
    ds  = H5::DataSpace::OneD(c_map.size());
    dst = m_maps->CreateDataSet("composite", tp, ds);
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
    H5::DataSetSharedPtr dst  = m_mesh->CreateDataSet("domain", tp, ds);
    dst->WriteVectorString(doms, tp);
}

void MeshGraphHDF5::WriteGeometry(
    std::string                          &outfilename,
    bool                                  defaultExp,
    const LibUtilities::FieldMetaDataMap &metadata)
{
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

    //this is serial IO so we will just override any exisiting file
    m_file = H5::File::Create(filenameHdf5, H5F_ACC_TRUNC);
    m_mesh = m_file->CreateGroup("mesh");
    m_maps = m_file->CreateGroup("maps");

    WriteGeometryMap(m_vertSet, "vert");
    WriteGeometryMap(m_segGeoms, "seg");
    if (m_meshDimension > 1)
    {
        WriteGeometryMap(m_triGeoms, "tri");
        WriteGeometryMap(m_quadGeoms, "quad");
    }
    if (m_meshDimension > 2)
    {
        WriteGeometryMap(m_tetGeoms, "tet");
        WriteGeometryMap(m_pyrGeoms, "pyr");
        WriteGeometryMap(m_prismGeoms, "prism");
        WriteGeometryMap(m_hexGeoms, "hex");
    }
    WriteCurves(m_curvedEdges, m_curvedFaces);
    WriteComposites(m_meshComposites);
    WriteDomain(m_domain);
}

void MeshGraphHDF5::WriteGeometry(std::string outname,
                                  std::vector<std::set<unsigned int>> elements,
                                  std::vector<unsigned int> partitions)
{
    //if we hit this function we need to be aware that the mesh
    //already exisits in the file so it only needs ammending
    //but it may have partition information already so need to delte it
    //first
}

}
}
