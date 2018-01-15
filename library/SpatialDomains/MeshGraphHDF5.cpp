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

    ReadVertices();
    ReadCurves();
    if (m_meshDimension >= 2)
    {
        ReadEdges();
        if (m_meshDimension == 3)
        {
            ReadFaces();
        }
    }
    ReadElements();
    ReadComposites();
    ReadDomain();
}

void MeshGraphHDF5::PartitionMesh(LibUtilities::SessionReaderSharedPtr session)
{
    // Don't do anything yet!
}

void MeshGraphHDF5::ReadVertices()
{
    string nm = "vert";

    H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
    H5::DataSpaceSharedPtr space = data->GetSpace();
    vector<hsize_t> dims = space->GetDims();

    vector<NekDouble> verts;
    data->Read(verts, space);

    H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
    H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
    vector<hsize_t> mdims = mspace->GetDims();

    vector<int> ids;
    mdata->Read(ids, mspace);

    ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

    for(int i = 0; i < dims[0]; i++)
    {
        PointGeomSharedPtr vert(
            MemoryManager<PointGeom>::AllocateSharedPtr(
                                m_spaceDimension, ids[i],
                                    verts[i*dims[1] + 0],
                                    verts[i*dims[1] + 1],
                                    verts[i*dims[1] + 2]));
        m_vertSet[ids[i]] = vert;
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

}

void MeshGraphHDF5::ReadEdges()
{
    string nm = "seg";

    H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
    H5::DataSpaceSharedPtr space = data->GetSpace();
    vector<hsize_t> dims = space->GetDims();

    vector<int> edges;
    data->Read(edges, space);

    H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
    H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
    vector<hsize_t> mdims = mspace->GetDims();

    vector<int> ids;
    mdata->Read(ids, mspace);

    ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

    for(int i = 0; i < dims[0]; i++)
    {
        PointGeomSharedPtr verts[2] = {GetVertex(edges[i*2+0]),
                                       GetVertex(edges[i*2+1])};

        SegGeomSharedPtr edge;

        auto it = m_curvedEdges.find(ids[i]);

        if(it == m_curvedEdges.end())
        {
            edge = MemoryManager<SegGeom>::AllocateSharedPtr(
                ids[i], m_spaceDimension, verts);
        }
        else
        {
            edge = MemoryManager<SegGeom>::AllocateSharedPtr(
                ids[i], m_spaceDimension, verts, it->second);
        }

        m_segGeoms[ids[i]] = edge;
    }
}

void MeshGraphHDF5::ReadFaces()
{
    //tris
    {
        string nm = "tri";

        if(m_mesh->ContainsDataSet(nm))
        {
            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> tris;
            data->Read(tris, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                SegGeomSharedPtr edges[TriGeom::kNedges] = {
                    GetSegGeom(tris[i*3 + 0]),
                    GetSegGeom(tris[i*3 + 1]),
                    GetSegGeom(tris[i*3 + 2])};

                TriGeomSharedPtr tri;

                auto it = m_curvedFaces.find(ids[i]);

                if(it == m_curvedFaces.end())
                {
                    tri = MemoryManager<TriGeom>::AllocateSharedPtr(
                        ids[i], edges);
                }
                else
                {
                    tri = MemoryManager<TriGeom>::AllocateSharedPtr(
                        ids[i], edges, it->second);
                }

                m_triGeoms[ids[i]] = tri;
            }
        }
    }

    //quads
    {
        string nm = "quad";

        if(m_mesh->ContainsDataSet(nm))
        {
            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> quads;
            data->Read(quads, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                SegGeomSharedPtr edges[QuadGeom::kNedges] = {
                    GetSegGeom(quads[i*4 + 0]),
                    GetSegGeom(quads[i*4 + 1]),
                    GetSegGeom(quads[i*4 + 2]),
                    GetSegGeom(quads[i*4 + 3])};

                QuadGeomSharedPtr quad;

                auto it = m_curvedFaces.find(ids[i]);

                if(it == m_curvedFaces.end())
                {
                    quad = MemoryManager<QuadGeom>::AllocateSharedPtr(
                        ids[i], edges);
                }
                else
                {
                    quad = MemoryManager<QuadGeom>::AllocateSharedPtr(
                        ids[i], edges, it->second);
                }

                m_quadGeoms[ids[i]] = quad;
            }
        }
    }
}

void MeshGraphHDF5::ReadElements()
{
    if(m_meshDimension == 1)
    {
        ReadEdges();
        return;
    }
    else if (m_meshDimension == 2)
    {
        ReadFaces();
        return;
    }

    //Tets
    {
        string nm = "tet";

        if(m_mesh->ContainsDataSet(nm))
        {

            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> tets;
            data->Read(tets, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                TriGeomSharedPtr faces[TetGeom::kNfaces];

                for(int j = 0; j < TetGeom::kNfaces; j++)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(tets[i*4+j]);
                    faces[j] = std::static_pointer_cast<TriGeom>(face);
                }

                TetGeomSharedPtr tet = MemoryManager<TetGeom>::AllocateSharedPtr(
                    ids[i], faces);

                m_tetGeoms[ids[i]] = tet;
                PopulateFaceToElMap(tet, TetGeom::kNfaces);
            }
        }
    }
    //Pyrs
    {
        string nm = "pyr";

        if(m_mesh->ContainsDataSet(nm))
        {

            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> pyrs;
            data->Read(pyrs, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                Geometry2DSharedPtr faces[PyrGeom::kNfaces];

                for(int j = 0; j < PyrGeom::kNfaces; j++)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(pyrs[i*5+j]);
                    faces[j] = face;
                }

                PyrGeomSharedPtr pyr = MemoryManager<PyrGeom>::AllocateSharedPtr(
                    ids[i], faces);

                m_pyrGeoms[ids[i]] = pyr;
                PopulateFaceToElMap(pyr, PyrGeom::kNfaces);
            }
        }
    }
    //Prism
    {
        string nm = "prism";

        if(m_mesh->ContainsDataSet(nm))
        {
            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> prisms;
            data->Read(prisms, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                Geometry2DSharedPtr faces[PrismGeom::kNfaces];

                for(int j = 0; j < PrismGeom::kNfaces; j++)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(prisms[i*5+j]);
                    faces[j] = face;
                }

                PrismGeomSharedPtr prism = MemoryManager<PrismGeom>::AllocateSharedPtr(
                    ids[i], faces);

                m_prismGeoms[ids[i]] = prism;
                PopulateFaceToElMap(prism, PrismGeom::kNfaces);
            }
        }
    }
    //Hex
    {
        string nm = "hex";

        if(m_mesh->ContainsDataSet(nm))
        {
            H5::DataSetSharedPtr data = m_mesh->OpenDataSet(nm);
            H5::DataSpaceSharedPtr space = data->GetSpace();
            vector<hsize_t> dims = space->GetDims();

            vector<int> hexs;
            data->Read(hexs, space);

            H5::DataSetSharedPtr mdata = m_maps->OpenDataSet(nm);
            H5::DataSpaceSharedPtr mspace = mdata->GetSpace();
            vector<hsize_t> mdims = mspace->GetDims();

            vector<int> ids;
            mdata->Read(ids, mspace);

            ASSERTL0(mdims[0] == dims[0], "map and data set lengths do not match");

            for(int i = 0; i < dims[0]; i++)
            {
                QuadGeomSharedPtr faces[HexGeom::kNfaces];

                for(int j = 0; j < HexGeom::kNfaces; j++)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(hexs[i*6+j]);
                    faces[j] = std::static_pointer_cast<QuadGeom>(face);
                }

                HexGeomSharedPtr hex = MemoryManager<HexGeom>::AllocateSharedPtr(
                    ids[i], faces);

                m_hexGeoms[ids[i]] = hex;
                PopulateFaceToElMap(hex, HexGeom::kNfaces);
            }
        }
    }

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

        ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);

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

void MeshGraphHDF5::WriteVertices(PointGeomMap &verts)
{
    vector<double> vertData(verts.size() * m_spaceDimension);
    //in the maps we write the id of the vertex and its location in
    //the hdf5 list, in serial the location is a 1 to 1 map
    //needed for parralell.
    vector<int> v_map(verts.size() * 2);

    int i = 0;
    for (auto it = verts.begin(); it != verts.end(); it++, i++)
    {
        v_map[i * 2 + 0] = it->first;
        v_map[i * 2 + 1] = i;
        Array<OneD, NekDouble> l(m_spaceDimension);
        it->second->GetCoords(l);
        for (int j = 0; j < m_spaceDimension; j++)
        {
            vertData[i * m_spaceDimension + j] = l[j];
        }
    }
    vector<hsize_t> dims;
    dims.push_back(verts.size());
    dims.push_back(m_spaceDimension);

    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(vertData[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("vert", tp, ds);
    dst->Write(vertData, ds);

    tp = H5::DataType::OfObject(v_map[0]);
    dims.clear();
    dims.push_back(verts.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("vert", tp, ds);
    dst->Write(v_map, ds);
}

void MeshGraphHDF5::WriteEdges(SegGeomMap &edges)
{
    vector<int> edgeData(edges.size() * 2);
    vector<int> e_map(edges.size() * 2);
    int i = 0;
    for (auto it = edges.begin(); it != edges.end(); it++, i++)
    {
        e_map[i * 2 + 0] = it->first;
        e_map[i * 2 + 1] = i;
        edgeData[i * 2 + 0] = it->second->GetVertex(0)->GetGlobalID();
        edgeData[i * 2 + 1] = it->second->GetVertex(1)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(edges.size());
    dims.push_back(2);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(edgeData[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("seg", tp, ds);
    dst->Write(edgeData, ds);

    tp = H5::DataType::OfObject(e_map[0]);
    dims.clear();
    dims.push_back(edges.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("seg", tp, ds);
    dst->Write(e_map, ds);
}

void MeshGraphHDF5::WriteTris(TriGeomMap &tris)
{
    if(tris.size() == 0)
    {
        return;
    }

    vector<int> triData(tris.size() * 3);
    vector<int> t_map(tris.size() * 2);
    int i = 0;
    for (auto it = tris.begin(); it != tris.end(); it++, i++)
    {
        t_map[i * 2 + 0] = it->first;
        t_map[i * 2 + 1] = i;
        triData[i * 3 + 0] = it->second->GetEdge(0)->GetGlobalID();
        triData[i * 3 + 1] = it->second->GetEdge(1)->GetGlobalID();
        triData[i * 3 + 2] = it->second->GetEdge(2)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(tris.size());
    dims.push_back(3);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(triData[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("tri", tp, ds);
    dst->Write(triData, ds);

    tp = H5::DataType::OfObject(t_map[0]);
    dims.clear();
    dims.push_back(tris.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("tri", tp, ds);
    dst->Write(t_map, ds);
}

void MeshGraphHDF5::WriteQuads(QuadGeomMap &quads)
{
    if(quads.size() == 0)
    {
        return;
    }

    vector<int> quadData(quads.size() * 4);
    vector<int> q_map(quads.size() * 2);
    int i = 0;
    for (auto it = quads.begin(); it != quads.end(); it++, i++)
    {
        q_map[i * 2 + 0] = it->first;
        q_map[i * 2 + 1] = i;
        quadData[i * 4 + 0] = it->second->GetEdge(0)->GetGlobalID();
        quadData[i * 4 + 1] = it->second->GetEdge(1)->GetGlobalID();
        quadData[i * 4 + 2] = it->second->GetEdge(2)->GetGlobalID();
        quadData[i * 4 + 3] = it->second->GetEdge(3)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(quads.size());
    dims.push_back(4);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(quadData[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("quad", tp, ds);
    dst->Write(quadData, ds);

    tp = H5::DataType::OfObject(q_map[0]);
    dims.clear();
    dims.push_back(quads.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("quad", tp, ds);
    dst->Write(q_map, ds);
}

void MeshGraphHDF5::WriteTets(TetGeomMap &tets)
{
    if(tets.size() == 0)
    {
        return;
    }

    vector<int> els(tets.size() * 4);
    vector<int> i_map(tets.size() * 2);
    int i = 0;
    for (auto it = tets.begin(); it != tets.end(); it++, i++)
    {
        i_map[i * 2 + 0] = it->first;
        i_map[i * 2 + 1] = i;
        els[i * 4 + 0] = it->second->GetFace(0)->GetGlobalID();
        els[i * 4 + 1] = it->second->GetFace(1)->GetGlobalID();
        els[i * 4 + 2] = it->second->GetFace(2)->GetGlobalID();
        els[i * 4 + 3] = it->second->GetFace(3)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(tets.size());
    dims.push_back(4);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(els[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("tet", tp, ds);
    dst->Write(els, ds);

    tp = H5::DataType::OfObject(i_map[0]);
    dims.clear();
    dims.push_back(tets.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("tet", tp, ds);
    dst->Write(i_map, ds);
}

void MeshGraphHDF5::WritePyrs(PyrGeomMap &pyrs)
{
    if(pyrs.size() == 0)
    {
        return;
    }

    vector<int> els(pyrs.size() * 5);
    vector<int> i_map(pyrs.size() * 2);
    int i = 0;
    for (auto it = pyrs.begin(); it != pyrs.end(); it++, i++)
    {
        i_map[i * 2 + 0] = it->first;
        i_map[i * 2 + 1] = i;
        els[i * 5 + 0] = it->second->GetFace(0)->GetGlobalID();
        els[i * 5 + 1] = it->second->GetFace(1)->GetGlobalID();
        els[i * 5 + 2] = it->second->GetFace(2)->GetGlobalID();
        els[i * 5 + 3] = it->second->GetFace(3)->GetGlobalID();
        els[i * 5 + 4] = it->second->GetFace(4)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(pyrs.size());
    dims.push_back(5);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(els[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("pyr", tp, ds);
    dst->Write(els, ds);

    tp = H5::DataType::OfObject(i_map[0]);
    dims.clear();
    dims.push_back(pyrs.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("pyr", tp, ds);
    dst->Write(i_map, ds);
}

void MeshGraphHDF5::WritePrisms(PrismGeomMap &prisms)
{
    if(prisms.size() == 0)
    {
        return;
    }

    vector<int> els(prisms.size() * 5);
    vector<int> i_map(prisms.size() * 2);
    int i = 0;
    for (auto it = prisms.begin(); it != prisms.end();
         it++, i++)
    {
        i_map[i * 2 + 0] = it->first;
        i_map[i * 2 + 1] = i;
        els[i * 5 + 0] = it->second->GetFace(0)->GetGlobalID();
        els[i * 5 + 1] = it->second->GetFace(1)->GetGlobalID();
        els[i * 5 + 2] = it->second->GetFace(2)->GetGlobalID();
        els[i * 5 + 3] = it->second->GetFace(3)->GetGlobalID();
        els[i * 5 + 4] = it->second->GetFace(4)->GetGlobalID();
    }
    vector<hsize_t> dims;
    dims.push_back(prisms.size());
    dims.push_back(5);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(els[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("prism", tp, ds);
    dst->Write(els, ds);

    tp = H5::DataType::OfObject(i_map[0]);
    dims.clear();
    dims.push_back(prisms.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("prism", tp, ds);
    dst->Write(i_map, ds);
}

void MeshGraphHDF5::WriteHexs(HexGeomMap &hexs)
{
    if(hexs.size() == 0)
    {
        return;
    }

    vector<int> els(hexs.size() * 6);
    vector<int> i_map(hexs.size() * 2);
    int i = 0;
    for (auto it = hexs.begin(); it != hexs.end(); it++, i++)
    {
        i_map[i * 2 + 0] = it->first;
        i_map[i * 2 + 1] = i;
        for(int j = 0; j < it->second->GetNumFaces(); j++)
        {
            els[i * 6 + j] = it->second->GetFace(j)->GetGlobalID();
        }
    }
    vector<hsize_t> dims;
    dims.push_back(hexs.size());
    dims.push_back(6);
    H5::DataTypeSharedPtr tp = H5::DataType::OfObject(els[0]);
    H5::DataSpaceSharedPtr ds =
        std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    H5::DataSetSharedPtr dst = m_mesh->CreateDataSet("hex", tp, ds);
    dst->Write(els, ds);

    tp = H5::DataType::OfObject(i_map[0]);
    dims.clear();
    dims.push_back(hexs.size());
    dims.push_back(2);
    ds  = std::shared_ptr<H5::DataSpace>(new H5::DataSpace(dims));
    dst = m_maps->CreateDataSet("hex", tp, ds);
    dst->Write(i_map, ds);
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

    std::vector<unsigned int> idxList;
    int i = 0;
    for (auto cIt = composites.begin(); cIt != composites.end();
         ++cIt, i++)
    {

        if (cIt->second->m_geomVec.size() == 0)
        {
            continue;
        }

        idxList.clear();

        for (int i = 0; i < cIt->second->m_geomVec.size(); ++i)
        {
            idxList.push_back(cIt->second->m_geomVec[i]->GetGlobalID());
        }

        comps.push_back(ParseUtils::GenerateSeqString(idxList));
        c_map.push_back(cIt->first);
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

    WriteVertices(m_vertSet);
    WriteEdges(m_segGeoms);
    if (m_meshDimension > 1)
    {
        WriteTris(m_triGeoms);
        WriteQuads(m_quadGeoms);
    }
    if (m_meshDimension > 2)
    {
        WriteTets(m_tetGeoms);
        WritePyrs(m_pyrGeoms);
        WritePrisms(m_prismGeoms);
        WriteHexs(m_hexGeoms);
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
