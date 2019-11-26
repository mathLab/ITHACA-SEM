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

#include "MeshGraphXmlCompressed.h"

#include <LibUtilities/BasicUtils/CompressData.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/CompressData.h>

#include <LibUtilities/Interpreter/Interpreter.h>
#include <SpatialDomains/MeshEntities.hpp>

// These are required for the Write(...) and Import(...) functions.
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/make_shared.hpp>

#include <tinyxml.h>
using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

std::string MeshGraphXmlCompressed::className =
    GetMeshGraphFactory().RegisterCreatorFunction(
        "XmlCompressed", MeshGraphXmlCompressed::create,
        "IO with Xml geometry");

void MeshGraphXmlCompressed::ReadVertices()
{
    // Now read the vertices
    TiXmlElement *element = m_xmlGeom->FirstChildElement("VERTEX");
    ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");

    NekDouble xscale, yscale, zscale;

    // check to see if any scaling parameters are in
    // attributes and determine these values
    LibUtilities::Interpreter expEvaluator;
    const char *xscal = element->Attribute("XSCALE");
    if (!xscal)
    {
        xscale = 1.0;
    }
    else
    {
        std::string xscalstr = xscal;
        int expr_id          = expEvaluator.DefineFunction("", xscalstr);
        xscale               = expEvaluator.Evaluate(expr_id);
    }

    const char *yscal = element->Attribute("YSCALE");
    if (!yscal)
    {
        yscale = 1.0;
    }
    else
    {
        std::string yscalstr = yscal;
        int expr_id          = expEvaluator.DefineFunction("", yscalstr);
        yscale               = expEvaluator.Evaluate(expr_id);
    }

    const char *zscal = element->Attribute("ZSCALE");
    if (!zscal)
    {
        zscale = 1.0;
    }
    else
    {
        std::string zscalstr = zscal;
        int expr_id          = expEvaluator.DefineFunction("", zscalstr);
        zscale               = expEvaluator.Evaluate(expr_id);
    }

    NekDouble xmove, ymove, zmove;

    // check to see if any moving parameters are in
    // attributes and determine these values

    const char *xmov = element->Attribute("XMOVE");
    if (!xmov)
    {
        xmove = 0.0;
    }
    else
    {
        std::string xmovstr = xmov;
        int expr_id         = expEvaluator.DefineFunction("", xmovstr);
        xmove               = expEvaluator.Evaluate(expr_id);
    }

    const char *ymov = element->Attribute("YMOVE");
    if (!ymov)
    {
        ymove = 0.0;
    }
    else
    {
        std::string ymovstr = ymov;
        int expr_id         = expEvaluator.DefineFunction("", ymovstr);
        ymove               = expEvaluator.Evaluate(expr_id);
    }

    const char *zmov = element->Attribute("ZMOVE");
    if (!zmov)
    {
        zmove = 0.0;
    }
    else
    {
        std::string zmovstr = zmov;
        int expr_id         = expEvaluator.DefineFunction("", zmovstr);
        zmove               = expEvaluator.Evaluate(expr_id);
    }

    string IsCompressed;
    element->QueryStringAttribute("COMPRESSED", &IsCompressed);

    if (boost::iequals(IsCompressed,
                       LibUtilities::CompressData::GetCompressString()))
    {
        // Extract the vertex body
        TiXmlNode *vertexChild = element->FirstChild();
        ASSERTL0(vertexChild, "Unable to extract the data from the compressed "
                              "vertex tag.");

        std::string vertexStr;
        if (vertexChild->Type() == TiXmlNode::TINYXML_TEXT)
        {
            vertexStr += vertexChild->ToText()->ValueStr();
        }

        std::vector<SpatialDomains::MeshVertex> vertData;
        LibUtilities::CompressData::ZlibDecodeFromBase64Str(vertexStr,
                                                            vertData);

        int indx;
        NekDouble xval, yval, zval;
        for (int i = 0; i < vertData.size(); ++i)
        {
            indx = vertData[i].id;
            xval = vertData[i].x;
            yval = vertData[i].y;
            zval = vertData[i].z;

            xval = xval * xscale + xmove;
            yval = yval * yscale + ymove;
            zval = zval * zscale + zmove;

            PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(
                m_spaceDimension, indx, xval, yval, zval));

            vert->SetGlobalID(indx);
            m_vertSet[indx] = vert;
        }
    }
    else
    {
        ASSERTL0(false, "Compressed formats do not match. Expected :" +
                            LibUtilities::CompressData::GetCompressString() +
                            " but got " + std::string(IsCompressed));
    }
}

void MeshGraphXmlCompressed::ReadCurves()
{
    // check to see if any scaling parameters are in
    // attributes and determine these values
    TiXmlElement *element = m_xmlGeom->FirstChildElement("VERTEX");
    ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");

    NekDouble xscale, yscale, zscale;

    LibUtilities::Interpreter expEvaluator;
    const char *xscal = element->Attribute("XSCALE");
    if (!xscal)
    {
        xscale = 1.0;
    }
    else
    {
        std::string xscalstr = xscal;
        int expr_id          = expEvaluator.DefineFunction("", xscalstr);
        xscale               = expEvaluator.Evaluate(expr_id);
    }

    const char *yscal = element->Attribute("YSCALE");
    if (!yscal)
    {
        yscale = 1.0;
    }
    else
    {
        std::string yscalstr = yscal;
        int expr_id          = expEvaluator.DefineFunction("", yscalstr);
        yscale               = expEvaluator.Evaluate(expr_id);
    }

    const char *zscal = element->Attribute("ZSCALE");
    if (!zscal)
    {
        zscale = 1.0;
    }
    else
    {
        std::string zscalstr = zscal;
        int expr_id          = expEvaluator.DefineFunction("", zscalstr);
        zscale               = expEvaluator.Evaluate(expr_id);
    }

    NekDouble xmove, ymove, zmove;

    // check to see if any moving parameters are in
    // attributes and determine these values

    const char *xmov = element->Attribute("XMOVE");
    if (!xmov)
    {
        xmove = 0.0;
    }
    else
    {
        std::string xmovstr = xmov;
        int expr_id         = expEvaluator.DefineFunction("", xmovstr);
        xmove               = expEvaluator.Evaluate(expr_id);
    }

    const char *ymov = element->Attribute("YMOVE");
    if (!ymov)
    {
        ymove = 0.0;
    }
    else
    {
        std::string ymovstr = ymov;
        int expr_id         = expEvaluator.DefineFunction("", ymovstr);
        ymove               = expEvaluator.Evaluate(expr_id);
    }

    const char *zmov = element->Attribute("ZMOVE");
    if (!zmov)
    {
        zmove = 0.0;
    }
    else
    {
        std::string zmovstr = zmov;
        int expr_id         = expEvaluator.DefineFunction("", zmovstr);
        zmove               = expEvaluator.Evaluate(expr_id);
    }

    /// Look for elements in CURVE block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("CURVED");

    if (!field) // return if no curved entities
    {
        return;
    }

    string IsCompressed;
    field->QueryStringAttribute("COMPRESSED", &IsCompressed);

    if (IsCompressed.size() == 0)
    {
        // this could be that the curved tag is empty
        // in this case we dont want to read it
        return;
    }

    ASSERTL0(boost::iequals(IsCompressed,
                            LibUtilities::CompressData::GetCompressString()),
             "Compressed formats do not match. Expected :" +
                 LibUtilities::CompressData::GetCompressString() + " but got " +
                 boost::lexical_cast<std::string>(IsCompressed));

    std::vector<SpatialDomains::MeshCurvedInfo> edginfo;
    std::vector<SpatialDomains::MeshCurvedInfo> facinfo;
    SpatialDomains::MeshCurvedPts cpts;

    // read edge, face info and curved poitns.
    TiXmlElement *x = field->FirstChildElement();
    while (x)
    {
        const char *entitytype = x->Value();
        // read in edge or face info
        if (boost::iequals(entitytype, "E"))
        {
            // read in data
            std::string elmtStr;
            TiXmlNode *child = x->FirstChild();

            if (child->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elmtStr += child->ToText()->ValueStr();
            }

            LibUtilities::CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                edginfo);
        }
        else if (boost::iequals(entitytype, "F"))
        {
            // read in data
            std::string elmtStr;
            TiXmlNode *child = x->FirstChild();

            if (child->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elmtStr += child->ToText()->ValueStr();
            }

            LibUtilities::CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                facinfo);
        }
        else if (boost::iequals(entitytype, "DATAPOINTS"))
        {
            NekInt id;
            ASSERTL0(x->Attribute("ID", &id),
                     "Failed to get ID from PTS section");
            cpts.id = id;

            // read in data
            std::string elmtStr;

            TiXmlElement *DataIdx = x->FirstChildElement("INDEX");
            ASSERTL0(DataIdx, "Cannot read data index tag in compressed "
                              "curved section");

            TiXmlNode *child = DataIdx->FirstChild();
            if (child->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elmtStr = child->ToText()->ValueStr();
            }

            LibUtilities::CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                cpts.index);

            TiXmlElement *DataPts = x->FirstChildElement("POINTS");
            ASSERTL0(DataPts, "Cannot read data pts tag in compressed "
                              "curved section");

            child = DataPts->FirstChild();
            if (child->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elmtStr = child->ToText()->ValueStr();
            }

            LibUtilities::CompressData::ZlibDecodeFromBase64Str(elmtStr,
                                                                cpts.pts);
        }
        else
        {
            ASSERTL0(false, "Unknown tag in curved section");
        }
        x = x->NextSiblingElement();
    }

    // rescale (x,y,z) points;
    for (int i = 0; i < cpts.pts.size(); ++i)
    {
        cpts.pts[i].x = xscale * cpts.pts[i].x + xmove;
        cpts.pts[i].y = yscale * cpts.pts[i].y + ymove;
        cpts.pts[i].z = zscale * cpts.pts[i].z + zmove;
    }

    for (int i = 0; i < edginfo.size(); ++i)
    {
        int edgeid = edginfo[i].entityid;
        LibUtilities::PointsType ptype;

        CurveSharedPtr curve(MemoryManager<Curve>::AllocateSharedPtr(
            edgeid, ptype = (LibUtilities::PointsType)edginfo[i].ptype));

        // load points
        int offset = edginfo[i].ptoffset;
        for (int j = 0; j < edginfo[i].npoints; ++j)
        {
            int idx = cpts.index[offset + j];

            PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(
                m_spaceDimension, edginfo[i].id, cpts.pts[idx].x,
                cpts.pts[idx].y, cpts.pts[idx].z));
            curve->m_points.push_back(vert);
        }

        m_curvedEdges[edgeid] = curve;
    }

    for (int i = 0; i < facinfo.size(); ++i)
    {
        int faceid = facinfo[i].entityid;
        LibUtilities::PointsType ptype;

        CurveSharedPtr curve(MemoryManager<Curve>::AllocateSharedPtr(
            faceid, ptype = (LibUtilities::PointsType)facinfo[i].ptype));

        int offset = facinfo[i].ptoffset;
        for (int j = 0; j < facinfo[i].npoints; ++j)
        {
            int idx = cpts.index[offset + j];

            PointGeomSharedPtr vert(MemoryManager<PointGeom>::AllocateSharedPtr(
                m_spaceDimension, facinfo[i].id, cpts.pts[idx].x,
                cpts.pts[idx].y, cpts.pts[idx].z));
            curve->m_points.push_back(vert);
        }

        m_curvedFaces[faceid] = curve;
    }
}

void MeshGraphXmlCompressed::ReadEdges()
{
    CurveMap::iterator it;

    /// Look for elements in ELEMENT block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("EDGE");

    ASSERTL0(field, "Unable to find EDGE tag in file.");

    string IsCompressed;
    field->QueryStringAttribute("COMPRESSED", &IsCompressed);

    ASSERTL0(boost::iequals(IsCompressed,
                            LibUtilities::CompressData::GetCompressString()),
             "Compressed formats do not match. Expected :" +
                 LibUtilities::CompressData::GetCompressString() + " but got " +
                 std::string(IsCompressed));
    // Extract the edge body
    TiXmlNode *edgeChild = field->FirstChild();
    ASSERTL0(edgeChild, "Unable to extract the data from "
                        "the compressed edge tag.");

    std::string edgeStr;
    if (edgeChild->Type() == TiXmlNode::TINYXML_TEXT)
    {
        edgeStr += edgeChild->ToText()->ValueStr();
    }

    std::vector<SpatialDomains::MeshEdge> edgeData;
    LibUtilities::CompressData::ZlibDecodeFromBase64Str(edgeStr, edgeData);

    int indx;
    for (int i = 0; i < edgeData.size(); ++i)
    {
        indx                           = edgeData[i].id;
        PointGeomSharedPtr vertices[2] = {GetVertex(edgeData[i].v0),
                                          GetVertex(edgeData[i].v1)};
        SegGeomSharedPtr edge;

        it = m_curvedEdges.find(indx);
        if (it == m_curvedEdges.end())
        {
            edge = MemoryManager<SegGeom>::AllocateSharedPtr(
                indx, m_spaceDimension, vertices);
        }
        else
        {
            edge = MemoryManager<SegGeom>::AllocateSharedPtr(
                indx, m_spaceDimension, vertices, it->second);
        }
        m_segGeoms[indx] = edge;
    }
}

void MeshGraphXmlCompressed::ReadFaces()
{
    /// Look for elements in FACE block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("FACE");

    ASSERTL0(field, "Unable to find FACE tag in file.");

    /// All faces are of the form: "<? ID="#"> ... </?>", with
    /// ? being an element type (either Q or T).
    /// They might be in compressed format and so then need upacking.

    TiXmlElement *element = field->FirstChildElement();
    CurveMap::iterator it;

    while (element)
    {
        std::string elementType(element->ValueStr());

        ASSERTL0(elementType == "Q" || elementType == "T",
                 (std::string("Unknown 3D face type: ") + elementType).c_str());

        string IsCompressed;
        element->QueryStringAttribute("COMPRESSED", &IsCompressed);

        ASSERTL0(
            boost::iequals(IsCompressed,
                           LibUtilities::CompressData::GetCompressString()),
            "Compressed formats do not match. Expected :" +
                LibUtilities::CompressData::GetCompressString() + " but got " +
                std::string(IsCompressed));

        // Extract the face body
        TiXmlNode *faceChild = element->FirstChild();
        ASSERTL0(faceChild, "Unable to extract the data from "
                            "the compressed face tag.");

        std::string faceStr;
        if (faceChild->Type() == TiXmlNode::TINYXML_TEXT)
        {
            faceStr += faceChild->ToText()->ValueStr();
        }

        int indx;
        if (elementType == "T")
        {
            std::vector<SpatialDomains::MeshTri> faceData;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                faceData);

            for (int i = 0; i < faceData.size(); ++i)
            {
                indx = faceData[i].id;

                /// See if this face has curves.
                it = m_curvedFaces.find(indx);

                /// Create a TriGeom to hold the new definition.
                SegGeomSharedPtr edges[TriGeom::kNedges] = {
                    GetSegGeom(faceData[i].e[0]), GetSegGeom(faceData[i].e[1]),
                    GetSegGeom(faceData[i].e[2])};

                TriGeomSharedPtr trigeom;
                if (it == m_curvedFaces.end())
                {
                    trigeom =
                        MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges);
                }
                else
                {
                    trigeom = MemoryManager<TriGeom>::AllocateSharedPtr(
                        indx, edges, it->second);
                }
                trigeom->SetGlobalID(indx);
                m_triGeoms[indx] = trigeom;
            }
        }
        else if (elementType == "Q")
        {
            std::vector<SpatialDomains::MeshQuad> faceData;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                faceData);

            for (int i = 0; i < faceData.size(); ++i)
            {
                indx = faceData[i].id;

                /// See if this face has curves.
                it = m_curvedFaces.find(indx);

                /// Create a QuadGeom to hold the new definition.
                SegGeomSharedPtr edges[QuadGeom::kNedges] = {
                    GetSegGeom(faceData[i].e[0]), GetSegGeom(faceData[i].e[1]),
                    GetSegGeom(faceData[i].e[2]), GetSegGeom(faceData[i].e[3])};

                QuadGeomSharedPtr quadgeom;
                if (it == m_curvedFaces.end())
                {
                    quadgeom =
                        MemoryManager<QuadGeom>::AllocateSharedPtr(indx, edges);
                }
                else
                {
                    quadgeom = MemoryManager<QuadGeom>::AllocateSharedPtr(
                        indx, edges, it->second);
                }
                quadgeom->SetGlobalID(indx);
                m_quadGeoms[indx] = quadgeom;
            }
        }
        /// Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXmlCompressed::ReadElements1D()
{
    TiXmlElement *field = NULL;

    /// Look for elements in ELEMENT block.
    field = m_xmlGeom->FirstChildElement("ELEMENT");

    ASSERTL0(field, "Unable to find ELEMENT tag in file.");

    /// All elements are of the form: "<S ID = n> ... </S>", with
    /// ? being the element type.

    TiXmlElement *segment = field->FirstChildElement("S");
    CurveMap::iterator it;

    while (segment)
    {
        string IsCompressed;
        segment->QueryStringAttribute("COMPRESSED", &IsCompressed);
        ASSERTL0(
            boost::iequals(IsCompressed,
                           LibUtilities::CompressData::GetCompressString()),
            "Compressed formats do not match. Expected :" +
                LibUtilities::CompressData::GetCompressString() + " but got " +
                std::string(IsCompressed));

        // Extract the face body
        TiXmlNode *child = segment->FirstChild();
        ASSERTL0(child, "Unable to extract the data from "
                        "the compressed face tag.");

        std::string str;
        if (child->Type() == TiXmlNode::TINYXML_TEXT)
        {
            str += child->ToText()->ValueStr();
        }

        int indx;

        std::vector<SpatialDomains::MeshEdge> data;
        LibUtilities::CompressData::ZlibDecodeFromBase64Str(str, data);

        for (int i = 0; i < data.size(); ++i)
        {
            indx = data[i].id;

            /// See if this face has curves.
            it = m_curvedEdges.find(indx);

            PointGeomSharedPtr vertices[2] = {GetVertex(data[i].v0),
                                              GetVertex(data[i].v1)};
            SegGeomSharedPtr seg;

            if (it == m_curvedEdges.end())
            {
                seg = MemoryManager<SegGeom>::AllocateSharedPtr(
                    indx, m_spaceDimension, vertices);
                seg->SetGlobalID(indx); // Set global mesh id
            }
            else
            {
                seg = MemoryManager<SegGeom>::AllocateSharedPtr(
                    indx, m_spaceDimension, vertices, it->second);
                seg->SetGlobalID(indx); // Set global mesh id
            }
            seg->SetGlobalID(indx);
            m_segGeoms[indx] = seg;
        }
        /// Keep looking for additional segments
        segment = segment->NextSiblingElement("S");
    }
}

void MeshGraphXmlCompressed::ReadElements2D()
{
    /// Look for elements in ELEMENT block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("ELEMENT");

    ASSERTL0(field, "Unable to find ELEMENT tag in file.");

    // Set up curve map for curved elements on an embedded manifold.
    CurveMap::iterator it;

    /// All elements are of the form: "<? ID="#"> ... </?>", with
    /// ? being the element type.

    TiXmlElement *element = field->FirstChildElement();

    while (element)
    {
        std::string elementType(element->ValueStr());

        ASSERTL0(
            elementType == "Q" || elementType == "T",
            (std::string("Unknown 2D element type: ") + elementType).c_str());

        string IsCompressed;
        element->QueryStringAttribute("COMPRESSED", &IsCompressed);

        ASSERTL0(
            boost::iequals(IsCompressed,
                           LibUtilities::CompressData::GetCompressString()),
            "Compressed formats do not match. Expected :" +
                LibUtilities::CompressData::GetCompressString() + " but got " +
                std::string(IsCompressed));

        // Extract the face body
        TiXmlNode *faceChild = element->FirstChild();
        ASSERTL0(faceChild, "Unable to extract the data from "
                            "the compressed face tag.");

        std::string faceStr;
        if (faceChild->Type() == TiXmlNode::TINYXML_TEXT)
        {
            faceStr += faceChild->ToText()->ValueStr();
        }

        int indx;
        if (elementType == "T")
        {
            std::vector<SpatialDomains::MeshTri> faceData;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                faceData);

            for (int i = 0; i < faceData.size(); ++i)
            {
                indx = faceData[i].id;

                /// See if this face has curves.
                it = m_curvedFaces.find(indx);

                /// Create a TriGeom to hold the new definition.
                SegGeomSharedPtr edges[TriGeom::kNedges] = {
                    GetSegGeom(faceData[i].e[0]), GetSegGeom(faceData[i].e[1]),
                    GetSegGeom(faceData[i].e[2])};

                TriGeomSharedPtr trigeom;
                if (it == m_curvedFaces.end())
                {
                    trigeom =
                        MemoryManager<TriGeom>::AllocateSharedPtr(indx, edges);
                }
                else
                {
                    trigeom = MemoryManager<TriGeom>::AllocateSharedPtr(
                        indx, edges, it->second);
                }
                trigeom->SetGlobalID(indx);
                m_triGeoms[indx] = trigeom;
            }
        }
        else if (elementType == "Q")
        {
            std::vector<SpatialDomains::MeshQuad> faceData;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(faceStr,
                                                                faceData);

            for (int i = 0; i < faceData.size(); ++i)
            {
                indx = faceData[i].id;

                /// See if this face has curves.
                it = m_curvedFaces.find(indx);

                /// Create a QuadGeom to hold the new definition.
                SegGeomSharedPtr edges[QuadGeom::kNedges] = {
                    GetSegGeom(faceData[i].e[0]), GetSegGeom(faceData[i].e[1]),
                    GetSegGeom(faceData[i].e[2]), GetSegGeom(faceData[i].e[3])};

                QuadGeomSharedPtr quadgeom;
                if (it == m_curvedFaces.end())
                {
                    quadgeom =
                        MemoryManager<QuadGeom>::AllocateSharedPtr(indx, edges);
                }
                else
                {
                    quadgeom = MemoryManager<QuadGeom>::AllocateSharedPtr(
                        indx, edges, it->second);
                }
                quadgeom->SetGlobalID(indx);
                m_quadGeoms[indx] = quadgeom;
            }
        }
        /// Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXmlCompressed::ReadElements3D()
{
    /// Look for elements in ELEMENT block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("ELEMENT");

    ASSERTL0(field, "Unable to find ELEMENT tag in file.");

    /// All elements are of the form: "<? ID="#"> ... </?>", with
    /// ? being the element type.

    TiXmlElement *element = field->FirstChildElement();

    while (element)
    {
        std::string elementType(element->ValueStr());

        // A - tet, P - pyramid, R - prism, H - hex
        ASSERTL0(
            elementType == "A" || elementType == "P" || elementType == "R" ||
                elementType == "H",
            (std::string("Unknown 3D element type: ") + elementType).c_str());

        string IsCompressed;
        element->QueryStringAttribute("COMPRESSED", &IsCompressed);

        ASSERTL0(
            boost::iequals(IsCompressed,
                           LibUtilities::CompressData::GetCompressString()),
            "Compressed formats do not match. Expected :" +
                LibUtilities::CompressData::GetCompressString() + " but got " +
                std::string(IsCompressed));

        // Extract the face body
        TiXmlNode *child = element->FirstChild();
        ASSERTL0(child, "Unable to extract the data from "
                        "the compressed face tag.");

        std::string str;
        if (child->Type() == TiXmlNode::TINYXML_TEXT)
        {
            str += child->ToText()->ValueStr();
        }

        int indx;
        if (elementType == "A")
        {
            std::vector<SpatialDomains::MeshTet> data;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(str, data);
            TriGeomSharedPtr tfaces[4];
            for (int i = 0; i < data.size(); ++i)
            {
                indx = data[i].id;
                for (int j = 0; j < 4; ++j)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(data[i].f[j]);
                    tfaces[j] = static_pointer_cast<TriGeom>(face);
                }

                TetGeomSharedPtr tetgeom(
                    MemoryManager<TetGeom>::AllocateSharedPtr(indx, tfaces));
                m_tetGeoms[indx] = tetgeom;
                PopulateFaceToElMap(tetgeom, 4);
            }
        }
        else if (elementType == "P")
        {
            std::vector<SpatialDomains::MeshPyr> data;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(str, data);
            Geometry2DSharedPtr faces[5];
            for (int i = 0; i < data.size(); ++i)
            {
                indx        = data[i].id;
                int Ntfaces = 0;
                int Nqfaces = 0;
                for (int j = 0; j < 5; ++j)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(data[i].f[j]);

                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << j;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        faces[j] = static_pointer_cast<TriGeom>(face);
                        Ntfaces++;
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        faces[j] = static_pointer_cast<QuadGeom>(face);
                        Nqfaces++;
                    }
                }
                ASSERTL0((Ntfaces == 4) && (Nqfaces == 1),
                         "Did not identify the correct number of "
                         "triangular and quadrilateral faces for a "
                         "pyramid");

                PyrGeomSharedPtr pyrgeom(
                    MemoryManager<PyrGeom>::AllocateSharedPtr(indx, faces));

                m_pyrGeoms[indx] = pyrgeom;
                PopulateFaceToElMap(pyrgeom, 5);
            }
        }
        else if (elementType == "R")
        {
            std::vector<SpatialDomains::MeshPrism> data;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(str, data);
            Geometry2DSharedPtr faces[5];
            for (int i = 0; i < data.size(); ++i)
            {
                indx        = data[i].id;
                int Ntfaces = 0;
                int Nqfaces = 0;
                for (int j = 0; j < 5; ++j)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(data[i].f[j]);
                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << j;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        faces[j] = static_pointer_cast<TriGeom>(face);
                        Ntfaces++;
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        faces[j] = static_pointer_cast<QuadGeom>(face);
                        Nqfaces++;
                    }
                }
                ASSERTL0((Ntfaces == 2) && (Nqfaces == 3),
                         "Did not identify the correct number of "
                         "triangular and quadrilateral faces for a "
                         "prism");

                PrismGeomSharedPtr prismgeom(
                    MemoryManager<PrismGeom>::AllocateSharedPtr(indx, faces));

                m_prismGeoms[indx] = prismgeom;
                PopulateFaceToElMap(prismgeom, 5);
            }
        }
        else if (elementType == "H")
        {
            std::vector<SpatialDomains::MeshHex> data;
            LibUtilities::CompressData::ZlibDecodeFromBase64Str(str, data);

            QuadGeomSharedPtr faces[6];
            for (int i = 0; i < data.size(); ++i)
            {
                indx = data[i].id;
                for (int j = 0; j < 6; ++j)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(data[i].f[j]);
                    faces[j] = static_pointer_cast<QuadGeom>(face);
                }

                HexGeomSharedPtr hexgeom(
                    MemoryManager<HexGeom>::AllocateSharedPtr(indx, faces));
                m_hexGeoms[indx] = hexgeom;
                PopulateFaceToElMap(hexgeom, 6);
            }
        }
        /// Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXmlCompressed::WriteVertices(TiXmlElement *geomTag,
                                           PointGeomMap &verts)
{
    if (verts.size() == 0)
    {
        return;
    }

    TiXmlElement *vertTag = new TiXmlElement("VERTEX");

    vector<MeshVertex> vertInfo;

    for (auto &i : verts)
    {
        MeshVertex v;
        v.id = i.first;
        v.x  = i.second->x();
        v.y  = i.second->y();
        v.z  = i.second->z();
        vertInfo.push_back(v);
    }

    vertTag->SetAttribute("COMPRESSED",
                          LibUtilities::CompressData::GetCompressString());
    vertTag->SetAttribute("BITSIZE",
                          LibUtilities::CompressData::GetBitSizeStr());

    string vertStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(vertInfo, vertStr);

    vertTag->LinkEndChild(new TiXmlText(vertStr));

    geomTag->LinkEndChild(vertTag);
}

void MeshGraphXmlCompressed::WriteEdges(TiXmlElement *geomTag,
                                        SegGeomMap &edges)
{
    if (edges.size() == 0)
    {
        return;
    }

    TiXmlElement *edgeTag =
        new TiXmlElement(m_meshDimension == 1 ? "S" : "EDGE");

    vector<MeshEdge> edgeInfo;

    for (auto &i : edges)
    {
        MeshEdge e;
        e.id = i.first;
        e.v0 = i.second->GetVid(0);
        e.v1 = i.second->GetVid(1);
        edgeInfo.push_back(e);
    }

    string edgeStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(edgeInfo, edgeStr);

    edgeTag->SetAttribute("COMPRESSED",
                          LibUtilities::CompressData::GetCompressString());
    edgeTag->SetAttribute("BITSIZE",
                          LibUtilities::CompressData::GetBitSizeStr());

    edgeTag->LinkEndChild(new TiXmlText(edgeStr));

    if (m_meshDimension == 1)
    {
        TiXmlElement *tmp = new TiXmlElement("ELEMENT");
        tmp->LinkEndChild(edgeTag);
        geomTag->LinkEndChild(tmp);
    }
    else
    {
        geomTag->LinkEndChild(edgeTag);
    }
}

void MeshGraphXmlCompressed::WriteTris(TiXmlElement *faceTag, TriGeomMap &tris)
{
    if (tris.size() == 0)
    {
        return;
    }

    string tag = "T";

    vector<MeshTri> triInfo;

    for (auto &i : tris)
    {
        MeshTri t;
        t.id   = i.first;
        t.e[0] = i.second->GetEid(0);
        t.e[1] = i.second->GetEid(1);
        t.e[2] = i.second->GetEid(2);
        triInfo.push_back(t);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string triStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(triInfo, triStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(triStr));

    faceTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WriteQuads(TiXmlElement *faceTag,
                                        QuadGeomMap &quads)
{
    if (quads.size() == 0)
    {
        return;
    }

    string tag = "Q";

    vector<MeshQuad> quadInfo;

    for (auto &i : quads)
    {
        MeshQuad q;
        q.id   = i.first;
        q.e[0] = i.second->GetEid(0);
        q.e[1] = i.second->GetEid(1);
        q.e[2] = i.second->GetEid(2);
        q.e[3] = i.second->GetEid(3);
        quadInfo.push_back(q);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string quadStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(quadInfo, quadStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(quadStr));

    faceTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WriteHexs(TiXmlElement *elmtTag, HexGeomMap &hexs)
{
    if (hexs.size() == 0)
    {
        return;
    }

    string tag = "H";

    vector<MeshHex> elementInfo;

    for (auto &i : hexs)
    {
        MeshHex e;
        e.id   = i.first;
        e.f[0] = i.second->GetFid(0);
        e.f[1] = i.second->GetFid(1);
        e.f[2] = i.second->GetFid(2);
        e.f[3] = i.second->GetFid(3);
        e.f[4] = i.second->GetFid(4);
        e.f[5] = i.second->GetFid(5);
        elementInfo.push_back(e);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string elStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(elementInfo, elStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(elStr));

    elmtTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WritePrisms(TiXmlElement *elmtTag,
                                         PrismGeomMap &pris)
{
    if (pris.size() == 0)
    {
        return;
    }

    string tag = "R";

    vector<MeshPrism> elementInfo;

    for (auto &i : pris)
    {
        MeshPrism e;
        e.id   = i.first;
        e.f[0] = i.second->GetFid(0);
        e.f[1] = i.second->GetFid(1);
        e.f[2] = i.second->GetFid(2);
        e.f[3] = i.second->GetFid(3);
        e.f[4] = i.second->GetFid(4);
        elementInfo.push_back(e);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string elStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(elementInfo, elStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(elStr));

    elmtTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WritePyrs(TiXmlElement *elmtTag, PyrGeomMap &pyrs)
{
    if (pyrs.size() == 0)
    {
        return;
    }

    string tag = "P";

    vector<MeshPyr> elementInfo;

    for (auto &i : pyrs)
    {
        MeshPyr e;
        e.id   = i.first;
        e.f[0] = i.second->GetFid(0);
        e.f[1] = i.second->GetFid(1);
        e.f[2] = i.second->GetFid(2);
        e.f[3] = i.second->GetFid(3);
        e.f[4] = i.second->GetFid(4);
        elementInfo.push_back(e);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string elStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(elementInfo, elStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(elStr));

    elmtTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WriteTets(TiXmlElement *elmtTag, TetGeomMap &tets)
{
    if (tets.size() == 0)
    {
        return;
    }

    string tag = "A";

    vector<MeshTet> elementInfo;

    for (auto &i : tets)
    {
        MeshTet e;
        e.id   = i.first;
        e.f[0] = i.second->GetFid(0);
        e.f[1] = i.second->GetFid(1);
        e.f[2] = i.second->GetFid(2);
        e.f[3] = i.second->GetFid(3);
        elementInfo.push_back(e);
    }

    TiXmlElement *x = new TiXmlElement(tag);
    string elStr;
    LibUtilities::CompressData::ZlibEncodeToBase64Str(elementInfo, elStr);

    x->SetAttribute("COMPRESSED",
                    LibUtilities::CompressData::GetCompressString());
    x->SetAttribute("BITSIZE", LibUtilities::CompressData::GetBitSizeStr());

    x->LinkEndChild(new TiXmlText(elStr));

    elmtTag->LinkEndChild(x);
}

void MeshGraphXmlCompressed::WriteCurves(TiXmlElement *geomTag, CurveMap &edges,
                                         CurveMap &faces)
{
    if (edges.size() == 0 && faces.size() == 0)
    {
        return;
    }

    TiXmlElement *curveTag = new TiXmlElement("CURVED");

    vector<MeshCurvedInfo> edgeInfo;
    vector<MeshCurvedInfo> faceInfo;
    MeshCurvedPts curvedPts;
    curvedPts.id = 0;
    int ptOffset = 0;
    int newIdx   = 0;
    int edgeCnt  = 0;
    int faceCnt  = 0;

    for (auto &i : edges)
    {
        MeshCurvedInfo cinfo;
        cinfo.id       = edgeCnt++;
        cinfo.entityid = i.first;
        cinfo.npoints  = i.second->m_points.size();
        cinfo.ptype    = i.second->m_ptype;
        cinfo.ptid     = 0;
        cinfo.ptoffset = ptOffset;

        edgeInfo.push_back(cinfo);

        for (int j = 0; j < i.second->m_points.size(); j++)
        {
            MeshVertex v;
            v.id = newIdx;
            v.x  = i.second->m_points[j]->x();
            v.y  = i.second->m_points[j]->y();
            v.z  = i.second->m_points[j]->z();
            curvedPts.pts.push_back(v);
            curvedPts.index.push_back(newIdx);
            newIdx++;
        }
        ptOffset += cinfo.npoints;
    }

    for (auto &i : faces)
    {
        MeshCurvedInfo cinfo;
        cinfo.id       = faceCnt++;
        cinfo.entityid = i.first;
        cinfo.npoints  = i.second->m_points.size();
        cinfo.ptype    = i.second->m_ptype;
        cinfo.ptid     = 0;
        cinfo.ptoffset = ptOffset;

        faceInfo.push_back(cinfo);

        for (int j = 0; j < i.second->m_points.size(); j++)
        {
            MeshVertex v;
            v.id = newIdx;
            v.x  = i.second->m_points[j]->x();
            v.y  = i.second->m_points[j]->y();
            v.z  = i.second->m_points[j]->z();
            curvedPts.pts.push_back(v);
            curvedPts.index.push_back(newIdx);
            newIdx++;
        }
        ptOffset += cinfo.npoints;
    }

    curveTag->SetAttribute("COMPRESSED",
                           LibUtilities::CompressData::GetCompressString());
    curveTag->SetAttribute("BITSIZE",
                           LibUtilities::CompressData::GetBitSizeStr());

    if (edgeInfo.size())
    {
        TiXmlElement *x = new TiXmlElement("E");
        string dataStr;
        LibUtilities::CompressData::ZlibEncodeToBase64Str(edgeInfo, dataStr);

        x->LinkEndChild(new TiXmlText(dataStr));
        curveTag->LinkEndChild(x);
    }

    if (faceInfo.size())
    {
        TiXmlElement *x = new TiXmlElement("F");
        string dataStr;
        LibUtilities::CompressData::ZlibEncodeToBase64Str(faceInfo, dataStr);

        x->LinkEndChild(new TiXmlText(dataStr));
        curveTag->LinkEndChild(x);
    }

    if (edgeInfo.size() || faceInfo.size())
    {
        TiXmlElement *x = new TiXmlElement("DATAPOINTS");
        x->SetAttribute("ID", curvedPts.id);
        TiXmlElement *subx = new TiXmlElement("INDEX");
        string dataStr;
        LibUtilities::CompressData::ZlibEncodeToBase64Str(curvedPts.index,
                                                          dataStr);
        subx->LinkEndChild(new TiXmlText(dataStr));
        x->LinkEndChild(subx);

        subx = new TiXmlElement("POINTS");
        LibUtilities::CompressData::ZlibEncodeToBase64Str(curvedPts.pts,
                                                          dataStr);
        subx->LinkEndChild(new TiXmlText(dataStr));
        x->LinkEndChild(subx);
        curveTag->LinkEndChild(x);
    }

    geomTag->LinkEndChild(curveTag);
}
} // namespace SpatialDomains
} // namespace Nektar
