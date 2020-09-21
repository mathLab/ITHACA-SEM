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

#include <iomanip>

#include <SpatialDomains/MeshGraphXml.h>
#include <SpatialDomains/MeshPartition.h>

#include <LibUtilities/Interpreter/Interpreter.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/FieldIOXml.h>

#include <boost/format.hpp>

#include <tinyxml.h>

using namespace std;

namespace Nektar
{
namespace SpatialDomains
{

std::string MeshGraphXml::className =
    GetMeshGraphFactory().RegisterCreatorFunction(
        "Xml", MeshGraphXml::create, "IO with Xml geometry");

void MeshGraphXml::PartitionMesh(
    const LibUtilities::SessionReaderSharedPtr session)
{
    // Get row of comm, or the whole comm if not split
    LibUtilities::CommSharedPtr comm     = session->GetComm();
    LibUtilities::CommSharedPtr commMesh = comm->GetRowComm();
    const bool                  isRoot   = comm->TreatAsRankZero();

    m_session = session;

    // Load file for root process only (since this is always needed)
    // and determine if the provided geometry has already been
    // partitioned. This will be the case if the user provides the
    // directory of mesh partitions as an input. Partitioned geometries
    // have the attribute
    //    PARTITION=X
    // where X is the number of the partition (and should match the
    // process rank). The result is shared with all other processes.
    int isPartitioned = 0;
    if (isRoot)
    {
        if (m_session->DefinesElement("Nektar/Geometry"))
        {
            if (m_session->GetElement("Nektar/Geometry")->Attribute("PARTITION"))
            {
                std::cout << "Using pre-partitioned mesh." << std::endl;
                isPartitioned = 1;
            }
        }
    }
    comm->Bcast(isPartitioned, 0);

    // If the mesh is already partitioned, we are done. Remaining
    // processes must load their partitions.
    if (isPartitioned)
    {
        if (!isRoot)
        {
            m_session->InitSession();
        }
    }
    else
    {
        // Default partitioner to use is Metis. Use Scotch as default if it is
        // installed. Override default with command-line flags if they are set.
        string partitionerName = "Metis";
        if (GetMeshPartitionFactory().ModuleExists("Scotch"))
        {
            partitionerName = "Scotch";
        }
        if (session->DefinesCmdLineArgument("use-metis"))
        {
            partitionerName = "Metis";
        }
        if (session->DefinesCmdLineArgument("use-scotch"))
        {
            partitionerName = "Scotch";
        }

        // Mesh has not been partitioned so do partitioning if required.  Note
        // in the serial case nothing is done as we have already loaded the
        // mesh.
        if (session->DefinesCmdLineArgument("part-only")||
            session->DefinesCmdLineArgument("part-only-overlapping"))
        {
            // Perform partitioning of the mesh only. For this we insist the
            // code is run in serial (parallel execution is pointless).
            ASSERTL0(comm->GetSize() == 1,
                     "The 'part-only' option should be used in serial.");

            // Read 'lite' geometry information
            ReadGeometry(NullDomainRangeShPtr, false);

            // Number of partitions is specified by the parameter.
            int nParts;
            auto comp = CreateCompositeDescriptor();

            MeshPartitionSharedPtr partitioner =
                GetMeshPartitionFactory().CreateInstance(
                    partitionerName, session, m_meshDimension,
                    CreateMeshEntities(), comp);

            if (session->DefinesCmdLineArgument("part-only"))
            {
                nParts = session->GetCmdLineArgument<int>("part-only");
                partitioner->PartitionMesh(nParts, true);
            }
            else
            {
                nParts = session->GetCmdLineArgument<int>("part-only-overlapping");
                partitioner->PartitionMesh(nParts, true, true);
            }

            vector<set<unsigned int>> elmtIDs;
            vector<unsigned int> parts(nParts);
            for (int i = 0; i < nParts; ++i)
            {
                vector<unsigned int> elIDs;
                set<unsigned int> tmp;
                partitioner->GetElementIDs(i, elIDs);
                tmp.insert(elIDs.begin(), elIDs.end());
                elmtIDs.push_back(tmp);
                parts[i] = i;
            }

            this->WriteXMLGeometry(m_session->GetSessionName(), elmtIDs, parts);

            if (isRoot && session->DefinesCmdLineArgument("part-info"))
            {
                partitioner->PrintPartInfo(std::cout);
            }

            session->Finalise();
            exit(0);
        }

        if (commMesh->GetSize() > 1)
        {
            int nParts = commMesh->GetSize();

            if (session->GetSharedFilesystem())
            {
                vector<unsigned int> keys, vals;
                int i;

                if (isRoot)
                {
                    // Read 'lite' geometry information
                    ReadGeometry(NullDomainRangeShPtr, false);

                    // Store composite ordering and boundary information.
                    m_compOrder = CreateCompositeOrdering();
                    auto comp = CreateCompositeDescriptor();

                    // Create mesh partitioner.
                    MeshPartitionSharedPtr partitioner =
                        GetMeshPartitionFactory().CreateInstance(
                            partitionerName, session, m_meshDimension,
                            CreateMeshEntities(), comp);

                    partitioner->PartitionMesh(nParts, true);

                    vector<set<unsigned int>> elmtIDs;
                    vector<unsigned int> parts(nParts);
                    for (i = 0; i < nParts; ++i)
                    {
                        vector<unsigned int> elIDs;
                        set<unsigned int> tmp;
                        partitioner->GetElementIDs(i, elIDs);
                        tmp.insert(elIDs.begin(), elIDs.end());
                        elmtIDs.push_back(tmp);
                        parts[i] = i;
                    }

                    // Call WriteGeometry to write out partition files. This
                    // will populate m_bndRegOrder.
                    this->WriteXMLGeometry(
                        m_session->GetSessionName(), elmtIDs, parts);

                    // Communicate orderings to the other processors.

                    // First send sizes of the orderings and boundary
                    // regions to allocate storage on the remote end.
                    keys.resize(2);
                    keys[0] = m_compOrder.size();
                    keys[1] = m_bndRegOrder.size();
                    comm->Bcast(keys, 0);

                    // Construct the keys and sizes of values for composite
                    // ordering
                    keys.resize(m_compOrder.size());
                    vals.resize(m_compOrder.size());

                    i = 0;
                    for (auto &cIt : m_compOrder)
                    {
                        keys[i  ] = cIt.first;
                        vals[i++] = cIt.second.size();
                    }

                    // Send across data.
                    comm->Bcast(keys, 0);
                    comm->Bcast(vals, 0);
                    for (auto &cIt : m_compOrder)
                    {
                        comm->Bcast(cIt.second, 0);
                    }

                    // Construct the keys and sizes of values for composite
                    // ordering
                    keys.resize(m_bndRegOrder.size());
                    vals.resize(m_bndRegOrder.size());

                    i = 0;
                    for (auto &bIt : m_bndRegOrder)
                    {
                        keys[i  ] = bIt.first;
                        vals[i++] = bIt.second.size();
                    }

                    // Send across data.
                    if (!keys.empty())
                    {
                        comm->Bcast(keys, 0);
                    }
                    if (!vals.empty())
                    {
                        comm->Bcast(vals, 0);
                    }
                    for (auto &bIt : m_bndRegOrder)
                    {
                        comm->Bcast(bIt.second, 0);
                    }

                    if (session->DefinesCmdLineArgument("part-info"))
                    {
                        partitioner->PrintPartInfo(std::cout);
                    }
                }
                else
                {
                    keys.resize(2);
                    comm->Bcast(keys, 0);

                    int cmpSize = keys[0];
                    int bndSize = keys[1];

                    keys.resize(cmpSize);
                    vals.resize(cmpSize);
                    comm->Bcast(keys, 0);
                    comm->Bcast(vals, 0);

                    for (int i = 0; i < keys.size(); ++i)
                    {
                        vector<unsigned int> tmp(vals[i]);
                        comm->Bcast(tmp, 0);
                        m_compOrder[keys[i]] = tmp;
                    }

                    keys.resize(bndSize);
                    vals.resize(bndSize);
                    if (!keys.empty())
                    {
                        comm->Bcast(keys, 0);
                    }
                    if (!vals.empty())
                    {
                        comm->Bcast(vals, 0);
                    }
                    for (int i = 0; i < keys.size(); ++i)
                    {
                        vector<unsigned int> tmp(vals[i]);
                        comm->Bcast(tmp, 0);
                        m_bndRegOrder[keys[i]] = tmp;
                    }
                }
            }
            else
            {
                m_session->InitSession();
                ReadGeometry(NullDomainRangeShPtr, false);

                m_compOrder = CreateCompositeOrdering();
                auto comp = CreateCompositeDescriptor();

                // Partitioner now operates in parallel. Each process receives
                // partitioning over interconnect and writes its own session
                // file to the working directory.
                MeshPartitionSharedPtr partitioner =
                    GetMeshPartitionFactory().CreateInstance(
                        partitionerName, session, m_meshDimension,
                        CreateMeshEntities(), comp);

                partitioner->PartitionMesh(nParts, false);

                vector<unsigned int> parts(1), tmp;
                parts[0] = commMesh->GetRank();
                vector<set<unsigned int>> elIDs(1);
                partitioner->GetElementIDs(parts[0], tmp);
                elIDs[0].insert(tmp.begin(), tmp.end());
                this->WriteXMLGeometry(session->GetSessionName(), elIDs, parts);

                if (m_session->DefinesCmdLineArgument("part-info") && isRoot)
                {
                    partitioner->PrintPartInfo(std::cout);
                }
            }

            // Wait for all processors to finish their writing activities.
            comm->Block();

            std::string  dirname = m_session->GetSessionName() + "_xml";
            fs::path    pdirname(dirname);
            boost::format pad("P%1$07d.xml");
            pad % comm->GetRowComm()->GetRank();
            fs::path    pFilename(pad.str());
            fs::path fullpath = pdirname / pFilename;

            std::vector<std::string> filenames = {
                LibUtilities::PortablePath(fullpath) };
            m_session->InitSession(filenames);
        }
        else if (!isRoot)
        {
            // No partitioning, non-root processors need to read the session
            // file -- handles case where --npz is the same as number of
            // processors.
            m_session->InitSession();
        }
    }
}

void MeshGraphXml::ReadGeometry(
    DomainRangeShPtr rng,
    bool             fillGraph)
{
    // Reset member variables.
    m_vertSet.clear();
    m_curvedEdges.clear();
    m_curvedFaces.clear();
    m_segGeoms.clear();
    m_triGeoms.clear();
    m_quadGeoms.clear();
    m_tetGeoms.clear();
    m_pyrGeoms.clear();
    m_prismGeoms.clear();
    m_hexGeoms.clear();
    m_meshComposites.clear();
    m_compositesLabels.clear();
    m_domain.clear();
    m_expansionMapShPtrMap.clear();
    m_geomInfo.clear();
    m_faceToElMap.clear();

    m_domainRange = rng;
    m_xmlGeom     = m_session->GetElement("NEKTAR/GEOMETRY");

    int err; /// Error value returned by TinyXML.

    TiXmlAttribute *attr = m_xmlGeom->FirstAttribute();

    // Initialize the mesh and space dimensions to 3 dimensions.
    // We want to do this each time we read a file, so it should
    // be done here and not just during class initialization.
    m_meshPartitioned = false;
    m_meshDimension   = 3;
    m_spaceDimension  = 3;

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
            err = attr->QueryIntValue(&m_partition);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read partition.");
            m_meshPartitioned = true;
        }
        else
        {
            std::string errstr("Unknown attribute: ");
            errstr += attrName;
            ASSERTL0(false, errstr.c_str());
        }

        // Get the next attribute.
        attr = attr->Next();
    }

    ASSERTL0(m_meshDimension <= m_spaceDimension,
             "Mesh dimension greater than space dimension");

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

    if (fillGraph)
    {
        MeshGraph::FillGraph();
    }
}

void MeshGraphXml::ReadVertices()
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

    TiXmlElement *vertex = element->FirstChildElement("V");

    int indx;
    int nextVertexNumber = -1;

    while (vertex)
    {
        nextVertexNumber++;

        TiXmlAttribute *vertexAttr = vertex->FirstAttribute();
        std::string attrName(vertexAttr->Name());

        ASSERTL0(attrName == "ID",
                 (std::string("Unknown attribute name: ") + attrName).c_str());

        int err = vertexAttr->QueryIntValue(&indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

        // Now read body of vertex
        std::string vertexBodyStr;

        TiXmlNode *vertexBody = vertex->FirstChild();

        while (vertexBody)
        {
            // Accumulate all non-comment body data.
            if (vertexBody->Type() == TiXmlNode::TINYXML_TEXT)
            {
                vertexBodyStr += vertexBody->ToText()->Value();
                vertexBodyStr += " ";
            }

            vertexBody = vertexBody->NextSibling();
        }

        ASSERTL0(!vertexBodyStr.empty(),
                 "Vertex definitions must contain vertex data.");

        // Get vertex data from the data string.
        NekDouble xval, yval, zval;
        std::istringstream vertexDataStrm(vertexBodyStr.c_str());

        try
        {
            while (!vertexDataStrm.fail())
            {
                vertexDataStrm >> xval >> yval >> zval;

                xval = xval * xscale + xmove;
                yval = yval * yscale + ymove;
                zval = zval * zscale + zmove;

                // Need to check it here because we may not be
                // good after the read indicating that there
                // was nothing to read.
                if (!vertexDataStrm.fail())
                {
                    PointGeomSharedPtr vert(
                        MemoryManager<PointGeom>::AllocateSharedPtr(
                            m_spaceDimension, indx, xval, yval, zval));
                    m_vertSet[indx] = vert;
                }
            }
        }
        catch (...)
        {
            ASSERTL0(false, "Unable to read VERTEX data.");
        }

        vertex = vertex->NextSiblingElement("V");
    }
}

void MeshGraphXml::ReadCurves()
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

    int err;

    /// Look for elements in CURVE block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("CURVED");

    if (!field) // return if no curved entities
    {
        return;
    }

    /// All curves are of the form: "<? ID="#" TYPE="GLL OR other
    /// points type" NUMPOINTS="#"> ... </?>", with ? being an
    /// element type (either E or F).

    TiXmlElement *edgelement = field->FirstChildElement("E");

    int edgeindx, edgeid;
    int nextEdgeNumber = -1;

    while (edgelement)
    {
        /// These should be ordered.
        nextEdgeNumber++;

        std::string edge(edgelement->ValueStr());
        ASSERTL0(edge == "E",
                 (std::string("Unknown 3D curve type:") + edge).c_str());

        /// Read id attribute.
        err = edgelement->QueryIntAttribute("ID", &edgeindx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");

        /// Read edge id attribute.
        err = edgelement->QueryIntAttribute("EDGEID", &edgeid);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Unable to read curve attribute EDGEID.");

        /// Read text edgelement description.
        std::string elementStr;
        TiXmlNode *elementChild = edgelement->FirstChild();

        while (elementChild)
        {
            // Accumulate all non-comment element data
            if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
                elementStr += " ";
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(), "Unable to read curve description body.");

        /// Parse out the element components corresponding to type of
        /// element.
        if (edge == "E")
        {
            int numPts = 0;
            // Determine the points type
            std::string typeStr = edgelement->Attribute("TYPE");
            ASSERTL0(!typeStr.empty(), "TYPE must be specified in "
                                       "points definition");

            LibUtilities::PointsType type;
            const std::string *begStr = LibUtilities::kPointsTypeStr;
            const std::string *endStr =
                LibUtilities::kPointsTypeStr + LibUtilities::SIZE_PointsType;
            const std::string *ptsStr = std::find(begStr, endStr, typeStr);

            ASSERTL0(ptsStr != endStr, "Invalid points type.");
            type = (LibUtilities::PointsType)(ptsStr - begStr);

            // Determine the number of points
            err = edgelement->QueryIntAttribute("NUMPOINTS", &numPts);
            ASSERTL0(err == TIXML_SUCCESS,
                     "Unable to read curve attribute NUMPOINTS.");
            CurveSharedPtr curve(
                MemoryManager<Curve>::AllocateSharedPtr(edgeid, type));

            // Read points (x, y, z)
            NekDouble xval, yval, zval;
            std::istringstream elementDataStrm(elementStr.c_str());
            try
            {
                while (!elementDataStrm.fail())
                {
                    elementDataStrm >> xval >> yval >> zval;

                    xval = xval * xscale + xmove;
                    yval = yval * yscale + ymove;
                    zval = zval * zscale + zmove;

                    // Need to check it here because we may not be
                    // good after the read indicating that there
                    // was nothing to read.
                    if (!elementDataStrm.fail())
                    {
                        PointGeomSharedPtr vert(
                            MemoryManager<PointGeom>::AllocateSharedPtr(
                                m_meshDimension, edgeindx, xval, yval, zval));

                        curve->m_points.push_back(vert);
                    }
                }
            }
            catch (...)
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unable to read curve data for EDGE: ") +
                          elementStr)
                             .c_str());
            }

            ASSERTL0(curve->m_points.size() == numPts,
                     "Number of points specificed by attribute "
                     "NUMPOINTS is different from number of points "
                     "in list (edgeid = " +
                         boost::lexical_cast<string>(edgeid));

            m_curvedEdges[edgeid] = curve;

            edgelement = edgelement->NextSiblingElement("E");

        } // end if-loop

    } // end while-loop

    TiXmlElement *facelement = field->FirstChildElement("F");
    int faceindx, faceid;

    while (facelement)
    {
        std::string face(facelement->ValueStr());
        ASSERTL0(face == "F",
                 (std::string("Unknown 3D curve type: ") + face).c_str());

        /// Read id attribute.
        err = facelement->QueryIntAttribute("ID", &faceindx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read curve attribute ID.");

        /// Read face id attribute.
        err = facelement->QueryIntAttribute("FACEID", &faceid);
        ASSERTL0(err == TIXML_SUCCESS,
                 "Unable to read curve attribute FACEID.");

        /// Read text face element description.
        std::string elementStr;
        TiXmlNode *elementChild = facelement->FirstChild();

        while (elementChild)
        {
            // Accumulate all non-comment element data
            if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
                elementStr += " ";
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(), "Unable to read curve description body.");

        /// Parse out the element components corresponding to type of
        /// element.
        if (face == "F")
        {
            std::string typeStr = facelement->Attribute("TYPE");
            ASSERTL0(!typeStr.empty(), "TYPE must be specified in "
                                       "points definition");
            LibUtilities::PointsType type;
            const std::string *begStr = LibUtilities::kPointsTypeStr;
            const std::string *endStr =
                LibUtilities::kPointsTypeStr + LibUtilities::SIZE_PointsType;
            const std::string *ptsStr = std::find(begStr, endStr, typeStr);

            ASSERTL0(ptsStr != endStr, "Invalid points type.");
            type = (LibUtilities::PointsType)(ptsStr - begStr);

            std::string numptsStr = facelement->Attribute("NUMPOINTS");
            ASSERTL0(!numptsStr.empty(),
                     "NUMPOINTS must be specified in points definition");
            int numPts = 0;
            std::stringstream s;
            s << numptsStr;
            s >> numPts;

            CurveSharedPtr curve(
                MemoryManager<Curve>::AllocateSharedPtr(faceid, type));

            ASSERTL0(numPts >= 3, "NUMPOINTS for face must be greater than 2");

            if (numPts == 3)
            {
                ASSERTL0(ptsStr != endStr, "Invalid points type.");
            }

            // Read points (x, y, z)
            NekDouble xval, yval, zval;
            std::istringstream elementDataStrm(elementStr.c_str());
            try
            {
                while (!elementDataStrm.fail())
                {
                    elementDataStrm >> xval >> yval >> zval;

                    // Need to check it here because we
                    // may not be good after the read
                    // indicating that there was nothing
                    // to read.
                    if (!elementDataStrm.fail())
                    {
                        PointGeomSharedPtr vert(
                            MemoryManager<PointGeom>::AllocateSharedPtr(
                                m_meshDimension, faceindx, xval, yval, zval));
                        curve->m_points.push_back(vert);
                    }
                }
            }
            catch (...)
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unable to read curve data for FACE: ") +
                          elementStr)
                             .c_str());
            }
            m_curvedFaces[faceid] = curve;

            facelement = facelement->NextSiblingElement("F");
        }
    }
}

void MeshGraphXml::ReadDomain()
{
    TiXmlElement *domain = NULL;
    /// Look for data in DOMAIN block.
    domain = m_xmlGeom->FirstChildElement("DOMAIN");

    ASSERTL0(domain, "Unable to find DOMAIN tag in file.");

    /// Elements are of the form: "<D ID = "N"> ... </D>".
    /// Read the ID field first.
    TiXmlElement *multidomains = domain->FirstChildElement("D");

    if (multidomains)
    {
        int nextDomainNumber = 0;
        while (multidomains)
        {
            int indx;
            int err = multidomains->QueryIntAttribute("ID", &indx);
            ASSERTL0(err == TIXML_SUCCESS,
                     "Unable to read attribute ID in Domain.");

            TiXmlNode *elementChild = multidomains->FirstChild();
            while (elementChild &&
                   elementChild->Type() != TiXmlNode::TINYXML_TEXT)
            {
                elementChild = elementChild->NextSibling();
            }

            ASSERTL0(elementChild, "Unable to read DOMAIN body.");
            std::string elementStr = elementChild->ToText()->ValueStr();

            elementStr = elementStr.substr(elementStr.find_first_not_of(" "));

            std::string::size_type indxBeg = elementStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = elementStr.find_last_of(']') - 1;
            std::string indxStr =
                elementStr.substr(indxBeg, indxEnd - indxBeg + 1);

            ASSERTL0(
                !indxStr.empty(),
                "Unable to read domain's composite index (index missing?).");

            // Read the domain composites.
            // Parse the composites into a list.
            map<int, CompositeSharedPtr> unrollDomain;
            GetCompositeList(indxStr, unrollDomain);
            m_domain.push_back(unrollDomain);

            ASSERTL0(!m_domain[nextDomainNumber++].empty(),
                     (std::string(
                          "Unable to obtain domain's referenced composite: ") +
                      indxStr)
                         .c_str());

            /// Keep looking
            multidomains = multidomains->NextSiblingElement("D");
        }
    }
    else // previous definition of just one composite
    {

        // find the non comment portion of the body.
        TiXmlNode *elementChild = domain->FirstChild();
        while (elementChild && elementChild->Type() != TiXmlNode::TINYXML_TEXT)
        {
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(elementChild, "Unable to read DOMAIN body.");
        std::string elementStr = elementChild->ToText()->ValueStr();

        elementStr = elementStr.substr(elementStr.find_first_not_of(" "));

        std::string::size_type indxBeg = elementStr.find_first_of('[') + 1;
        std::string::size_type indxEnd = elementStr.find_last_of(']') - 1;
        std::string indxStr = elementStr.substr(indxBeg, indxEnd - indxBeg + 1);

        ASSERTL0(!indxStr.empty(),
                 "Unable to read domain's composite index (index missing?).");

        // Read the domain composites.
        // Parse the composites into a list.
        map<int, CompositeSharedPtr> fullDomain;
        GetCompositeList(indxStr, fullDomain);
        m_domain.push_back(fullDomain);

        ASSERTL0(
            !m_domain[0].empty(),
            (std::string("Unable to obtain domain's referenced composite: ") +
             indxStr)
                .c_str());
    }
}

void MeshGraphXml::ReadEdges()
{
    CurveMap::iterator it;

    /// Look for elements in ELEMENT block.
    TiXmlElement *field = m_xmlGeom->FirstChildElement("EDGE");

    ASSERTL0(field, "Unable to find EDGE tag in file.");

    /// All elements are of the form: "<E ID="#"> ... </E>", with
    /// ? being the element type.
    /// Read the ID field first.
    TiXmlElement *edge = field->FirstChildElement("E");

    /// Since all edge data is one big text block, we need to accumulate
    /// all TINYXML_TEXT data and then parse it.  This approach effectively
    /// skips
    /// all comments or other node types since we only care about the
    /// edge list.  We cannot handle missing edge numbers as we could
    /// with missing element numbers due to the text block format.
    std::string edgeStr;
    int indx;

    while (edge)
    {
        int err = edge->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read edge attribute ID.");

        TiXmlNode *child = edge->FirstChild();
        edgeStr.clear();
        if (child->Type() == TiXmlNode::TINYXML_TEXT)
        {
            edgeStr += child->ToText()->ValueStr();
        }

        /// Now parse out the edges, three fields at a time.
        int vertex1, vertex2;
        std::istringstream edgeDataStrm(edgeStr.c_str());

        try
        {
            while (!edgeDataStrm.fail())
            {
                edgeDataStrm >> vertex1 >> vertex2;

                // Must check after the read because we
                // may be at the end and not know it.  If
                // we are at the end we will add a
                // duplicate of the last entry if we don't
                // check here.
                if (!edgeDataStrm.fail())
                {
                    PointGeomSharedPtr vertices[2] = {GetVertex(vertex1),
                                                      GetVertex(vertex2)};
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
        }
        catch (...)
        {
            NEKERROR(
                ErrorUtil::efatal,
                (std::string("Unable to read edge data: ") + edgeStr).c_str());
        }

        edge = edge->NextSiblingElement("E");
    }
}

void MeshGraphXml::ReadFaces()
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

        /// Read id attribute.
        int indx;
        int err = element->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read face attribute ID.");

        /// See if this face has curves.
        it = m_curvedFaces.find(indx);

        /// Read text element description.
        TiXmlNode *elementChild = element->FirstChild();
        std::string elementStr;
        while (elementChild)
        {
            if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(), "Unable to read face description body.");

        /// Parse out the element components corresponding to type of
        /// element.
        if (elementType == "T")
        {
            // Read three edge numbers
            int edge1, edge2, edge3;
            std::istringstream elementDataStrm(elementStr.c_str());

            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;

                ASSERTL0(
                    !elementDataStrm.fail(),
                    (std::string("Unable to read face data for TRIANGLE: ") +
                     elementStr)
                        .c_str());

                /// Create a TriGeom to hold the new definition.
                SegGeomSharedPtr edges[TriGeom::kNedges] = {
                    GetSegGeom(edge1), GetSegGeom(edge2), GetSegGeom(edge3)};

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
            catch (...)
            {
                NEKERROR(
                    ErrorUtil::efatal,
                    (std::string("Unable to read face data for TRIANGLE: ") +
                     elementStr)
                        .c_str());
            }
        }
        else if (elementType == "Q")
        {
            // Read four edge numbers
            int edge1, edge2, edge3, edge4;
            std::istringstream elementDataStrm(elementStr.c_str());

            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;
                elementDataStrm >> edge4;

                ASSERTL0(!elementDataStrm.fail(),
                         (std::string("Unable to read face data for QUAD: ") +
                          elementStr)
                             .c_str());

                /// Create a QuadGeom to hold the new definition.
                SegGeomSharedPtr edges[QuadGeom::kNedges] = {
                    GetSegGeom(edge1), GetSegGeom(edge2), GetSegGeom(edge3),
                    GetSegGeom(edge4)};

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
            catch (...)
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unable to read face data for QUAD: ") +
                          elementStr)
                             .c_str());
            }
        }
        // Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXml::ReadElements()
{
    switch (m_meshDimension)
    {
        case 1:
            ReadElements1D();
            break;
        case 2:
            ReadElements2D();
            break;
        case 3:
            ReadElements3D();
            break;
    }
}

void MeshGraphXml::ReadElements1D()
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
        int indx;
        int err = segment->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");

        TiXmlNode *elementChild = segment->FirstChild();
        while (elementChild && elementChild->Type() != TiXmlNode::TINYXML_TEXT)
        {
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(elementChild, "Unable to read element description body.");
        std::string elementStr = elementChild->ToText()->ValueStr();

        /// Parse out the element components corresponding to type of
        /// element.
        /// Read two vertex numbers
        int vertex1, vertex2;
        std::istringstream elementDataStrm(elementStr.c_str());

        try
        {
            elementDataStrm >> vertex1;
            elementDataStrm >> vertex2;

            ASSERTL0(!elementDataStrm.fail(),
                     (std::string("Unable to read element data for SEGMENT: ") +
                      elementStr)
                         .c_str());

            PointGeomSharedPtr vertices[2] = {GetVertex(vertex1),
                                              GetVertex(vertex2)};
            SegGeomSharedPtr seg;
            it = m_curvedEdges.find(indx);

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
        catch (...)
        {
            NEKERROR(ErrorUtil::efatal,
                     (std::string("Unable to read element data for segment: ") +
                      elementStr)
                         .c_str());
        }
        /// Keep looking for additional segments
        segment = segment->NextSiblingElement("S");
    }
}

void MeshGraphXml::ReadElements2D()
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

        /// Read id attribute.
        int indx;
        int err = element->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");

        it = m_curvedFaces.find(indx);

        /// Read text element description.
        TiXmlNode *elementChild = element->FirstChild();
        std::string elementStr;
        while (elementChild)
        {
            if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(),
                 "Unable to read element description body.");

        /// Parse out the element components corresponding to type of
        /// element.
        if (elementType == "T")
        {
            // Read three edge numbers
            int edge1, edge2, edge3;
            std::istringstream elementDataStrm(elementStr.c_str());

            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;

                ASSERTL0(
                    !elementDataStrm.fail(),
                    (std::string("Unable to read element data for TRIANGLE: ") +
                     elementStr)
                        .c_str());

                /// Create a TriGeom to hold the new definition.
                SegGeomSharedPtr edges[TriGeom::kNedges] = {
                    GetSegGeom(edge1), GetSegGeom(edge2), GetSegGeom(edge3)};

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
            catch (...)
            {
                NEKERROR(
                    ErrorUtil::efatal,
                    (std::string("Unable to read element data for TRIANGLE: ") +
                     elementStr)
                        .c_str());
            }
        }
        else if (elementType == "Q")
        {
            // Read four edge numbers
            int edge1, edge2, edge3, edge4;
            std::istringstream elementDataStrm(elementStr.c_str());

            try
            {
                elementDataStrm >> edge1;
                elementDataStrm >> edge2;
                elementDataStrm >> edge3;
                elementDataStrm >> edge4;

                ASSERTL0(
                    !elementDataStrm.fail(),
                    (std::string("Unable to read element data for QUAD: ") +
                     elementStr)
                        .c_str());

                /// Create a QuadGeom to hold the new definition.
                SegGeomSharedPtr edges[QuadGeom::kNedges] = {
                    GetSegGeom(edge1), GetSegGeom(edge2), GetSegGeom(edge3),
                    GetSegGeom(edge4)};

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
            catch (...)
            {
                NEKERROR(
                    ErrorUtil::efatal,
                    (std::string("Unable to read element data for QUAD: ") +
                     elementStr)
                        .c_str());
            }
        }
        /// Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXml::ReadElements3D()
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

        /// Read id attribute.
        int indx;
        int err = element->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read element attribute ID.");

        /// Read text element description.
        TiXmlNode *elementChild = element->FirstChild();
        std::string elementStr;
        while (elementChild)
        {
            if (elementChild->Type() == TiXmlNode::TINYXML_TEXT)
            {
                elementStr += elementChild->ToText()->ValueStr();
            }
            elementChild = elementChild->NextSibling();
        }

        ASSERTL0(!elementStr.empty(),
                 "Unable to read element description body.");

        std::istringstream elementDataStrm(elementStr.c_str());

        /// Parse out the element components corresponding to type of
        /// element.

        // Tetrahedral
        if (elementType == "A")
        {
            try
            {
                /// Create arrays for the tri and quad faces.
                const int kNfaces  = TetGeom::kNfaces;
                const int kNtfaces = TetGeom::kNtfaces;
                const int kNqfaces = TetGeom::kNqfaces;
                TriGeomSharedPtr tfaces[kNtfaces];
                int Ntfaces = 0;
                int Nqfaces = 0;

                /// Fill the arrays and make sure there aren't too many
                /// faces.
                std::stringstream errorstring;
                errorstring << "Element " << indx << " must have " << kNtfaces
                            << " triangle face(s), and " << kNqfaces
                            << " quadrilateral face(s).";
                for (int i = 0; i < kNfaces; i++)
                {
                    int faceID;
                    elementDataStrm >> faceID;
                    Geometry2DSharedPtr face = GetGeometry2D(faceID);
                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << faceID;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
                        tfaces[Ntfaces++] =
                            static_pointer_cast<TriGeom>(face);
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
                    }
                }

                /// Make sure all of the face indicies could be read, and
                /// that there weren't too few.
                ASSERTL0(!elementDataStrm.fail(),
                         (std::string(
                              "Unable to read element data for TETRAHEDRON: ") +
                          elementStr)
                             .c_str());
                ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
                ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                TetGeomSharedPtr tetgeom(
                    MemoryManager<TetGeom>::AllocateSharedPtr(indx, tfaces));

                m_tetGeoms[indx] = tetgeom;
                PopulateFaceToElMap(tetgeom, kNfaces);
            }
            catch (...)
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string(
                              "Unable to read element data for TETRAHEDRON: ") +
                          elementStr)
                             .c_str());
            }
        }
        // Pyramid
        else if (elementType == "P")
        {
            try
            {
                /// Create arrays for the tri and quad faces.
                const int kNfaces  = PyrGeom::kNfaces;
                const int kNtfaces = PyrGeom::kNtfaces;
                const int kNqfaces = PyrGeom::kNqfaces;
                Geometry2DSharedPtr faces[kNfaces];
                int Nfaces  = 0;
                int Ntfaces = 0;
                int Nqfaces = 0;

                /// Fill the arrays and make sure there aren't too many
                /// faces.
                std::stringstream errorstring;
                errorstring << "Element " << indx << " must have " << kNtfaces
                            << " triangle face(s), and " << kNqfaces
                            << " quadrilateral face(s).";
                for (int i = 0; i < kNfaces; i++)
                {
                    int faceID;
                    elementDataStrm >> faceID;
                    Geometry2DSharedPtr face = GetGeometry2D(faceID);
                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << faceID;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
                        faces[Nfaces++] =
                            static_pointer_cast<TriGeom>(face);
                        Ntfaces++;
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
                        faces[Nfaces++] =
                            static_pointer_cast<QuadGeom>(face);
                        Nqfaces++;
                    }
                }

                /// Make sure all of the face indicies could be read, and
                /// that there weren't too few.
                ASSERTL0(
                    !elementDataStrm.fail(),
                    (std::string("Unable to read element data for PYRAMID: ") +
                     elementStr)
                        .c_str());
                ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
                ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                PyrGeomSharedPtr pyrgeom(
                    MemoryManager<PyrGeom>::AllocateSharedPtr(indx, faces));

                m_pyrGeoms[indx] = pyrgeom;
                PopulateFaceToElMap(pyrgeom, kNfaces);
            }
            catch (...)
            {
                NEKERROR(
                    ErrorUtil::efatal,
                    (std::string("Unable to read element data for PYRAMID: ") +
                     elementStr)
                        .c_str());
            }
        }
        // Prism
        else if (elementType == "R")
        {
            try
            {
                /// Create arrays for the tri and quad faces.
                const int kNfaces  = PrismGeom::kNfaces;
                const int kNtfaces = PrismGeom::kNtfaces;
                const int kNqfaces = PrismGeom::kNqfaces;
                Geometry2DSharedPtr faces[kNfaces];
                int Ntfaces = 0;
                int Nqfaces = 0;
                int Nfaces  = 0;

                /// Fill the arrays and make sure there aren't too many
                /// faces.
                std::stringstream errorstring;
                errorstring << "Element " << indx << " must have " << kNtfaces
                            << " triangle face(s), and " << kNqfaces
                            << " quadrilateral face(s).";

                for (int i = 0; i < kNfaces; i++)
                {
                    int faceID;
                    elementDataStrm >> faceID;
                    Geometry2DSharedPtr face = GetGeometry2D(faceID);
                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << faceID;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
                        faces[Nfaces++] =
                            std::static_pointer_cast<TriGeom>(face);
                        Ntfaces++;
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
                        faces[Nfaces++] =
                            std::static_pointer_cast<QuadGeom>(face);
                        Nqfaces++;
                    }
                }

                /// Make sure all of the face indicies could be read, and
                /// that there weren't too few.
                ASSERTL0(
                    !elementDataStrm.fail(),
                    (std::string("Unable to read element data for PRISM: ") +
                     elementStr)
                        .c_str());
                ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
                ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                PrismGeomSharedPtr prismgeom(
                    MemoryManager<PrismGeom>::AllocateSharedPtr(indx, faces));

                m_prismGeoms[indx] = prismgeom;
                PopulateFaceToElMap(prismgeom, kNfaces);
            }
            catch (...)
            {
                NEKERROR(
                    ErrorUtil::efatal,
                    (std::string("Unable to read element data for PRISM: ") +
                     elementStr)
                        .c_str());
            }
        }
        // Hexahedral
        else if (elementType == "H")
        {
            try
            {
                /// Create arrays for the tri and quad faces.
                const int kNfaces  = HexGeom::kNfaces;
                const int kNtfaces = HexGeom::kNtfaces;
                const int kNqfaces = HexGeom::kNqfaces;
                // TriGeomSharedPtr tfaces[kNtfaces];
                QuadGeomSharedPtr qfaces[kNqfaces];
                int Ntfaces = 0;
                int Nqfaces = 0;

                /// Fill the arrays and make sure there aren't too many
                /// faces.
                std::stringstream errorstring;
                errorstring << "Element " << indx << " must have " << kNtfaces
                            << " triangle face(s), and " << kNqfaces
                            << " quadrilateral face(s).";
                for (int i = 0; i < kNfaces; i++)
                {
                    int faceID;
                    elementDataStrm >> faceID;
                    Geometry2DSharedPtr face = GetGeometry2D(faceID);
                    if (face == Geometry2DSharedPtr() ||
                        (face->GetShapeType() != LibUtilities::eTriangle &&
                         face->GetShapeType() != LibUtilities::eQuadrilateral))
                    {
                        std::stringstream errorstring;
                        errorstring << "Element " << indx
                                    << " has invalid face: " << faceID;
                        ASSERTL0(false, errorstring.str().c_str());
                    }
                    else if (face->GetShapeType() == LibUtilities::eTriangle)
                    {
                        ASSERTL0(Ntfaces < kNtfaces, errorstring.str().c_str());
                        // tfaces[Ntfaces++] =
                        // boost::static_pointer_cast<TriGeom>(face);
                    }
                    else if (face->GetShapeType() ==
                             LibUtilities::eQuadrilateral)
                    {
                        ASSERTL0(Nqfaces < kNqfaces, errorstring.str().c_str());
                        qfaces[Nqfaces++] =
                            std::static_pointer_cast<QuadGeom>(face);
                    }
                }

                /// Make sure all of the face indicies could be read, and
                /// that there weren't too few.
                ASSERTL0(!elementDataStrm.fail(),
                         (std::string(
                              "Unable to read element data for HEXAHEDRAL: ") +
                          elementStr)
                             .c_str());
                ASSERTL0(Ntfaces == kNtfaces, errorstring.str().c_str());
                ASSERTL0(Nqfaces == kNqfaces, errorstring.str().c_str());

                HexGeomSharedPtr hexgeom(
                    MemoryManager<HexGeom>::AllocateSharedPtr(indx, qfaces));

                m_hexGeoms[indx] = hexgeom;
                PopulateFaceToElMap(hexgeom, kNfaces);
            }
            catch (...)
            {
                NEKERROR(ErrorUtil::efatal,
                         (std::string(
                              "Unable to read element data for HEXAHEDRAL: ") +
                          elementStr)
                             .c_str());
            }
        }
        /// Keep looking
        element = element->NextSiblingElement();
    }
}

void MeshGraphXml::ReadComposites()
{
    TiXmlElement *field = NULL;

    /// Look for elements in ELEMENT block.
    field = m_xmlGeom->FirstChildElement("COMPOSITE");

    ASSERTL0(field, "Unable to find COMPOSITE tag in file.");

    TiXmlElement *node = field->FirstChildElement("C");

    // Sequential counter for the composite numbers.
    int nextCompositeNumber = -1;

    while (node)
    {
        /// All elements are of the form: "<? ID="#"> ... </?>", with
        /// ? being the element type.

        nextCompositeNumber++;

        int indx;
        int err = node->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
        // ASSERTL0(indx == nextCompositeNumber, "Composite IDs must begin with
        // zero and be sequential.");

        TiXmlNode *compositeChild = node->FirstChild();
        // This is primarily to skip comments that may be present.
        // Comments appear as nodes just like elements.
        // We are specifically looking for text in the body
        // of the definition.
        while (compositeChild &&
               compositeChild->Type() != TiXmlNode::TINYXML_TEXT)
        {
            compositeChild = compositeChild->NextSibling();
        }

        ASSERTL0(compositeChild, "Unable to read composite definition body.");
        std::string compositeStr = compositeChild->ToText()->ValueStr();

        /// Parse out the element components corresponding to type of element.

        std::istringstream compositeDataStrm(compositeStr.c_str());

        try
        {
            bool first = true;
            std::string prevCompositeElementStr;

            while (!compositeDataStrm.fail())
            {
                std::string compositeElementStr;
                compositeDataStrm >> compositeElementStr;

                if (!compositeDataStrm.fail())
                {
                    if (first)
                    {
                        first = false;

                        CompositeSharedPtr curVector =
                            MemoryManager<Composite>::AllocateSharedPtr();
                        m_meshComposites[indx] = curVector;
                    }

                    if (compositeElementStr.length() > 0)
                    {
                        ResolveGeomRef(prevCompositeElementStr,
                                       compositeElementStr,
                                       m_meshComposites[indx]);
                    }
                    prevCompositeElementStr = compositeElementStr;
                }
            }
        }
        catch (...)
        {
            NEKERROR(
                ErrorUtil::efatal,
                (std::string("Unable to read COMPOSITE data for composite: ") +
                 compositeStr)
                    .c_str());
        }

        /// Keep looking for additional composite definitions.
        node = node->NextSiblingElement("C");
    }

    ASSERTL0(nextCompositeNumber >= 0,
             "At least one composite must be specified.");
}

void MeshGraphXml::ResolveGeomRef(const std::string &prevToken,
                                  const std::string &token,
                                  CompositeSharedPtr &composite)
{
    switch (m_meshDimension)
    {
        case 1:
            ResolveGeomRef1D(prevToken, token, composite);
            break;
        case 2:
            ResolveGeomRef2D(prevToken, token, composite);
            break;
        case 3:
            ResolveGeomRef3D(prevToken, token, composite);
            break;
    }
}

void MeshGraphXml::ResolveGeomRef1D(const std::string &prevToken,
                                    const std::string &token,
                                    CompositeSharedPtr &composite)
{
    try
    {
        std::istringstream tokenStream(token);
        std::istringstream prevTokenStream(prevToken);

        char type;
        char prevType;

        tokenStream >> type;

        std::string::size_type indxBeg = token.find_first_of('[') + 1;
        std::string::size_type indxEnd = token.find_last_of(']') - 1;

        ASSERTL0(
            indxBeg <= indxEnd,
            (std::string("Error reading index definition:") + token).c_str());

        std::string indxStr = token.substr(indxBeg, indxEnd - indxBeg + 1);

        typedef vector<unsigned int> SeqVectorType;
        SeqVectorType seqVector;

        if (!ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector))
        {
            NEKERROR(ErrorUtil::efatal,
                     (std::string("Ill-formed sequence definition: ") + indxStr)
                         .c_str());
        }

        prevTokenStream >> prevType;

        // All composites must be of the same dimension.
        bool validSequence =
            (prevToken.empty() || // No previous, then current is just fine.
             (type == 'V' && prevType == 'V') ||
             (type == 'S' && prevType == 'S'));

        ASSERTL0(validSequence,
                 (std::string("Invalid combination of composite items: ") +
                  type + " and " + prevType + ".")
                     .c_str());

        switch (type)
        {
            case 'V': // Vertex
                for (SeqVectorType::iterator iter = seqVector.begin();
                     iter != seqVector.end(); ++iter)
                {
                    if (m_vertSet.find(*iter) == m_vertSet.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *iter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown vertex index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_vertSet[*iter]);
                    }
                }
                break;

            case 'S': // Segment
                for (SeqVectorType::iterator iter = seqVector.begin();
                     iter != seqVector.end(); ++iter)
                {
                    if (m_segGeoms.find(*iter) == m_segGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *iter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown segment index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_segGeoms[*iter]);
                    }
                }
                break;

            default:
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unrecognized composite token: ") + token)
                             .c_str());
        }
    }
    catch (...)
    {
        NEKERROR(ErrorUtil::efatal,
                 (std::string("Problem processing composite token: ") + token)
                     .c_str());
    }

    return;
}

void MeshGraphXml::ResolveGeomRef2D(const std::string &prevToken,
                                    const std::string &token,
                                    CompositeSharedPtr &composite)
{
    try
    {
        std::istringstream tokenStream(token);
        std::istringstream prevTokenStream(prevToken);

        char type;
        char prevType;

        tokenStream >> type;

        std::string::size_type indxBeg = token.find_first_of('[') + 1;
        std::string::size_type indxEnd = token.find_last_of(']') - 1;

        ASSERTL0(
            indxBeg <= indxEnd,
            (std::string("Error reading index definition:") + token).c_str());

        std::string indxStr = token.substr(indxBeg, indxEnd - indxBeg + 1);
        std::vector<unsigned int> seqVector;
        std::vector<unsigned int>::iterator seqIter;

        bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);

        ASSERTL0(err,
                 (std::string("Error reading composite elements: ") + indxStr)
                     .c_str());

        prevTokenStream >> prevType;

        // All composites must be of the same dimension.
        bool validSequence =
            (prevToken.empty() || // No previous, then current is just fine.
             (type == 'V' && prevType == 'V') ||
             (type == 'E' && prevType == 'E') ||
             ((type == 'T' || type == 'Q') &&
              (prevType == 'T' || prevType == 'Q')));

        ASSERTL0(validSequence,
                 (std::string("Invalid combination of composite items: ") +
                  type + " and " + prevType + ".")
                     .c_str());

        switch (type)
        {
            case 'E': // Edge
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_segGeoms.find(*seqIter) == m_segGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown edge index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_segGeoms[*seqIter]);
                    }
                }
                break;

            case 'T': // Triangle
            {
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_triGeoms.count(*seqIter) == 0)
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown triangle index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_triGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_triGeoms[*seqIter]);
                        }
                    }
                }
            }
            break;

            case 'Q': // Quad
            {
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_quadGeoms.count(*seqIter) == 0)
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown quad index: ") + errStr +
                                  std::string(" in Composite section"))
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_quadGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_quadGeoms[*seqIter]);
                        }
                    }
                }
            }
            break;

            case 'V': // Vertex
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (*seqIter >= m_vertSet.size())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown vertex index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_vertSet[*seqIter]);
                    }
                }
                break;

            default:
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unrecognized composite token: ") + token)
                             .c_str());
        }
    }
    catch (...)
    {
        NEKERROR(ErrorUtil::efatal,
                 (std::string("Problem processing composite token: ") + token)
                     .c_str());
    }

    return;
}

void MeshGraphXml::ResolveGeomRef3D(const std::string &prevToken,
                                    const std::string &token,
                                    CompositeSharedPtr &composite)
{
    try
    {
        std::istringstream tokenStream(token);
        std::istringstream prevTokenStream(prevToken);

        char type;
        char prevType;

        tokenStream >> type;

        std::string::size_type indxBeg = token.find_first_of('[') + 1;
        std::string::size_type indxEnd = token.find_last_of(']') - 1;

        ASSERTL0(
            indxBeg <= indxEnd,
            (std::string("Error reading index definition:") + token).c_str());

        std::string indxStr = token.substr(indxBeg, indxEnd - indxBeg + 1);

        std::vector<unsigned int> seqVector;
        std::vector<unsigned int>::iterator seqIter;

        bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);

        ASSERTL0(err,
                 (std::string("Error reading composite elements: ") + indxStr)
                     .c_str());

        prevTokenStream >> prevType;

        // All composites must be of the same dimension.  This map makes things
        // clean to compare.
        map<char, int> typeMap;
        typeMap['V'] = 1; // Vertex
        typeMap['E'] = 1; // Edge
        typeMap['T'] = 2; // Triangle
        typeMap['Q'] = 2; // Quad
        typeMap['A'] = 3; // Tet
        typeMap['P'] = 3; // Pyramid
        typeMap['R'] = 3; // Prism
        typeMap['H'] = 3; // Hex

        // Make sure only geoms of the same dimension are combined.
        bool validSequence =
            (prevToken.empty() || (typeMap[type] == typeMap[prevType]));

        ASSERTL0(validSequence,
                 (std::string("Invalid combination of composite items: ") +
                  type + " and " + prevType + ".")
                     .c_str());

        switch (type)
        {
            case 'V': // Vertex
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_vertSet.find(*seqIter) == m_vertSet.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown vertex index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_vertSet[*seqIter]);
                    }
                }
                break;

            case 'E': // Edge
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_segGeoms.find(*seqIter) == m_segGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown edge index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        composite->m_geomVec.push_back(m_segGeoms[*seqIter]);
                    }
                }
                break;

            case 'F': // Face
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    Geometry2DSharedPtr face = GetGeometry2D(*seqIter);
                    if (face == Geometry2DSharedPtr())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown face index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*face))
                        {
                            composite->m_geomVec.push_back(face);
                        }
                    }
                }
                break;

            case 'T': // Triangle
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_triGeoms.find(*seqIter) == m_triGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown triangle index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_triGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_triGeoms[*seqIter]);
                        }
                    }
                }
                break;

            case 'Q': // Quad
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_quadGeoms.find(*seqIter) == m_quadGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown quad index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_quadGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_quadGeoms[*seqIter]);
                        }
                    }
                }
                break;

            // Tetrahedron
            case 'A':
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_tetGeoms.find(*seqIter) == m_tetGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown tet index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_tetGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_tetGeoms[*seqIter]);
                        }
                    }
                }
                break;

            // Pyramid
            case 'P':
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_pyrGeoms.find(*seqIter) == m_pyrGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(
                            ErrorUtil::ewarning,
                            (std::string("Unknown pyramid index: ") + errStr)
                                .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_pyrGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_pyrGeoms[*seqIter]);
                        }
                    }
                }
                break;

            // Prism
            case 'R':
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_prismGeoms.find(*seqIter) == m_prismGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown prism index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_prismGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_prismGeoms[*seqIter]);
                        }
                    }
                }
                break;

            // Hex
            case 'H':
                for (seqIter = seqVector.begin(); seqIter != seqVector.end();
                     ++seqIter)
                {
                    if (m_hexGeoms.find(*seqIter) == m_hexGeoms.end())
                    {
                        char errStr[16] = "";
                        ::sprintf(errStr, "%d", *seqIter);
                        NEKERROR(ErrorUtil::ewarning,
                                 (std::string("Unknown hex index: ") + errStr)
                                     .c_str());
                    }
                    else
                    {
                        if (CheckRange(*m_hexGeoms[*seqIter]))
                        {
                            composite->m_geomVec.push_back(
                                m_hexGeoms[*seqIter]);
                        }
                    }
                }
                break;

            default:
                NEKERROR(ErrorUtil::efatal,
                         (std::string("Unrecognized composite token: ") + token)
                             .c_str());
        }
    }
    catch (...)
    {
        NEKERROR(ErrorUtil::efatal,
                 (std::string("Problem processing composite token: ") + token)
                     .c_str());
    }

    return;
}

void MeshGraphXml::WriteVertices(TiXmlElement *geomTag, PointGeomMap &verts)
{
    TiXmlElement *vertTag = new TiXmlElement("VERTEX");

    for (auto &i : verts)
    {
        stringstream s;
        s << scientific << setprecision(8) << (*i.second)(0) << " "
          << (*i.second)(1) << " " << (*i.second)(2);
        TiXmlElement *v = new TiXmlElement("V");
        v->SetAttribute("ID", i.second->GetGlobalID());
        v->LinkEndChild(new TiXmlText(s.str()));
        vertTag->LinkEndChild(v);
    }

    geomTag->LinkEndChild(vertTag);
}

void MeshGraphXml::WriteEdges(TiXmlElement *geomTag, SegGeomMap &edges)
{
    TiXmlElement *edgeTag =
        new TiXmlElement(m_meshDimension == 1 ? "ELEMENT" : "EDGE");
    string tag = m_meshDimension == 1 ? "S" : "E";

    for (auto &i : edges)
    {
        stringstream s;
        SegGeomSharedPtr seg = i.second;
        s << seg->GetVid(0) << " " << seg->GetVid(1);
        TiXmlElement *e = new TiXmlElement(tag);
        e->SetAttribute("ID", i.first);
        e->LinkEndChild(new TiXmlText(s.str()));
        edgeTag->LinkEndChild(e);
    }

    geomTag->LinkEndChild(edgeTag);
}

void MeshGraphXml::WriteTris(TiXmlElement *faceTag, TriGeomMap &tris)
{
    string tag = "T";

    for (auto &i : tris)
    {
        stringstream s;
        TriGeomSharedPtr tri = i.second;
        s << tri->GetEid(0) << " " << tri->GetEid(1) << " " << tri->GetEid(2);
        TiXmlElement *t = new TiXmlElement(tag);
        t->SetAttribute("ID", i.first);
        t->LinkEndChild(new TiXmlText(s.str()));
        faceTag->LinkEndChild(t);
    }
}

void MeshGraphXml::WriteQuads(TiXmlElement *faceTag, QuadGeomMap &quads)
{
    string tag = "Q";

    for (auto &i : quads)
    {
        stringstream s;
        QuadGeomSharedPtr quad = i.second;
        s << quad->GetEid(0) << " " << quad->GetEid(1) << " " << quad->GetEid(2)
          << " " << quad->GetEid(3);
        TiXmlElement *q = new TiXmlElement(tag);
        q->SetAttribute("ID", i.first);
        q->LinkEndChild(new TiXmlText(s.str()));
        faceTag->LinkEndChild(q);
    }
}

void MeshGraphXml::WriteHexs(TiXmlElement *elmtTag, HexGeomMap &hexs)
{
    string tag = "H";

    for (auto &i : hexs)
    {
        stringstream s;
        HexGeomSharedPtr hex = i.second;
        s << hex->GetFid(0) << " " << hex->GetFid(1) << " " << hex->GetFid(2)
          << " " << hex->GetFid(3) << " " << hex->GetFid(4) << " "
          << hex->GetFid(5) << " ";
        TiXmlElement *h = new TiXmlElement(tag);
        h->SetAttribute("ID", i.first);
        h->LinkEndChild(new TiXmlText(s.str()));
        elmtTag->LinkEndChild(h);
    }
}

void MeshGraphXml::WritePrisms(TiXmlElement *elmtTag, PrismGeomMap &pris)
{
    string tag = "R";

    for (auto &i : pris)
    {
        stringstream s;
        PrismGeomSharedPtr prism = i.second;
        s << prism->GetFid(0) << " " << prism->GetFid(1) << " "
          << prism->GetFid(2) << " " << prism->GetFid(3) << " "
          << prism->GetFid(4) << " ";
        TiXmlElement *p = new TiXmlElement(tag);
        p->SetAttribute("ID", i.first);
        p->LinkEndChild(new TiXmlText(s.str()));
        elmtTag->LinkEndChild(p);
    }
}

void MeshGraphXml::WritePyrs(TiXmlElement *elmtTag, PyrGeomMap &pyrs)
{
    string tag = "P";

    for (auto &i : pyrs)
    {
        stringstream s;
        PyrGeomSharedPtr pyr = i.second;
        s << pyr->GetFid(0) << " " << pyr->GetFid(1) << " " << pyr->GetFid(2)
          << " " << pyr->GetFid(3) << " " << pyr->GetFid(4) << " ";
        TiXmlElement *p = new TiXmlElement(tag);
        p->SetAttribute("ID", i.first);
        p->LinkEndChild(new TiXmlText(s.str()));
        elmtTag->LinkEndChild(p);
    }
}

void MeshGraphXml::WriteTets(TiXmlElement *elmtTag, TetGeomMap &tets)
{
    string tag = "A";

    for (auto &i : tets)
    {
        stringstream s;
        TetGeomSharedPtr tet = i.second;
        s << tet->GetFid(0) << " " << tet->GetFid(1) << " " << tet->GetFid(2)
          << " " << tet->GetFid(3) << " ";
        TiXmlElement *t = new TiXmlElement(tag);
        t->SetAttribute("ID", i.first);
        t->LinkEndChild(new TiXmlText(s.str()));
        elmtTag->LinkEndChild(t);
    }
}

void MeshGraphXml::WriteCurves(TiXmlElement *geomTag, CurveMap &edges,
                               CurveMap &faces)
{
    TiXmlElement *curveTag = new TiXmlElement("CURVED");
    CurveMap::iterator curveIt;
    int curveId = 0;

    for (curveIt = edges.begin(); curveIt != edges.end(); ++curveIt)
    {
        CurveSharedPtr curve = curveIt->second;
        TiXmlElement *c      = new TiXmlElement("E");
        stringstream s;
        s.precision(8);

        for (int j = 0; j < curve->m_points.size(); ++j)
        {
            SpatialDomains::PointGeomSharedPtr p = curve->m_points[j];
            s << scientific << (*p)(0) << " " << (*p)(1) << " " << (*p)(2)
              << "   ";
        }

        c->SetAttribute("ID", curveId++);
        c->SetAttribute("EDGEID", curve->m_curveID);
        c->SetAttribute("NUMPOINTS", curve->m_points.size());
        c->SetAttribute("TYPE", LibUtilities::kPointsTypeStr[curve->m_ptype]);
        c->LinkEndChild(new TiXmlText(s.str()));
        curveTag->LinkEndChild(c);
    }

    for (curveIt = faces.begin(); curveIt != faces.end(); ++curveIt)
    {
        CurveSharedPtr curve = curveIt->second;
        TiXmlElement *c      = new TiXmlElement("F");
        stringstream s;
        s.precision(8);

        for (int j = 0; j < curve->m_points.size(); ++j)
        {
            SpatialDomains::PointGeomSharedPtr p = curve->m_points[j];
            s << scientific << (*p)(0) << " " << (*p)(1) << " " << (*p)(2)
              << "   ";
        }

        c->SetAttribute("ID", curveId++);
        c->SetAttribute("FACEID", curve->m_curveID);
        c->SetAttribute("NUMPOINTS", curve->m_points.size());
        c->SetAttribute("TYPE", LibUtilities::kPointsTypeStr[curve->m_ptype]);
        c->LinkEndChild(new TiXmlText(s.str()));
        curveTag->LinkEndChild(c);
    }

    geomTag->LinkEndChild(curveTag);
}

void MeshGraphXml::WriteComposites(TiXmlElement *geomTag, CompositeMap &comps)
{
    TiXmlElement *compTag = new TiXmlElement("COMPOSITE");

    for (auto &cIt : comps)
    {
        if (cIt.second->m_geomVec.size() == 0)
        {
            continue;
        }

        TiXmlElement *c = new TiXmlElement("C");
        c->SetAttribute("ID", cIt.first);
        c->LinkEndChild(new TiXmlText(GetCompositeString(cIt.second)));
        compTag->LinkEndChild(c);
    }

    geomTag->LinkEndChild(compTag);
}

void MeshGraphXml::WriteDomain(TiXmlElement *geomTag,
                               vector<CompositeMap> &domain)
{
    TiXmlElement *domTag = new TiXmlElement("DOMAIN");
    stringstream domString;

    // @todo Fix this to accomodate multi domain output
    vector<unsigned int> idxList;
    for (auto cIt = domain[0].begin(); cIt != domain[0].end(); ++cIt)
    {
        idxList.push_back(cIt->first);
    }

    domString << " C[" << ParseUtils::GenerateSeqString(idxList) << "] ";
    domTag->LinkEndChild(new TiXmlText(domString.str()));
    geomTag->LinkEndChild(domTag);
}

void MeshGraphXml::WriteDefaultExpansion(TiXmlElement *root)
{
    TiXmlElement *expTag = new TiXmlElement("EXPANSIONS");

    for (auto it = m_meshComposites.begin(); it != m_meshComposites.end(); it++)
    {
        if (it->second->m_geomVec[0]->GetShapeDim() == m_meshDimension)
        {
            TiXmlElement *exp = new TiXmlElement("E");
            exp->SetAttribute("COMPOSITE",
                              "C[" + boost::lexical_cast<string>(it->first) +
                                  "]");
            exp->SetAttribute("NUMMODES", 4);
            exp->SetAttribute("TYPE", "MODIFIED");
            exp->SetAttribute("FIELDS", "u");

            expTag->LinkEndChild(exp);
        }
    }
    root->LinkEndChild(expTag);
}

/**
 * @brief Write out an XML file containing the GEOMETRY block
 * representing this MeshGraph instance inside a NEKTAR tag.
 */
void MeshGraphXml::WriteGeometry(
    std::string                          &outfilename,
    bool                                  defaultExp,
    const LibUtilities::FieldMetaDataMap &metadata)
{
    // Create empty TinyXML document.
    TiXmlDocument doc;
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);
    TiXmlElement *geomTag = new TiXmlElement("GEOMETRY");
    root->LinkEndChild(geomTag);

    // Add provenance information using FieldIO library.
    LibUtilities::FieldIO::AddInfoTag(
        LibUtilities::XmlTagWriterSharedPtr(
            new LibUtilities::XmlTagWriter(root)),
        metadata);

    // Update attributes with dimensions.
    geomTag->SetAttribute("DIM", m_meshDimension);
    geomTag->SetAttribute("SPACE", m_spaceDimension);

    // Clear existing elements.
    geomTag->Clear();

    // Write out informatio
    WriteVertices(geomTag, m_vertSet);
    WriteEdges(geomTag, m_segGeoms);
    if (m_meshDimension > 1)
    {
        TiXmlElement *faceTag =
            new TiXmlElement(m_meshDimension == 2 ? "ELEMENT" : "FACE");

        WriteTris(faceTag, m_triGeoms);
        WriteQuads(faceTag, m_quadGeoms);
        geomTag->LinkEndChild(faceTag);
    }
    if (m_meshDimension > 2)
    {
        TiXmlElement *elmtTag = new TiXmlElement("ELEMENT");

        WriteHexs(elmtTag, m_hexGeoms);
        WritePyrs(elmtTag, m_pyrGeoms);
        WritePrisms(elmtTag, m_prismGeoms);
        WriteTets(elmtTag, m_tetGeoms);

        geomTag->LinkEndChild(elmtTag);
    }
    WriteCurves(geomTag, m_curvedEdges, m_curvedFaces);
    WriteComposites(geomTag, m_meshComposites);
    WriteDomain(geomTag, m_domain);

    if (defaultExp)
    {
        WriteDefaultExpansion(root);
    }

    // Save file.
    doc.SaveFile(outfilename);
}

void MeshGraphXml::WriteXMLGeometry(std::string outname,
                                    vector<set<unsigned int>> elements,
                                    vector<unsigned int> partitions)
{
    // so in theory this function is used by the mesh partitioner
    // giving instructions on how to write out a paritioned mesh.
    // the theory goes that the elements stored in the meshgraph are the
    // "whole" mesh so based on the information from the elmements list
    // we can filter the mesh entities and write some individual files
    // hopefully

    // this is xml so we are going to write a directory with lots of
    // xml files
    string dirname = outname + "_xml";
    boost::filesystem::path pdirname(dirname);

    if (!boost::filesystem::is_directory(dirname))
    {
        boost::filesystem::create_directory(dirname);
    }

    ASSERTL0(elements.size() == partitions.size(),
             "error in partitioned information");

    for (int i = 0; i < partitions.size(); i++)
    {
        TiXmlDocument doc;
        TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
        doc.LinkEndChild(decl);

        TiXmlElement *root = doc.FirstChildElement("NEKTAR");
        TiXmlElement *geomTag;

        // Try to find existing NEKTAR tag.
        if (!root)
        {
            root = new TiXmlElement("NEKTAR");
            doc.LinkEndChild(root);

            geomTag = new TiXmlElement("GEOMETRY");
            root->LinkEndChild(geomTag);
        }
        else
        {
            // Try to find existing GEOMETRY tag.
            geomTag = root->FirstChildElement("GEOMETRY");

            if (!geomTag)
            {
                geomTag = new TiXmlElement("GEOMETRY");
                root->LinkEndChild(geomTag);
            }
        }

        geomTag->SetAttribute("DIM", m_meshDimension);
        geomTag->SetAttribute("SPACE", m_spaceDimension);
        geomTag->SetAttribute("PARTITION", partitions[i]);

        // Add Mesh //
        // Get the elements
        HexGeomMap localHex;
        PyrGeomMap localPyr;
        PrismGeomMap localPrism;
        TetGeomMap localTet;
        TriGeomMap localTri;
        QuadGeomMap localQuad;
        SegGeomMap localEdge;
        PointGeomMap localVert;
        CurveMap localCurveEdge;
        CurveMap localCurveFace;

        vector<set<unsigned int>> entityIds(4);
        entityIds[m_meshDimension] = elements[i];

        switch (m_meshDimension)
        {
            case 3:
            {
                for (auto &j : entityIds[3])
                {
                    GeometrySharedPtr g;
                    if (m_hexGeoms.count(j))
                    {
                        g           = m_hexGeoms[j];
                        localHex[j] = m_hexGeoms[j];
                    }
                    else if (m_pyrGeoms.count(j))
                    {
                        g           = m_pyrGeoms[j];
                        localPyr[j] = m_pyrGeoms[j];
                    }
                    else if (m_prismGeoms.count(j))
                    {
                        g             = m_prismGeoms[j];
                        localPrism[j] = m_prismGeoms[j];
                    }
                    else if (m_tetGeoms.count(j))
                    {
                        g           = m_tetGeoms[j];
                        localTet[j] = m_tetGeoms[j];
                    }
                    else
                    {
                        ASSERTL0(false, "element in partition not found");
                    }

                    for (int k = 0; k < g->GetNumFaces(); k++)
                    {
                        entityIds[2].insert(g->GetFid(k));
                    }
                    for (int k = 0; k < g->GetNumEdges(); k++)
                    {
                        entityIds[1].insert(g->GetEid(k));
                    }
                    for (int k = 0; k < g->GetNumVerts(); k++)
                    {
                        entityIds[0].insert(g->GetVid(k));
                    }
                }
            }
            break;
            case 2:
            {
                for (auto &j : entityIds[2])
                {
                    GeometrySharedPtr g;
                    if (m_triGeoms.count(j))
                    {
                        g           = m_triGeoms[j];
                        localTri[j] = m_triGeoms[j];
                    }
                    else if (m_quadGeoms.count(j))
                    {
                        g            = m_quadGeoms[j];
                        localQuad[j] = m_quadGeoms[j];
                    }
                    else
                    {
                        ASSERTL0(false, "element in partition not found");
                    }

                    for (int k = 0; k < g->GetNumEdges(); k++)
                    {
                        entityIds[1].insert(g->GetEid(k));
                    }
                    for (int k = 0; k < g->GetNumVerts(); k++)
                    {
                        entityIds[0].insert(g->GetVid(k));
                    }
                }
            }
            break;
            case 1:
            {
                for (auto &j : entityIds[1])
                {
                    GeometrySharedPtr g;
                    if (m_segGeoms.count(j))
                    {
                        g            = m_segGeoms[j];
                        localEdge[j] = m_segGeoms[j];
                    }
                    else
                    {
                        ASSERTL0(false, "element in partition not found");
                    }

                    for (int k = 0; k < g->GetNumVerts(); k++)
                    {
                        entityIds[0].insert(g->GetVid(k));
                    }
                }
            }
        }

        if (m_meshDimension > 2)
        {
            for (auto &j : entityIds[2])
            {
                if (m_triGeoms.count(j))
                {
                    localTri[j] = m_triGeoms[j];
                }
                else if (m_quadGeoms.count(j))
                {
                    localQuad[j] = m_quadGeoms[j];
                }
                else
                {
                    ASSERTL0(false, "face not found");
                }
            }
        }

        if (m_meshDimension > 1)
        {
            for (auto &j : entityIds[1])
            {
                if (m_segGeoms.count(j))
                {
                    localEdge[j] = m_segGeoms[j];
                }
                else
                {
                    ASSERTL0(false, "edge not found");
                }
            }
        }

        for (auto &j : entityIds[0])
        {
            if (m_vertSet.count(j))
            {
                localVert[j] = m_vertSet[j];
            }
            else
            {
                ASSERTL0(false, "vert not found");
            }
        }

        WriteVertices(geomTag, localVert);
        WriteEdges(geomTag, localEdge);
        if (m_meshDimension > 1)
        {
            TiXmlElement *faceTag =
                new TiXmlElement(m_meshDimension == 2 ? "ELEMENT" : "FACE");

            WriteTris(faceTag, localTri);
            WriteQuads(faceTag, localQuad);
            geomTag->LinkEndChild(faceTag);
        }
        if (m_meshDimension > 2)
        {
            TiXmlElement *elmtTag = new TiXmlElement("ELEMENT");

            WriteHexs(elmtTag, localHex);
            WritePyrs(elmtTag, localPyr);
            WritePrisms(elmtTag, localPrism);
            WriteTets(elmtTag, localTet);

            geomTag->LinkEndChild(elmtTag);
        }

        for (auto &j : localTri)
        {
            if (m_curvedFaces.count(j.first))
            {
                localCurveFace[j.first] = m_curvedFaces[j.first];
            }
        }
        for (auto &j : localQuad)
        {
            if (m_curvedFaces.count(j.first))
            {
                localCurveFace[j.first] = m_curvedFaces[j.first];
            }
        }
        for (auto &j : localEdge)
        {
            if (m_curvedEdges.count(j.first))
            {
                localCurveEdge[j.first] = m_curvedEdges[j.first];
            }
        }

        WriteCurves(geomTag, localCurveEdge, localCurveFace);

        CompositeMap localComp;

        for (auto &j : m_meshComposites)
        {
            CompositeSharedPtr comp = CompositeSharedPtr(new Composite);
            int dim                 = j.second->m_geomVec[0]->GetShapeDim();

            for (int k = 0; k < j.second->m_geomVec.size(); k++)
            {
                if (entityIds[dim].count(j.second->m_geomVec[k]->GetGlobalID()))
                {
                    comp->m_geomVec.push_back(j.second->m_geomVec[k]);
                }
            }

            if (comp->m_geomVec.size())
            {
                localComp[j.first] = comp;
            }
        }

        WriteComposites(geomTag, localComp);

        vector<CompositeMap> domain;
        CompositeMap domMap;
        for (auto &j : localComp)
        {
            if (j.second->m_geomVec[0]->GetShapeDim() == m_meshDimension)
            {
                domMap[j.first] = j.second;
            }
        }
        domain.push_back(domMap);

        WriteDomain(geomTag, domain);

        if (m_session->DefinesElement("NEKTAR/CONDITIONS"))
        {
            std::set<int> vBndRegionIdList;
            TiXmlElement *vConditions =
                new TiXmlElement(*m_session->GetElement("Nektar/Conditions"));
            TiXmlElement *vBndRegions =
                vConditions->FirstChildElement("BOUNDARYREGIONS");
            TiXmlElement *vBndConditions =
                vConditions->FirstChildElement("BOUNDARYCONDITIONS");
            TiXmlElement *vItem;

            if (vBndRegions)
            {
                TiXmlElement *vNewBndRegions =
                    new TiXmlElement("BOUNDARYREGIONS");
                vItem = vBndRegions->FirstChildElement();
                while (vItem)
                {
                    std::string vSeqStr =
                        vItem->FirstChild()->ToText()->Value();
                    std::string::size_type indxBeg =
                        vSeqStr.find_first_of('[') + 1;
                    std::string::size_type indxEnd =
                        vSeqStr.find_last_of(']') - 1;
                    vSeqStr = vSeqStr.substr(indxBeg, indxEnd - indxBeg + 1);
                    std::vector<unsigned int> vSeq;
                    ParseUtils::GenerateSeqVector(vSeqStr.c_str(), vSeq);

                    vector<unsigned int> idxList;

                    for (unsigned int i = 0; i < vSeq.size(); ++i)
                    {
                        if (localComp.find(vSeq[i]) != localComp.end())
                        {
                            idxList.push_back(vSeq[i]);
                        }
                    }
                    int p = atoi(vItem->Attribute("ID"));

                    std::string vListStr =
                        ParseUtils::GenerateSeqString(idxList);

                    if (vListStr.length() == 0)
                    {
                        TiXmlElement *tmp = vItem;
                        vItem             = vItem->NextSiblingElement();
                        vBndRegions->RemoveChild(tmp);
                    }
                    else
                    {
                        vListStr                  = "C[" + vListStr + "]";
                        TiXmlText *vList          = new TiXmlText(vListStr);
                        TiXmlElement *vNewElement = new TiXmlElement("B");
                        vNewElement->SetAttribute("ID", p);
                        vNewElement->LinkEndChild(vList);
                        vNewBndRegions->LinkEndChild(vNewElement);
                        vBndRegionIdList.insert(p);
                        vItem = vItem->NextSiblingElement();
                    }

                    // store original bnd region order
                    m_bndRegOrder[p] = vSeq;
                }
                vConditions->ReplaceChild(vBndRegions, *vNewBndRegions);
            }

            if (vBndConditions)
            {
                vItem = vBndConditions->FirstChildElement();
                while (vItem)
                {
                    std::set<int>::iterator x;
                    if ((x = vBndRegionIdList.find(atoi(vItem->Attribute(
                             "REF")))) != vBndRegionIdList.end())
                    {
                        vItem->SetAttribute("REF", *x);
                        vItem = vItem->NextSiblingElement();
                    }
                    else
                    {
                        TiXmlElement *tmp = vItem;
                        vItem             = vItem->NextSiblingElement();
                        vBndConditions->RemoveChild(tmp);
                    }
                }
            }
            root->LinkEndChild(vConditions);
        }

        // Distribute other sections of the XML to each process as is.
        TiXmlElement *vSrc =
            m_session->GetElement("Nektar")->FirstChildElement();
        while (vSrc)
        {
            std::string vName = boost::to_upper_copy(vSrc->ValueStr());
            if (vName != "GEOMETRY" && vName != "CONDITIONS")
            {
                root->LinkEndChild(new TiXmlElement(*vSrc));
            }
            vSrc = vSrc->NextSiblingElement();
        }

        // Save Mesh

        boost::format pad("P%1$07d.xml");
        pad % partitions[i];
        boost::filesystem::path pFilename(pad.str());

        boost::filesystem::path fullpath = pdirname / pFilename;
        doc.SaveFile(LibUtilities::PortablePath(fullpath));
    }
}

CompositeOrdering MeshGraphXml::CreateCompositeOrdering()
{
    CompositeOrdering ret;

    for (auto &c : m_meshComposites)
    {
        std::vector<unsigned int> ids;
        for (auto &elmt : c.second->m_geomVec)
        {
            ids.push_back(elmt->GetGlobalID());
        }
        ret[c.first] = ids;
    }

    return ret;
}

}
}
