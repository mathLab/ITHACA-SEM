////////////////////////////////////////////////////////////////////////////////
//
// File: CouplingCwipi.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: CWIPI Exchange class
//
////////////////////////////////////////////////////////////////////////////////

#include "CouplingCwipi.h"

#include <LibUtilities/BasicUtils/CsvIO.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>
#include <boost/functional/hash.hpp>

#define OUTPUT_FREQ 0

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

std::string CouplingCwipi::className =
    GetCouplingFactory().RegisterCreatorFunction(
        "Cwipi", CouplingCwipi::create, "Cwipi Coupling");

void CouplingCwipi::InterpCallback(
    const int entities_dim,
    const int n_local_vertex,
    const int n_local_element,
    const int n_local_polyhedra,
    const int n_distant_point,
    const double local_coordinates[],
    const int local_connectivity_index[],
    const int local_connectivity[],
    const int local_polyhedra_face_index[],
    const int local_polyhedra_cell_to_face_connectivity[],
    const int local_polyhedra_face_connectivity_index[],
    const int local_polyhedra_face_connectivity[],
    const double distant_points_coordinates[],
    const int distant_points_location[],
    const float distant_points_distance[],
    const int distant_points_barycentric_coordinates_index[],
    const double distant_points_barycentric_coordinates[],
    const int stride,
    const cwipi_solver_type_t solver_type,
    const void *local_field,
    void *distant_field)
{
    boost::ignore_unused(n_local_element, n_local_polyhedra, local_coordinates,
                         local_connectivity_index, local_connectivity,
                         local_polyhedra_face_index,
                         local_polyhedra_cell_to_face_connectivity,
                         local_polyhedra_face_connectivity_index,
                         local_polyhedra_face_connectivity,
                         distant_points_location,
                         distant_points_distance,
                         distant_points_barycentric_coordinates_index,
                         distant_points_barycentric_coordinates,
                         solver_type,
                         local_field);

    Array<OneD, Array<OneD, NekDouble> > interpField(stride);

    Array<OneD, Array<OneD, NekDouble> > distCoords(n_distant_point);
    for (int i = 0; i < n_distant_point; ++i)
    {
        distCoords[i] = Array<OneD, NekDouble>(3);
        for (int j = 0; j < 3; ++j)
        {
            distCoords[i][j] = distant_points_coordinates[3 * i + j];
        }
    }

    std::stringstream sst;
    sst << entities_dim << "," << n_local_vertex << "," << stride;
    SendCallbackMap[sst.str()](interpField, distCoords);

    ASSERTL0(interpField.size() == stride, "size mismatch");
    ASSERTL0(interpField[0].size() == n_distant_point, "size mismatch");

    for (int i = 0; i < n_distant_point; i++)
    {
        for (int j = 0; j < stride; ++j)
        {
            ((double *)distant_field)[i * stride + j] = interpField[j][i];
        }
    }
}

CouplingCwipi::CouplingCwipi(MultiRegions::ExpListSharedPtr field)
    : Coupling(field), m_sendHandle(-1), m_recvHandle(-1), m_lastSend(-1E6),
      m_lastReceive(-1E6), m_points(NULL), m_coords(NULL), m_connecIdx(NULL),
      m_connec(NULL), m_rValsInterl(NULL), m_sValsInterl(NULL)
{
    // defaults
    m_config["GEOMTOL"]      = "0.1";
    m_config["LOCALNAME"]    = "nektar";
    m_config["REMOTENAME"]   = "precise";
    m_config["OVERSAMPLE"]   = "0";
    m_config["FILTERWIDTH"]  = "-1";
    m_config["DUMPRAW"]      = "0";
    m_config["SENDMETHOD"]   = "EVALUATE";
    m_config["NOTLOCMETHOD"] = "KEEP";
}

void CouplingCwipi::v_Init()
{
    Coupling::v_Init();

    ReadConfig(m_evalField->GetSession());

    cwipi_add_local_int_control_parameter("nSendVars", m_nSendVars);
    cwipi_add_local_int_control_parameter("nRecvVars", m_nRecvVars);
    cwipi_add_local_string_control_parameter(
        "recvFieldNames", m_config["RECEIVEVARIABLES"].c_str());
    cwipi_add_local_string_control_parameter("sendFieldNames",
                                             m_config["SENDVARIABLES"].c_str());

    // MPI_TAG_UB is guaranteed to be at least 32767, so we make
    // sure m_recvTag < 32767. Only caveat: m_recvTag is not guaranteed to be
    // unique.
    m_recvTag = boost::hash<std::string>()(m_couplingName +
                    m_config["REMOTENAME"] + m_config["LOCALNAME"]) % 32767;

    cwipi_add_local_int_control_parameter("receiveTag", m_recvTag);

    m_spacedim = m_evalField->GetGraph()->GetSpaceDimension();

    if (m_filtWidth > 0)
    {
        m_filtWidth = 2 * M_PI / m_filtWidth;
        m_filtWidth = m_filtWidth * m_filtWidth;
    }

    //  Init Coupling
    cwipi_solver_type_t solver_type = CWIPI_SOLVER_CELL_VERTEX;
    NekDouble geom_tol = boost::lexical_cast<NekDouble>(m_config["GEOMTOL"]);
    cwipi_create_coupling(m_couplingName.c_str(),
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          m_config["REMOTENAME"].c_str(),
                          m_spacedim,
                          geom_tol,
                          CWIPI_STATIC_MESH,
                          solver_type,
                          OUTPUT_FREQ,
                          "Ensight Gold",
                          "text");
    cwipi_synchronize_control_parameter(m_config["REMOTENAME"].c_str());

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cwipi_dump_application_properties();
    }

    m_sendTag = cwipi_get_distant_int_control_parameter(
        m_config["REMOTENAME"].c_str(), "receiveTag");

    if (cwipi_has_int_parameter(m_config["REMOTENAME"].c_str(), "nRecvVars"))
    {
        int remoteNRecvVars = cwipi_get_distant_int_control_parameter(
            m_config["REMOTENAME"].c_str(), "nRecvVars");
        ASSERTL0(remoteNRecvVars == m_nSendVars,
                 "Number of local send vars different to remote received vars");
    }

    if (cwipi_has_int_parameter(m_config["REMOTENAME"].c_str(), "nSendVars"))
    {
        int remoteNSendVars = cwipi_get_distant_int_control_parameter(
            m_config["REMOTENAME"].c_str(), "nSendVars");
        ASSERTL0(remoteNSendVars == m_nRecvVars,
                 "Number of local receive vars different to remote sent vars");
    }

    AnnounceMesh();

    if (m_nRecvVars > 0 && m_recvSteps > 0)
    {
        SetupReceive();
    }

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "locating..." << endl;
    }
    cwipi_locate(m_couplingName.c_str());

    if (m_nSendVars > 0 && m_sendSteps > 0)
    {
        SetupSend();
    }

    if (m_nRecvVars > 0 && m_recvSteps > 0)
    {
        ReceiveStart();
    }
}

CouplingCwipi::~CouplingCwipi()
{
    free(m_coords);
    free(m_points);
    free(m_connec);
    free(m_connecIdx);
    free(m_rValsInterl);
    free(m_sValsInterl);
}

void CouplingCwipi::ReadConfig(LibUtilities::SessionReaderSharedPtr session)
{
    ASSERTL0(session->DefinesElement("Nektar/Coupling"),
             "No Coupling config found");

    m_config["LOCALNAME"] = session->GetCmdLineArgument<std::string>("cwipi");

    TiXmlElement *vCoupling = session->GetElement("Nektar/Coupling");
    ASSERTL0(vCoupling, "Invalid Coupling config");

    m_filtWidth = boost::lexical_cast<NekDouble>(m_config["FILTERWIDTH"]);
}

void CouplingCwipi::SetupReceive()
{
    int oversamp = boost::lexical_cast<int>(m_config["OVERSAMPLE"]);

    SpatialDomains::MeshGraphSharedPtr recvGraph =
        SpatialDomains::MeshGraph::Read(m_evalField->GetSession());
    recvGraph->SetExpansionsToPointOrder(
        oversamp + m_evalField->GetExp(0)->GetNumPoints(0));

    // TODO: DeclareCoeffPhysArrays
    switch (m_spacedim)
    {
        case 1:
        {
            m_recvField =
                MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph, "DefaultVar");
            break;
        }

        case 2:
        {

            m_recvField =
                MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph);
            break;
        }

        case 3:
        {

            m_recvField =
                MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(
                    m_evalField->GetSession(), recvGraph);
            break;
        }

        default:
        {
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
        }
    }

    m_oldFields = Array<OneD, Array<OneD, NekDouble> >(m_nRecvVars);
    m_newFields = Array<OneD, Array<OneD, NekDouble> >(m_nRecvVars);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        m_oldFields[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints());
        m_newFields[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints());
    }

    // define the quadrature points at which we want to receive data
    m_nPoints = m_recvField->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    coords[0] = Array<OneD, NekDouble>(m_nPoints);
    coords[1] = Array<OneD, NekDouble>(m_nPoints);
    coords[2] = Array<OneD, NekDouble>(m_nPoints);
    m_recvField->GetCoords(coords[0], coords[1], coords[2]);

    m_points = (double *)malloc(sizeof(double) * 3 * m_nPoints);
    ASSERTL1(m_points != NULL, "malloc failed for m_points");

    for (int i = 0; i < m_nPoints; ++i)
    {
        m_points[3 * i + 0] = double(coords[0][i]);

        if (m_spacedim > 1)
        {
            m_points[3 * i + 1] = double(coords[1][i]);
        }
        else
        {
            m_points[3 * i + 1] = 0.0;
        }

        if (m_spacedim > 2)
        {
            m_points[3 * i + 2] = double(coords[2][i]);
        }
        else
        {
            m_points[3 * i + 2] = 0.0;
        }
    }

    cwipi_set_points_to_locate(m_couplingName.c_str(), m_nPoints, m_points);

    m_rValsInterl = (double *)malloc(sizeof(double) * m_nPoints * m_nRecvVars);
    ASSERTL1(m_rValsInterl != NULL, "malloc failed for m_rValsInterl");
}

void CouplingCwipi::SetupSend()
{
    // this array is never used because of our send callback method
    m_sValsInterl = (double *)malloc(
        sizeof(double) * m_evalField->GetGraph()->GetNvertices() * m_nSendVars);
    ASSERTL1(m_sValsInterl != NULL, "malloc failed for m_sValsInterl");
    for (int i = 0; i < m_evalField->GetGraph()->GetNvertices() * m_nSendVars;
         ++i)
    {
        m_sValsInterl[i] = i;
    }

    // register this coupling as sender
    std::stringstream sst;
    sst << m_spacedim << "," << m_evalField->GetGraph()->GetNvertices() << ","
        << m_nSendVars;
    SendCallbackMap[sst.str()] = std::bind(&CouplingCwipi::SendCallback,
                                           this,
                                           std::placeholders::_1,
                                           std::placeholders::_2);
    cwipi_set_interpolation_function(m_couplingName.c_str(),
                                     CouplingCwipi::InterpCallback);
}

void CouplingCwipi::EvaluateFields(
    Array<OneD, Array<OneD, NekDouble> > interpField,
    Array<OneD, Array<OneD, NekDouble> > distCoords)
{
    int nOutPts = distCoords.size();

    Array<OneD, NekDouble> Lcoords(m_spacedim, 0.0);
    for (int i = 0; i < nOutPts; ++i)
    {
        // Obtain Element and LocalCoordinate to interpolate
        int elmtid    = -1;
        NekDouble tol = NekConstants::kNekZeroTol;
        while (elmtid < 0 && tol <= 1E3 * NekConstants::kNekZeroTol)
        {
            elmtid = m_evalField->GetExpIndex(distCoords[i], Lcoords, tol);
            tol *= 2;
        }
        if (tol > 2 * NekConstants::kNekZeroTol)
        {
            for (int j = 0; j < m_spacedim; ++j)
            {
                if (Lcoords[j] < -1 - 0.75 * NekConstants::kNekZeroTol)
                {
                    Lcoords[j] = -1;
                }
                if (Lcoords[j] > 1 + 0.75 * NekConstants::kNekZeroTol)
                {
                    Lcoords[j] = 1;
                }
            }
        }

        ASSERTL0(elmtid >= 0,
                 "no element found for (" +
                     boost::lexical_cast<string>(distCoords[i][0]) + ", " +
                     boost::lexical_cast<string>(distCoords[i][1]) + ", " +
                     boost::lexical_cast<string>(distCoords[i][2]) + ")");

        int offset = m_evalField->GetPhys_Offset(elmtid);

        for (int f = 0; f < m_nSendVars; ++f)
        {
            NekDouble value = m_evalField->GetExp(elmtid)->StdPhysEvaluate(
                Lcoords, m_sendField[f] + offset);

            ASSERTL0(!(boost::math::isnan)(value), "new value is not a number");
            interpField[f][i] = value;
        }
    }
}

void CouplingCwipi::SetupSendInterpolation()
{
    const double *distCoords =
        cwipi_get_distant_coordinates(m_couplingName.c_str());
    int nPts = cwipi_get_n_distant_points(m_couplingName.c_str());

    Array<OneD, Array<OneD, NekDouble> > local(3);
    for (int i = 0; i < 3; ++i)
    {
        local[i] = Array<OneD, NekDouble>(m_evalField->GetTotPoints(), 0.0);
    }
    m_evalField->GetCoords(local[0], local[1], local[2]);
    LibUtilities::PtsFieldSharedPtr locatPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, local);

    Array<OneD, Array<OneD, NekDouble> > dist(3);
    for (int i = 0; i < 3; ++i)
    {
        dist[i] = Array<OneD, NekDouble>(nPts);
        for (int j = 0; j < nPts; ++j)
        {
            dist[i][j] = distCoords[3 * j + i];
        }
    }
    LibUtilities::PtsFieldSharedPtr distPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, dist);

    LibUtilities::InterpMethod method = LibUtilities::eNearestNeighbour;
    if (boost::to_upper_copy(m_config["SENDMETHOD"]) == "SHEPARD")
    {
        method = LibUtilities::eShepard;
    }
    m_sendInterpolator =
        MemoryManager<FieldUtils::Interpolator>::AllocateSharedPtr(method);
    m_sendInterpolator->CalcWeights(locatPts, distPts);
    m_sendInterpolator->PrintStatistics();
}

void CouplingCwipi::AnnounceMesh()
{
    SpatialDomains::MeshGraphSharedPtr graph = m_evalField->GetGraph();

    // get Elements
    SpatialDomains::SegGeomMap seggeom;
    SpatialDomains::TriGeomMap trigeom;
    SpatialDomains::QuadGeomMap quadgeom;
    SpatialDomains::TetGeomMap tetgeom;
    SpatialDomains::PyrGeomMap pyrgeom;
    SpatialDomains::PrismGeomMap prismgeom;
    SpatialDomains::HexGeomMap hexgeom;
    if (m_spacedim == 1)
    {
        seggeom = graph->GetAllSegGeoms();
    }
    else if (m_spacedim == 2)
    {
        trigeom  = graph->GetAllTriGeoms();
        quadgeom = graph->GetAllQuadGeoms();
    }
    else if (m_spacedim == 3)
    {
        tetgeom   = graph->GetAllTetGeoms();
        pyrgeom   = graph->GetAllPyrGeoms();
        prismgeom = graph->GetAllPrismGeoms();
        hexgeom   = graph->GetAllHexGeoms();
    };

    int nVerts = graph->GetNvertices();
    int nElts  = seggeom.size() + trigeom.size() + quadgeom.size() +
                tetgeom.size() + pyrgeom.size() + prismgeom.size() +
                hexgeom.size();

    // allocate CWIPI arrays
    m_coords = (double *)malloc(sizeof(double) * 3 * nVerts);
    ASSERTL1(m_coords != NULL, "malloc failed for m_coords");
    int tmp = 2 * seggeom.size() + 3 * trigeom.size() + 4 * quadgeom.size() +
              4 * tetgeom.size() + 5 * pyrgeom.size() + 6 * prismgeom.size() +
              8 * hexgeom.size();
    m_connec = (int *)malloc(sizeof(int) * tmp);
    ASSERTL1(m_connec != NULL, "malloc failed for m_connec");
    m_connecIdx = (int *)malloc(sizeof(int) * (nElts + 1));
    ASSERTL1(m_connecIdx != NULL, "malloc failed for m_connecIdx");

    m_connecIdx[0] = 0;
    int coordsPos  = 0;
    int connecPos  = 0;
    int conidxPos  = 0;

    AddElementsToMesh(seggeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(trigeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(quadgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(tetgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(pyrgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(prismgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(hexgeom, coordsPos, connecPos, conidxPos);

    // output the mesh in tecplot format. If this works, CWIPI will be able
    // to process it, too
    /*
    cout << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\"" <<  endl;
    cout << "ZONE NODES=" <<  nVerts << ",ELEMENTS=" <<  nElts;
    cout << ",DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" <<  endl;
    for (int i = 0; i < nVerts; ++i)
    {
        cout << m_coords[3*i + 0] << " " << m_coords[3*i + 1] << " ";
        cout << m_coords[3*i + 2] << " " <<  1.0 << endl;
    }
    for (int i = 0; i < nElts; ++i)
    {
        cout << m_connec[i*4 + 0] << " " << m_connec[i*4 + 1] << " ";
        cout << m_connec[i*4 + 2] <<  " " <<  m_connec[i*4 + 3] <<  endl;
    }
    */

    cwipi_define_mesh(
        m_couplingName.c_str(), nVerts, nElts, m_coords, m_connecIdx, m_connec);
}

void CouplingCwipi::v_Finalize(void)
{
    cwipi_delete_coupling(m_couplingName.c_str());
}

template <typename T>
void CouplingCwipi::AddElementsToMesh(T geom,
                                      int &coordsPos,
                                      int &connecPos,
                                      int &conidxPos)
{
    // helper variables
    Array<OneD, NekDouble> x(3);
    SpatialDomains::PointGeomSharedPtr vert;
    int vertID;

    int kNverts = T::mapped_type::element_type::kNverts;

    // iterate over all elements
    for (auto it = geom.begin(); it != geom.end(); it++)
    {
        //  iterate over the elements vertices
        for (int j = 0; j < kNverts; ++j)
        {
            vert   = it->second->GetVertex(j);
            vertID = vert->GetVid();

            // check if we already stored the vertex
            if (m_vertMap.count(vertID) == 0)
            {
                //  store the vertex
                vert->GetCoords(x[0], x[1], x[2]);
                m_coords[3 * coordsPos + 0] = double(x[0]);
                m_coords[3 * coordsPos + 1] = double(x[1]);
                m_coords[3 * coordsPos + 2] = double(x[2]);

                // store the vertex position in the m_coords array
                m_vertMap[vertID] = coordsPos;
                coordsPos++;
            }

            m_connec[connecPos] = m_vertMap[vertID] + 1;
            connecPos++;
        }

        m_connecIdx[conidxPos + 1] = m_connecIdx[conidxPos] + kNverts;
        conidxPos++;
    }
}

void CouplingCwipi::SendCallback(
    Array<OneD, Array<OneD, NekDouble> > &interpField,
    Array<OneD, Array<OneD, NekDouble> > &distCoords)
{
    ASSERTL0(interpField.size() == m_nSendVars, "size mismatch");

    for (int i = 0; i < m_nSendVars; ++i)
    {
        interpField[i] = Array<OneD, NekDouble>(distCoords.size());
    }

    if (boost::to_upper_copy(m_config["SENDMETHOD"]) == "NEARESTNEIGHBOUR" ||
        boost::to_upper_copy(m_config["SENDMETHOD"]) == "SHEPARD")
    {
        if (!m_sendInterpolator)
        {
            SetupSendInterpolation();
        }

        LibUtilities::PtsFieldSharedPtr ptsIn =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                0, m_sendField);
        LibUtilities::PtsFieldSharedPtr ptsOut =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                0, interpField);
        m_sendInterpolator->Interpolate(ptsIn, ptsOut);
    }
    else
    {
        EvaluateFields(interpField, distCoords);
    }
}

void CouplingCwipi::v_Send(
    const int step,
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble> > &field,
    vector<string> &varNames)
{
    if (m_nSendVars < 1 || m_sendSteps < 1)
    {
        return;
    }

    if (step >= m_lastSend + m_sendSteps)
    {
        SendComplete();

        m_lastSend = step;

        vector<int> sendVarsToVars =
            GenerateVariableMapping(varNames, m_sendFieldNames);
        m_sendField = Array<OneD, Array<OneD, NekDouble> >(m_nSendVars);
        for (int i = 0; i < sendVarsToVars.size(); ++i)
        {
            m_sendField[i] = field[sendVarsToVars[i]];
        }

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "sending fields at i = " << step << ", t = " << time
                 << endl;
        }

        LibUtilities::Timer timer1;
        timer1.Start();

        char sendFN[10];
        strcpy(sendFN, "dummyName");

        cwipi_issend(m_couplingName.c_str(),
                     "ex1",
                     m_sendTag,
                     m_nSendVars,
                     step,
                     time,
                     sendFN,
                     m_sValsInterl,
                     &m_sendHandle);
        timer1.Stop();

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "Send total time: " << timer1.TimePerTest(1) << endl;
        }
    }
}

void CouplingCwipi::SendComplete()
{
    if (m_sendHandle < 0)
    {
        return;
    }

    LibUtilities::Timer timer1;
    timer1.Start();
    cwipi_wait_issend(m_couplingName.c_str(), m_sendHandle);
    timer1.Stop();

    // set to -1 so we dont try finishing a send before a new one was started
    m_sendHandle = -1;

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Send waiting time: " << timer1.TimePerTest(1) << endl;
    }
}

void CouplingCwipi::ReceiveStart()
{
    if (m_recvHandle >= 0)
    {
        return;
    }

    LibUtilities::Timer timer1;
    timer1.Start();
    // workaround a bug in cwipi: receiving_field_name should be const char* but
    // is char*
    char recFN[10];
    strcpy(recFN, "dummyName");

    cwipi_irecv(m_couplingName.c_str(),
                "ex1",
                m_recvTag,
                m_nRecvVars,
                0,
                0.0,
                recFN,
                m_rValsInterl,
                &m_recvHandle);
    timer1.Stop();

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Receive start time: " << timer1.TimePerTest(1) << endl;
    }
}

void CouplingCwipi::v_Receive(const int step,
                              const NekDouble time,
                              Array<OneD, Array<OneD, NekDouble> > &field,
                              vector<string> &varNames)
{
    if (m_nRecvVars < 1 || m_recvSteps < 1)
    {
        return;
    }

    Array<OneD, Array<OneD, NekDouble> > recvFields(m_nRecvVars);
    vector<int> recvVarsToVars =
        GenerateVariableMapping(varNames, m_recvFieldNames);
    ASSERTL1(m_nRecvVars == recvVarsToVars.size(), "field size mismatch");
    for (int i = 0; i < recvVarsToVars.size(); ++i)
    {
        recvFields[i] = field[recvVarsToVars[i]];
    }

    int nq = m_evalField->GetTotPoints();

    // make sure we have sensible data in old/new recvFields the first time this
    // method is called
    if (m_lastReceive < 0)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, recvFields[i], 1, m_oldFields[i], 1);
            Vmath::Vcopy(nq, recvFields[i], 1, m_newFields[i], 1);
        }
    }

    if (step >= m_lastReceive + m_recvSteps)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, m_newFields[i], 1, m_oldFields[i], 1);
        }

        ReceiveCwipi(step, time, m_newFields);
    }

    NekDouble fact =
        NekDouble(step - m_lastReceive + 1) / NekDouble(m_recvSteps);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        Vmath::Svtsvtp(nq,
                       fact,
                       m_newFields[i],
                       1,
                       (1 - fact),
                       m_oldFields[i],
                       1,
                       recvFields[i],
                       1);
    }
}

void CouplingCwipi::ReceiveCwipi(const int step,
                                 const NekDouble time,
                                 Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL1(m_nRecvVars == field.size(), "field size mismatch");

    if (m_nRecvVars < 1 || m_recvSteps < 1)
    {
        return;
    }

    if (step >= m_lastReceive + m_recvSteps)
    {
        m_lastReceive = step;

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "waiting for receive at i = " << step << ", t = " << time
                 << endl;
        }

        LibUtilities::Timer timer1, timer2, timer3;
        timer1.Start();

        Array<OneD, NekDouble> tmpC(m_recvField->GetNcoeffs());
        Array<OneD, Array<OneD, NekDouble> > rVals(m_nRecvVars);
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            rVals[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints());
        }

        timer2.Start();
        cwipi_wait_irecv(m_couplingName.c_str(), m_recvHandle);
        timer2.Stop();

        // set to -1 so we know we can start receiving again
        m_recvHandle = -1;

        int nNotLoc = cwipi_get_n_not_located_points(m_couplingName.c_str());
        Array<OneD, int> notLoc;
        if (nNotLoc != 0)
        {
            cout << "WARNING: relocating " << nNotLoc << " of " << m_nPoints
                 << " points" << endl;

            const int *tmp =
                cwipi_get_not_located_points(m_couplingName.c_str());
            notLoc = Array<OneD, int>(nNotLoc);
            for (int i = 0; i < nNotLoc; ++i)
            {
                notLoc[i] = tmp[i] - 1;
            }

            if (boost::to_upper_copy(m_config["NOTLOCMETHOD"]) == "KEEP")
            {
                // interpolate from m_evalField to m_recvField
                for (int i = 0; i < m_nRecvVars; ++i)
                {
                    m_evalField->FwdTrans(field[i], tmpC);
                    m_recvField->BwdTrans(tmpC, rVals[i]);
                }
            }
        }

        for (int i = 0, locPos = 0, intPos = 0; i < m_nPoints; ++i)
        {
            if (locPos < nNotLoc && notLoc[locPos] == i)
            {
                // keep the original value of field[j][i]
                locPos++;
            }
            else
            {
                for (int j = 0; j < m_nRecvVars; ++j)
                {
                    rVals[j][i] = m_rValsInterl[intPos * m_nRecvVars + j];
                }
                intPos++;
            }
        }

        if (boost::to_upper_copy(m_config["NOTLOCMETHOD"]) == "EXTRAPOLATE")
        {
            int doExtrapolate = 0;
            if (nNotLoc != 0)
            {
                doExtrapolate = 1;
            }
            m_evalField->GetSession()->GetComm()->AllReduce(
                doExtrapolate, LibUtilities::ReduceMax);
            if (doExtrapolate > 0)
            {
                ExtrapolateFields(rVals, notLoc);
            }
        }

        if (m_config["DUMPRAW"] != "0")
        {
            DumpRawFields(time, rVals);
        }

        if (m_filtWidth > 0)
        {
            for (int i = 0; i < m_nRecvVars; ++i)
            {
                timer3.Start();

                Array<OneD, NekDouble> forcing(m_nPoints);

                Array<OneD, Array<OneD, NekDouble> > Velocity(m_spacedim);
                for (int j = 0; j < m_spacedim; ++j)
                {
                    Velocity[j] = Array<OneD, NekDouble>(m_nPoints, 0.0);
                }

                Vmath::Smul(m_nPoints, -m_filtWidth, rVals[i], 1, forcing, 1);

                // Note we are using the
                // LinearAdvectionDiffusionReaction solver here
                // instead of HelmSolve since m_filtWidth is negative and
                // so matrices are not positive definite. Ideally
                // should allow for negative m_filtWidth coefficient in
                // HelmSolve
                m_recvField->LinearAdvectionDiffusionReactionSolve(
                    Velocity, forcing, tmpC, -m_filtWidth);

                m_evalField->BwdTrans(tmpC, field[i]);
                timer3.Stop();

                if (m_evalField->GetComm()->GetRank() == 0 &&
                    m_evalField->GetSession()->DefinesCmdLineArgument(
                        "verbose"))
                {
                    cout << "Smoother time (" << m_recvFieldNames[i]
                         << "): " << timer3.TimePerTest(1) << endl;
                }
            }
        }
        else
        {
            for (int i = 0; i < m_nRecvVars; ++i)
            {
                m_recvField->FwdTrans(rVals[i], tmpC);
                m_evalField->BwdTrans(tmpC, field[i]);
            }
        }
        timer1.Stop();

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "Receive total time: " << timer1.TimePerTest(1) << ", ";
            cout << "Receive waiting time: " << timer2.TimePerTest(1) << endl;
        }

        ReceiveStart();
    }
}

void CouplingCwipi::ExtrapolateFields(
    Array<OneD, Array<OneD, NekDouble> > &rVals, Array<OneD, int> &notLoc)
{
    LibUtilities::Timer timer1, timer2;
    timer1.Start();

    int totvars = 3 + m_nRecvVars;
    int nNotLoc = notLoc.size();
    int nranks  = m_evalField->GetSession()->GetComm()->GetSize();

    Array<OneD, int> thisNLoc(1, m_nPoints - nNotLoc);
    Array<OneD, int> sizeMap(nranks);
    m_evalField->GetSession()->GetComm()->AllGather(thisNLoc, sizeMap);

    Array<OneD, int> offsetMap(nranks);
    offsetMap[0] = 0;
    for (int i = 1; i < nranks; ++i)
    {
        offsetMap[i] = offsetMap[i - 1] + sizeMap[i - 1];
    }
    int totNLoc = offsetMap[nranks - 1] + sizeMap[nranks - 1];

    Array<OneD, Array<OneD, NekDouble> > allVals(totvars);
    for (int i = 0; i < 3; ++i)
    {
        allVals[i] = Array<OneD, NekDouble>(m_nPoints);
    }
    m_recvField->GetCoords(allVals[0], allVals[1], allVals[2]);

    if (m_spacedim < 3)
    {
        Vmath::Zero(m_nPoints, allVals[2], 1);
    }
    if (m_spacedim < 2)
    {
        Vmath::Zero(m_nPoints, allVals[1], 1);
    }

    for (int i = 0; i < m_nRecvVars; ++i)
    {
        allVals[3 + i] = rVals[i];
    }

    Array<OneD, Array<OneD, NekDouble> > locatedVals(totvars);
    for (int i = 0; i < totvars; ++i)
    {
        locatedVals[i] = Array<OneD, NekDouble>(totNLoc, -42.0);
    }

    // only copy points from allVals to locatedVals that were located
    int offset = offsetMap[m_evalField->GetSession()->GetComm()->GetRank()];
    for (int i = 0, intPos = 0, locPos = 0; i < m_nPoints; ++i)
    {
        if (locPos < nNotLoc && notLoc[locPos] == i)
        {
            // do nothing
            locPos++;
        }
        else
        {
            for (int k = 0; k < totvars; ++k)
            {
                locatedVals[k][offset + intPos] = allVals[k][i];
            }
            intPos++;
        }
    }

    // send all located points to all ranks. This is probably horribly expensive
    timer2.Start();
    for (int i = 0; i < totvars; ++i)
    {
        m_evalField->GetSession()->GetComm()->AllGatherv(
            locatedVals[i], sizeMap, offsetMap);
    }
    timer2.Stop();

    if (nNotLoc > 0)
    {
        LibUtilities::PtsFieldSharedPtr locatedPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                3, locatedVals);

        Array<OneD, Array<OneD, NekDouble> > notLocVals(totvars);
        for (int j = 0; j < totvars; ++j)
        {
            notLocVals[j] = Array<OneD, NekDouble>(nNotLoc);
            for (int i = 0; i < nNotLoc; ++i)
            {
                notLocVals[j][i] = allVals[j][notLoc[i]];
            }
        }
        LibUtilities::PtsFieldSharedPtr notlocPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                3, notLocVals);

        // perform a nearest neighbour interpolation from locatedVals to the not
        // located rVals
        if (!m_extrapInterpolator)
        {
            m_extrapInterpolator =
                MemoryManager<FieldUtils::Interpolator>::AllocateSharedPtr(
                    LibUtilities::eNearestNeighbour);
            m_extrapInterpolator->CalcWeights(locatedPts, notlocPts);
            m_extrapInterpolator->PrintStatistics();
        }
        m_extrapInterpolator->Interpolate(locatedPts, notlocPts);

        for (int j = 3; j < totvars; ++j)
        {
            for (int i = 0; i < nNotLoc; ++i)
            {
                allVals[j][notLoc[i]] = notLocVals[j][i];
            }
        }
    }

    timer1.Stop();
    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "ExtrapolateFields total time: " << timer1.TimePerTest(1);
        cout << " (AllGathervI: " << timer2.TimePerTest(1) << ")" << endl;
    }
}

void CouplingCwipi::DumpRawFields(const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > rVals)
{
    LibUtilities::Timer timer1;
    timer1.Start();

#if (defined _WIN32 && _MSC_VER < 1900)
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    filename            = boost::str(boost::format(filename) % m_time);
    _set_output_format(old_exponent_format);
#else
    std::string filename =
        boost::str(boost::format(m_config["DUMPRAW"]) % time);
#endif

    Array<OneD, Array<OneD, NekDouble> > tmp(m_nRecvVars + m_spacedim);
    for (int i = 0; i < 3; ++i)
    {
        tmp[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
    }
    m_recvField->GetCoords(tmp[0], tmp[1], tmp[2]);

    for (int i = 0; i < m_nRecvVars; ++i)
    {
        tmp[m_spacedim + i] = rVals[i];
    }

    LibUtilities::CsvIO csvIO(m_evalField->GetSession()->GetComm());
    LibUtilities::PtsFieldSharedPtr rvPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            m_spacedim, m_recvFieldNames, tmp);
    csvIO.Write(filename, rvPts);

    timer1.Stop();
    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "DumpRawFields total time: " << timer1.TimePerTest(1) << endl;
    }
}
}
}
