////////////////////////////////////////////////////////////////////////////////
//
// File: CwipiExchange.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
//
// License for the specific language governing rights and limitations under
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

#include "CwipiExchange.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

namespace Nektar
{
namespace SolverUtils
{

using namespace std;

static void InterpCallback(
    const int entities_dim,
    const int n_local_vertex,
    const int n_local_element,
    const int n_local_polhyedra,
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
    SenderCouplings[sst.str()](interpField, distCoords);

    ASSERTL0(interpField.num_elements() == stride, "size mismatch");
    ASSERTL0(interpField[0].num_elements() == n_distant_point, "size mismatch");

    for (int i = 0; i < n_distant_point; i++)
    {
        for (int j = 0; j < stride; ++j)
        {
            ((double *)distant_field)[i * stride + j] = interpField[j][i];
        }
    }
}

CwipiCoupling::CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                             std::string name,
                             int outputFreq,
                             double geomTol)
    : m_couplingName(name), m_evalField(field), m_lastSend(-1E6),
      m_lastReceive(-1E6), m_points(NULL), m_coords(NULL), m_connecIdx(NULL),
      m_connec(NULL), m_rValsInterl(NULL), m_sValsInterl(NULL)
{
    ReadConfig(m_evalField->GetSession());

    cwipi_add_local_int_control_parameter("nSendVars", m_nSendVars);
    cwipi_add_local_int_control_parameter("nRecvVars", m_nRecvVars);
    cwipi_add_local_string_control_parameter(
        "recvFieldNames", m_config["RECEIVEVARIABLES"].c_str());
    cwipi_add_local_string_control_parameter("sendFieldNames",
                                             m_config["SENDVARIABLES"].c_str());

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cwipi_dump_application_properties();
    }

    m_spacedim = m_evalField->GetGraph()->GetSpaceDimension();

    if (m_filtWidth > 0)
    {
        m_filtWidth = 2 * M_PI / m_filtWidth;
        m_filtWidth = m_filtWidth * m_filtWidth;
    }

    //  Init Coupling
    cwipi_solver_type_t solver_type = CWIPI_SOLVER_CELL_VERTEX;
    cwipi_create_coupling(m_couplingName.c_str(),
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          m_config["REMOTENAME"].c_str(),
                          m_spacedim,
                          geomTol,
                          CWIPI_STATIC_MESH,
                          solver_type,
                          outputFreq,
                          "Ensight Gold",
                          "text");

    AnnounceMesh();

    if (m_nRecvVars > 0 and m_recvSteps > 0)
    {
        SetupReceive();
    }

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "locating..." << endl;
    }
    cwipi_locate(m_couplingName.c_str());

    if (m_nSendVars > 0 and m_sendSteps > 0)
    {
        SetupSend();
    }
}

CwipiCoupling::~CwipiCoupling()
{
    free(m_coords);
    free(m_points);
    free(m_connec);
    free(m_connecIdx);
    free(m_rValsInterl);
    free(m_sValsInterl);
}

void CwipiCoupling::ReadConfig(LibUtilities::SessionReaderSharedPtr session)
{
    // defaults
    m_config["REMOTENAME"]       = "precise";
    m_config["OVERSAMPLE"]       = "0";
    m_config["FILTERWIDTH"]      = "-1";
    m_config["RECEIVESTEPS"]     = "0";
    m_config["RECEIVEVARIABLES"] = "";
    m_config["SENDSTEPS"]        = "0";
    m_config["SENDVARIABLES"]    = "";
    m_config["DUMPRAW"]          = "0";
    m_config["SENDMETHOD"]       = "NEARESTNEIGHBOUR";

    ASSERTL0(session->DefinesElement("Nektar/Coupling"),
             "No Coupling config found");

    TiXmlElement *vCoupling = session->GetElement("Nektar/Coupling");
    ASSERTL0(vCoupling, "Invalid Coupling config");

    // TODO: set m_name here instead of in the constructor
    string nName;
    vCoupling->QueryStringAttribute("NAME", &nName);
    ASSERTL0(nName.size(), "No Coupling NAME attribute set");
    ASSERTL0(m_couplingName == nName, "Wrong Coupling name");

    TiXmlElement *element = vCoupling->FirstChildElement("I");
    while (element)
    {
        std::stringstream tagcontent;
        tagcontent << *element;
        // read the property name
        ASSERTL0(element->Attribute("PROPERTY"),
                 "Missing PROPERTY attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");
        std::string property = element->Attribute("PROPERTY");
        ASSERTL0(!property.empty(),
                 "PROPERTY attribute must be non-empty in XML "
                 "element: \n\t'" +
                     tagcontent.str() + "'");

        // make sure that solver property is capitalised
        std::string propertyUpper = boost::to_upper_copy(property);

        CouplingConfigMap::const_iterator x = m_config.find(propertyUpper);
        ASSERTL0(x != m_config.end(),
                 "Invalid PROPERTY attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");

        // read the value
        ASSERTL0(element->Attribute("VALUE"),
                 "Missing VALUE attribute in Coupling section "
                 "XML element: \n\t'" +
                     tagcontent.str() + "'");
        std::string value = element->Attribute("VALUE");
        ASSERTL0(!value.empty(),
                 "VALUE attribute must be non-empty in XML "
                 "element: \n\t'" +
                     tagcontent.str() + "'");

        // Set Variable
        m_config[propertyUpper] = value;

        element = element->NextSiblingElement("I");
    }

    // mangle config into variables. This is ugly
    ParseUtils::GenerateOrderedStringVector(
        m_config["RECEIVEVARIABLES"].c_str(), m_recvFieldNames);
    m_nRecvVars = m_recvFieldNames.size();

    ParseUtils::GenerateOrderedStringVector(m_config["SENDVARIABLES"].c_str(),
                                            m_sendFieldNames);
    m_nSendVars = m_sendFieldNames.size();

    m_recvSteps = boost::lexical_cast<int>(m_config["RECEIVESTEPS"]);
    m_sendSteps = boost::lexical_cast<int>(m_config["SENDSTEPS"]);

    m_filtWidth = boost::lexical_cast<NekDouble>(m_config["FILTERWIDTH"]);

    if (session->GetComm()->GetRank() == 0 &&
        session->DefinesCmdLineArgument("verbose") && m_config.size() > 0)
    {
        cout << "Coupling Config:" << endl;
        CouplingConfigMap::iterator x;
        for (x = m_config.begin(); x != m_config.end(); ++x)
        {
            cout << "\t" << x->first << " = " << x->second << endl;
        }

        cout << "\tRECEIVEVARIABLES (parsed) = [";
        std::vector<string>::const_iterator i;
        for (i = m_recvFieldNames.begin(); i != m_recvFieldNames.end(); ++i)
        {
            cout << *i << ", ";
        }
        cout << "]" << endl;

        cout << "\tSENDVARIABLES (parsed) = [";
        for (i = m_sendFieldNames.begin(); i != m_sendFieldNames.end(); ++i)
        {
            cout << *i << ", ";
        }
        cout << "]" << endl;
    }
}

void CwipiCoupling::SetupReceive()
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

void CwipiCoupling::SetupSend()
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
    SenderCouplings[sst.str()] =
        boost::bind(&CwipiCoupling::SendCallback, this, _1, _2);
    cwipi_set_interpolation_function(m_couplingName.c_str(), InterpCallback);
}

void CwipiCoupling::EvaluateFields(
    Array<OneD, Array<OneD, NekDouble> > interpField,
    Array<OneD, Array<OneD, NekDouble> > distCoords)
{
    int nOutPts = distCoords.num_elements();

    Array<OneD, NekDouble> Lcoords(m_spacedim, 0.0);
    for (int i = 0; i < nOutPts; ++i)
    {
        // Obtain Element and LocalCoordinate to interpolate
        int elmtid    = -1;
        NekDouble tol = NekConstants::kNekZeroTol;
        while (elmtid < 0 and tol <= 1E3 * NekConstants::kNekZeroTol)
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

        int offset =
            m_evalField->GetPhys_Offset(m_evalField->GetOffset_Elmt_Id(elmtid));

        for (int f = 0; f < m_nSendVars; ++f)
        {
            NekDouble value = m_evalField->GetExp(elmtid)->StdPhysEvaluate(
                Lcoords, m_sendField[f] + offset);

            ASSERTL0(!(boost::math::isnan)(value), "new value is not a number");
            interpField[f][i] = value;
        }
    }
}

void CwipiCoupling::SetupSendInterpolation()
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

    m_sendInterpolator =
        MemoryManager<SolverUtils::Interpolator>::AllocateSharedPtr(
            SolverUtils::eNearestNeighbour);
    m_sendInterpolator->CalcWeights(locatPts, distPts);
    m_sendInterpolator->PrintStatistics();
}

void CwipiCoupling::AnnounceMesh()
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
    int nElts = seggeom.size() + trigeom.size() + quadgeom.size() +
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

void CwipiCoupling::v_FinalizeCoupling(void)
{
    cwipi_delete_coupling(m_couplingName.c_str());
}

template <typename T>
void CwipiCoupling::AddElementsToMesh(T geom,
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
    typename T::iterator it;

    for (it = geom.begin(); it != geom.end(); it++)
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

void CwipiCoupling::PrintProgressbar(const int position, const int goal) const
{
    // print only every 2 percent
    if (int(100 * position / goal) % 2 == 0)
    {
        cout << "." << flush;
    }
}

void CwipiCoupling::SendCallback(
    Array<OneD, Array<OneD, NekDouble> > &interpField,
    Array<OneD, Array<OneD, NekDouble> > distCoords)
{
    ASSERTL0(interpField.num_elements() == m_nSendVars, "size mismatch");

    for (int i = 0; i < m_nSendVars; ++i)
    {
        interpField[i] = Array<OneD, NekDouble>(distCoords.num_elements());
    }

    if (m_config["SENDMETHOD"] == "NEARESTNEIGHBOUR")
    {
        if (not m_sendInterpolator)
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

void CwipiCoupling::SendFields(
    const int step,
    const NekDouble time,
    const Array<OneD, const Array<OneD, NekDouble> > &field)
{
    if (m_nSendVars < 1 or m_sendSteps < 1)
    {
        return;
    }

    ASSERTL1(m_nSendVars == field.num_elements(), "field size mismatch");

    if (step >= m_lastSend + m_sendSteps)
    {
        m_lastSend = step;

        m_sendField = field;

        if (m_evalField->GetComm()->GetRank() == 0 &&
            m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "sending fields at i = " << step << ", t = " << time
                 << endl;
        }

        Timer timer1;
        timer1.Start();

        int nNotLoc = 0;
        // workaround a bug in cwipi: receiving_field_name should be const char*
        // but
        // is char*
        char sendFN[10];
        strcpy(sendFN, "dummyName");

        cwipi_exchange(m_couplingName.c_str(),
                       "ex1",
                       m_nSendVars,
                       step,
                       time,
                       sendFN,
                       m_sValsInterl,
                       "",
                       NULL,
                       &nNotLoc);

        timer1.Stop();

        if (m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
        {
            cout << "Send total time: " << timer1.TimePerTest(1) << endl;
        }
    }
}

void CwipiCoupling::ReceiveFields(const int step,
                                  const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > &field)
{
    if (m_nRecvVars < 1 or m_recvSteps < 1)
    {
        return;
    }

    ASSERTL1(m_nRecvVars == field.num_elements(), "field size mismatch");

    int nq = m_evalField->GetTotPoints();

    // make sure we have sensible data in old/new field the first time this
    // method is called
    if (m_lastReceive < 0)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, field[i], 1, m_oldFields[i], 1);
            Vmath::Vcopy(nq, field[i], 1, m_newFields[i], 1);
        }
    }

    if (step >= m_lastReceive + m_recvSteps)
    {
        m_lastReceive = step;

        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, m_newFields[i], 1, m_oldFields[i], 1);
        }

        FetchFields(step, time, m_newFields);

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
                       field[i],
                       1);
    }
}

void CwipiCoupling::FetchFields(const int step,
                                const NekDouble time,
                                Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL1(m_nRecvVars == field.num_elements(), "field size mismatch");

    if (m_evalField->GetComm()->GetRank() == 0 &&
        m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "receiving fields at i = " << step << ", t = " << time << endl;
    }

    Timer timer1, timer2, timer3;
    timer1.Start();

    Array<OneD, NekDouble> tmpC(m_recvField->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble> > rVals(m_nRecvVars);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        rVals[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints());
    }

    int nNotLoc = 0;

    // workaround a bug in cwipi: receiving_field_name should be const char* but
    // is char*
    char recFN[10];
    strcpy(recFN, "dummyName");

    timer2.Start();
    cwipi_exchange(m_couplingName.c_str(),
                   "ex1",
                   m_nRecvVars,
                   step,
                   time,
                   "",
                   NULL,
                   recFN,
                   m_rValsInterl,
                   &nNotLoc);
    timer2.Stop();

    int tmp           = -1;
    const int *notLoc = &tmp;
    if (nNotLoc != 0)
    {
        cout << "WARNING: relocating " << nNotLoc << " of " << m_nPoints
             << " points" << endl;
        notLoc = cwipi_get_not_located_points(m_couplingName.c_str());

        // interpolate from m_evalField to m_recvField
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            m_evalField->FwdTrans(field[i], tmpC);
            m_recvField->BwdTrans(tmpC, rVals[i]);
        }
    }

    int locPos = 0;
    int intPos = 0;
    for (int j = 0; j < m_nRecvVars; ++j)
    {
        locPos = 0;
        intPos = 0;

        for (int i = 0; i < m_nPoints; ++i)
        {
            // cwipi indices start from 1
            if (notLoc[locPos] - 1 == i)
            {
                // keep the original value of field[j][i]
                locPos++;
            }
            else
            {
                rVals[j][i] = m_rValsInterl[intPos * m_nRecvVars + j];
                intPos++;
            }
        }
    }

    OverrrideFields(rVals);

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

            if (m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
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

    if (m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Receive total time: " << timer1.TimePerTest(1) << ", ";
        cout << "CWIPI time: " << timer2.TimePerTest(1) << endl;
    }
}

void CwipiCoupling::OverrrideFields(Array<OneD, Array<OneD, NekDouble> > &rVals)
{
    Timer timer1;
    timer1.Start();

    // HACK
    Array<OneD, Array<OneD, NekDouble> > x(3);
    x[0] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
    x[1] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
    x[2] = Array<OneD, NekDouble>(m_recvField->GetTotPoints(), 0.0);
    m_recvField->GetCoords(x[0], x[1], x[2]);
    NekDouble distsq = 1E23;
    NekDouble distsqNew = distsq;
    int refID = -1;
    for (int i = 0; i < m_recvField->GetTotPoints(); ++i)
    {
        distsqNew = pow(x[0][i] - 0.135, NekDouble(2)) + pow(x[1][i], NekDouble(2)) + pow(x[2][i], NekDouble(2));
        if (distsqNew < distsq)
        {
            distsq = distsqNew;
            refID = i;
        }
    }
    ASSERTL0(refID >= 0, "refID < 0");

    int ownerRank = -1;
    distsqNew = distsq;
    m_evalField->GetSession()->GetComm()->AllReduce(distsqNew, LibUtilities::ReduceMin);
    if(abs(distsqNew - distsq) <= NekConstants::kNekSqrtTol)
    {
        // this rank has the closest point
        ownerRank = m_evalField->GetSession()->GetComm()->GetRank();
    }
    m_evalField->GetSession()->GetComm()->AllReduce(ownerRank, LibUtilities::ReduceMax);
    ASSERTL0(ownerRank >= 0, "ownerRank < 0");


    Array<OneD, NekDouble> data(m_nRecvVars, 0.0);
    if (ownerRank == m_evalField->GetSession()->GetComm()->GetRank())
    {
        for (int j = 0; j < m_nRecvVars; ++j)
        {
            data[j] = rVals[j][refID];
        }
    }
    m_evalField->GetSession()->GetComm()->AllReduce(data, LibUtilities::ReduceSum);


    for (int i = 0; i < m_recvField->GetTotPoints(); ++i)
    {
        NekDouble tmp = pow(x[0][i] - 0.173, NekDouble(2)) + pow(x[1][i], NekDouble(2)) + pow(x[2][i], NekDouble(2));
        if (tmp <= pow(0.173-0.135, 2))
        {
            for (int j = 0; j < m_nRecvVars; ++j)
            {
                rVals[j][i] = data[j];
            }
        }
    }
}

void CwipiCoupling::DumpRawFields(const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > rVals)
{
    Timer timer1;
    timer1.Start();

#ifdef _WIN32
    // We need this to make sure boost::format has always
    // two digits in the exponents of Scientific notation.
    unsigned int old_exponent_format;
    old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
    filename = boost::str(boost::format(filename) % m_time);
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

    LibUtilities::PtsIO ptsIO(m_evalField->GetSession()->GetComm());
    LibUtilities::PtsFieldSharedPtr rvPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
            m_spacedim, m_recvFieldNames, tmp);
    ptsIO.Write(filename, rvPts);

    timer1.Stop();
    if (m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "DumpRawFields total time: " << timer1.TimePerTest(1) << endl;
    }
}
}
}
