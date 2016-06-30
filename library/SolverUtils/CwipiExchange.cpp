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

#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/Timer.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

#include <cwipi.h>

namespace Nektar
{
namespace SolverUtils
{

CwipiCoupling::CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                             string name,
                             int outputFreq,
                             double geomTol)
    : m_couplingName(name), m_evalField(field), m_nSendVars(0), m_nRecvVars(6),
      m_lastUpdate(-1E23), m_points(NULL), m_coords(NULL), m_connecIdx(NULL),
      m_connec(NULL), m_rValsInterl(NULL), m_sValsInterl(NULL)
{
    // HACK: m_nSendVars(0), m_nRecvVars(6)
    ReadConfig(m_evalField->GetSession());

    cwipi_dump_application_properties();

    m_spacedim = m_evalField->GetGraph()->GetSpaceDimension();

    m_filtWidth = boost::lexical_cast<NekDouble>(m_config["FILTERWIDTH"]);
    if (m_filtWidth > 0)
    {
        m_filtWidth = 2 * M_PI / m_filtWidth;
        m_filtWidth = m_filtWidth * m_filtWidth;
    }

    m_recvSteps = boost::lexical_cast<int>(m_config["RECEIVESTEPS"]);

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

    SetupRecvFields();

    AnnounceMesh();

    AnnounceRecvPoints();
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
    m_config["REMOTENAME"]   = "precise";
    m_config["OVERSAMPLE"]   = "0";
    m_config["FILTERWIDTH"]  = "-1";
    m_config["RECEIVESTEPS"] = "1";
    m_config["DUMPRAW"]      = "0";

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

    if (session->GetComm()->GetRank() == 0 && m_config.size() > 0)
    {
        cout << "Coupling Config:" << endl;
        CouplingConfigMap::iterator x;
        for (x = m_config.begin(); x != m_config.end(); ++x)
        {
            cout << "\t" << x->first << " = " << x->second << endl;
        }
        cout << endl;
    }
}

void CwipiCoupling::SetupRecvFields()
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

void CwipiCoupling::AnnounceRecvPoints()
{
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

    if (m_nRecvVars > 0)
    {
        m_rValsInterl =
            (double *)malloc(sizeof(double) * m_nPoints * m_nRecvVars);
        ASSERTL1(m_rValsInterl != NULL, "malloc failed for m_rValsInterl");
    }
    if (m_nSendVars > 0)
    {
        m_sValsInterl =
            (double *)malloc(sizeof(double) * m_nPoints * m_nSendVars);
        ASSERTL1(m_sValsInterl != NULL, "malloc failed for m_sValsInterl");
    }
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

void CwipiCoupling::SendFields(const int step,
                               const NekDouble time,
                               Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL0(false, "not implemented yet")
}

void CwipiCoupling::ReceiveFields(const int step,
                                  const NekDouble time,
                                  const NekDouble timestep,
                                  Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL1(m_nRecvVars == field.num_elements(), "field size mismatch");

    int nq = m_evalField->GetTotPoints();

    // make sure we have sensible data in old/new field the first time this
    // method is called
    if (m_lastUpdate < 0.0)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, field[i], 1, m_oldFields[i], 1);
            Vmath::Vcopy(nq, field[i], 1, m_newFields[i], 1);
        }
    }

    if (time >= m_lastUpdate + m_recvSteps * timestep)
    {
        m_lastUpdate = time;

        for (int i = 0; i < m_nRecvVars; ++i)
        {
            Vmath::Vcopy(nq, m_newFields[i], 1, m_oldFields[i], 1);
        }

        FetchFields(step, time, m_newFields);
    }

    NekDouble fact =
        (time - m_lastUpdate + timestep) / (m_recvSteps * timestep);
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

    cout << "receiving fields at i = " << step << ", t = " << time << endl;

    Timer timer1, timer2, timer3;
    timer1.Start();

    Array<OneD, NekDouble> tmpC(m_recvField->GetNcoeffs());
    Array<OneD, Array<OneD, NekDouble> > rVals(m_nRecvVars);
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        rVals[i] = Array<OneD, NekDouble>(m_recvField->GetTotPoints());
    }

    // interpolate from m_evalField to m_recvField
    for (int i = 0; i < m_nRecvVars; ++i)
    {
        m_evalField->FwdTrans(field[i], tmpC);
        m_recvField->BwdTrans(tmpC, rVals[i]);
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

    timer3.Start();
    if (boost::lexical_cast<bool>(m_config["DUMPRAW"]))
    {
        DumpRawFields(time, rVals);
    }

    if (m_filtWidth > 0)
    {
        for (int i = 0; i < m_nRecvVars; ++i)
        {
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
    timer3.Stop();
    timer1.Stop();

    if (m_evalField->GetSession()->DefinesCmdLineArgument("verbose"))
    {
        cout << "Receive total time: " << timer1.TimePerTest(1) << ", ";
        cout << "CWIPI time: " << timer2.TimePerTest(1) << ", ";
        cout << "Smoother time: " << timer3.TimePerTest(1) << endl;
    }
}

void CwipiCoupling::DumpRawFields(const NekDouble time,
                                  Array<OneD, Array<OneD, NekDouble> > rVals)
{
    std::string filename =
        "rawFields_" + boost::lexical_cast<std::string>(time) + ".pts";

    // HACK
    std::vector<std::string> recvFieldNames;
    recvFieldNames.push_back("u0");
    recvFieldNames.push_back("v0");
    recvFieldNames.push_back("w0");
    recvFieldNames.push_back("p0");
    recvFieldNames.push_back("rho0");
    recvFieldNames.push_back("S");

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
            m_spacedim, recvFieldNames, tmp);
    ptsIO.Write(filename, rvPts);
}
}
}
