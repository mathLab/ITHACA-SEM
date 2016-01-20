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

#include <cwipi.h>


namespace Nektar
{

CwipiCoupling::CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                                   string name, int outputFreq, double geomTol) :
    Coupling(field, name),
    m_outputFormat("Ensight Gold"),
    m_outputFormatOption("text"),
    m_outputFreq(outputFreq),
    m_geomTol(geomTol)
{
    cwipi_dump_application_properties();

    //  Init Coupling
    cwipi_solver_type_t solver_type = CWIPI_SOLVER_CELL_VERTEX;

    cwipi_create_coupling(m_name.c_str(),
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          "precise",
                          3,
                          m_geomTol,
                          CWIPI_STATIC_MESH,
                          solver_type,
                          m_outputFreq,
                          m_outputFormat.c_str(),
                          m_outputFormatOption.c_str());

    SpatialDomains::MeshGraphSharedPtr graph = m_field->GetGraph();

    // get Elements
    SpatialDomains::TetGeomMap   tetgeom   = graph->GetAllTetGeoms();
    SpatialDomains::PyrGeomMap   pyrgeom   = graph->GetAllPyrGeoms();
    SpatialDomains::PrismGeomMap prismgeom = graph->GetAllPrismGeoms();
    SpatialDomains::HexGeomMap   hexgeom   = graph->GetAllHexGeoms();

    int nVerts = graph->GetNvertices();
    int nElts = tetgeom.size() + pyrgeom.size() + prismgeom.size() + hexgeom.size();
    m_nPoints = m_field->GetTotPoints();


    // allocate CWIPI arrays
    m_coords = (double *) malloc(sizeof(double) * 3 * nVerts);
    ASSERTL1(m_coords != NULL, "malloc failed for m_coords");
    m_points = (double *) malloc(sizeof(double) * 3 * m_nPoints);
    ASSERTL1(m_points != NULL, "malloc failed for m_points");
    int tmp = 4 * tetgeom.size() + 5 * pyrgeom.size() +
              6 * prismgeom.size() + 8 * hexgeom.size();
    m_connec = (int *) malloc(sizeof(int) * tmp);
    ASSERTL1(m_connec != NULL, "malloc failed for m_connec");
    m_connecIdx = (int *) malloc(sizeof(int) * (nElts + 1));
    ASSERTL1(m_connecIdx != NULL, "malloc failed for m_connecIdx");

    m_connecIdx[0] = 0;
    int coordsPos = 0;
    int connecPos = 0;
    int conidxPos = 0;

    AddElementsToMesh(tetgeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(pyrgeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(prismgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(hexgeom,   coordsPos, connecPos, conidxPos);


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

    cwipi_define_mesh(m_name.c_str(),
                      nVerts,
                      nElts,
                      m_coords,
                      m_connecIdx,
                      m_connec);

    Array<OneD, NekDouble> x0(m_nPoints);
    Array<OneD, NekDouble> x1(m_nPoints);
    Array<OneD, NekDouble> x2(m_nPoints);
    m_field->GetCoords(x0, x1, x2);

    for (int i = 0; i < m_nPoints; ++i)
    {
        m_points[3 * i + 0] = double(x0[i]);
        m_points[3 * i + 1] = double(x1[i]);
        m_points[3 * i + 2] = double(x2[i]);
    }

    cwipi_set_points_to_locate(m_name.c_str(), m_nPoints, m_points);
}


CwipiCoupling::~CwipiCoupling()
{
    free(m_coords);
    free(m_points);
    free(m_connec);
    free(m_connecIdx);

    cwipi_delete_coupling(m_name.c_str());
}


template <typename T>
void CwipiCoupling::AddElementsToMesh(T geom, int &coordsPos, int &connecPos,
        int &conidxPos)
{
    // helper variables
    Array<OneD, NekDouble>  x(3);
    SpatialDomains::PointGeomSharedPtr vert;
    int vertID;

    // iterate over all elements
    typename T::iterator it;

    for (it = geom.begin(); it != geom.end(); it++)
    {
        //  iterate over the elements vertices
        for (int j = 0; j < it->second->kNverts; ++j)
        {
            vert = it->second->GetVertex(j);
            vertID = vert->GetVid();

            // check if we already stored the vertex
            if (m_vertMap.count(vertID) ==  0)
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

        m_connecIdx[conidxPos + 1] = m_connecIdx[conidxPos] + it->second->kNverts;
        conidxPos++;
    }
}








CwipiExchange::CwipiExchange(SolverUtils::CouplingSharedPointer coupling,
                                   string name, int nEVars) :
    Exchange(coupling,  name),
    m_nEVars(nEVars)
{
    int nPoints = m_coupling->GetNPoints();

    m_rValsInterl = (double *) malloc(sizeof(double) * nPoints * m_nEVars);
    ASSERTL1(m_rValsInterl != NULL, "malloc failed for m_rValsInterl");
}


CwipiExchange::~CwipiExchange()
{
    free(m_rValsInterl);
}


void CwipiExchange::v_SendFields(const int step, const NekDouble time,
                                    Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL0(false, "not implemented yet")
}


void CwipiExchange::v_ReceiveFields(const int step, const NekDouble time,
                                       Array<OneD, Array<OneD, NekDouble> > &field)
{
    static int lastUdate = -1;

    if (step <= lastUdate)
    {
        cout <<  "returning" << endl;
        return;
    }

    lastUdate = step;

    int nPoints = m_coupling->GetNPoints();
    ASSERTL1(nPoints ==  field[0].num_elements(), "field size mismatch");
    ASSERTL1(m_nEVars ==  field.num_elements(), "field size mismatch");

    cout << "receiving fields at i = " << step << ", t = " << time << endl;

    int nNotLoc;
    cwipi_exchange(m_coupling->GetName().c_str(),
                   m_name.c_str(),
                   m_nEVars,
                   step,
                   time,
                   "",
                   NULL,
                   m_recvFieldName.c_str(),
                   m_rValsInterl,
                   &nNotLoc);

    if (nNotLoc !=  0)
    {
        cout << "WARNING: " << nNotLoc << " not located points found" <<  endl;

        int *notLoc;
        notLoc = (int *) malloc(sizeof(int) * nNotLoc);
        ASSERTL1(notLoc != NULL, "malloc failed for notLoc");

//         retVal = PCW_Get_not_located_points(cp_name, nNotLoc, notLoc);
//         ASSERTL0(retVal == 0, "ERROR: PCW_Get_not_located_points failed with code" + retVal);
//         retVal = PCW_Reorder(m_rValsInterl, nPoints, m_nEVars, 0.0, notLoc, nNotLoc);
//         ASSERTL0(retVal == 0, "ERROR: PCW_Reorder failed with code" + retVal);

        free(notLoc);
    }

    // retain proper fields from the interlaced array
    for (int j = 0; j < m_nEVars; ++j)
    {
        for (int i = 0; i < nPoints; ++i)
        {
            field[j][i] = m_rValsInterl[i * m_nEVars + j];
        }
    }
}



}
