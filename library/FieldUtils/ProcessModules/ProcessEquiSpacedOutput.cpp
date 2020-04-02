////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessEquiSpacedOutput.cpp
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
//  Description: Set up fields as interpolation to equispaced output
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/Interp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>

#include "ProcessEquiSpacedOutput.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessEquiSpacedOutput::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "equispacedoutput"),
        ProcessEquiSpacedOutput::create,
        "Write data as equi-spaced output using simplices to represent the "
        "data for connecting points");

ProcessEquiSpacedOutput::ProcessEquiSpacedOutput(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["tetonly"] =
        ConfigOption(true, "0", "Only process tetrahedral elements");

    m_config["modalenergy"] =
        ConfigOption(true, "0", "Write output as modal energy");
}

ProcessEquiSpacedOutput::~ProcessEquiSpacedOutput()
{
}

void ProcessEquiSpacedOutput::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nel = m_f->m_exp[0]->GetExpSize();
    if (!nel)
    {
        // Create empty PtsField
        int nfields = m_f->m_variables.size();
        int coordim = 3;

        Array<OneD, Array<OneD, NekDouble> > pts(nfields + coordim);
        for (int i = 0; i < nfields + coordim; ++i)
        {
            pts[i] = Array<OneD, NekDouble>(0);
        }
        vector<Array<OneD, int> > ptsConn;

        m_f->m_fieldPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
                coordim, m_f->m_variables, pts);
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTetBlock);
        m_f->m_fieldPts->SetConnectivity(ptsConn);
        return;
    }

    int coordim  = m_f->m_exp[0]->GetCoordim(0);
    int shapedim = m_f->m_exp[0]->GetExp(0)->GetShapeDimension();
    int npts     = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > coords(3);

    // Check if we have a homogeneous expansion
    bool homogeneous1D = false;
    if (m_f->m_numHomogeneousDir == 1)
    {
        coordim++;
        shapedim++;
        homogeneous1D = true;
    }
    else if (m_f->m_numHomogeneousDir == 2)
    {
        ASSERTL0(false, "Homegeneous2D case not supported");
    }

    // set up the number of points in each element
    int newpoints    = 0;
    int newtotpoints = 0;

    Array<OneD, int> conn;
    int prevNcoeffs = 0;
    int prevNpoints = 0;
    int cnt         = 0;

    // identify face 1 connectivity for prisms
    map<int, StdRegions::Orientation> face0orient;
    set<int> prismorient;
    LocalRegions::ExpansionSharedPtr e;

    // prepare PtsField
    vector<int> ppe;
    vector<Array<OneD, int> > ptsConn;
    int nfields;

    for (int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if (e->DetShapeType() == LibUtilities::ePrism)
        {
            StdRegions::Orientation forient = e->GetForient(0);
            int fid                         = e->GetGeom()->GetFid(0);
            if (face0orient.count(fid))
            { // face 1 meeting face 1 so reverse this id
                prismorient.insert(i);
            }
            else
            {
                // just store if Dir 1 is fwd or bwd
                if ((forient == StdRegions::eDir1BwdDir1_Dir2FwdDir2) ||
                    (forient == StdRegions::eDir1BwdDir1_Dir2BwdDir2) ||
                    (forient == StdRegions::eDir1BwdDir2_Dir2FwdDir1) ||
                    (forient == StdRegions::eDir1BwdDir2_Dir2BwdDir1))
                {
                    face0orient[fid] = StdRegions::eBwd;
                }
                else
                {
                    face0orient[fid] = StdRegions::eFwd;
                }
            }
        }
    }

    for (int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if (e->DetShapeType() == LibUtilities::ePrism)
        {
            int fid = e->GetGeom()->GetFid(2);
            // check to see if face 2 meets face 1
            if (face0orient.count(fid))
            {
                // check to see how face 2 is orientated
                StdRegions::Orientation forient2 = e->GetForient(2);
                StdRegions::Orientation forient0 = face0orient[fid];

                // If dir 1 or forient2 is bwd then check agains
                // face 1 value
                if ((forient2 == StdRegions::eDir1BwdDir1_Dir2FwdDir2) ||
                    (forient2 == StdRegions::eDir1BwdDir1_Dir2BwdDir2) ||
                    (forient2 == StdRegions::eDir1BwdDir2_Dir2FwdDir1) ||
                    (forient2 == StdRegions::eDir1BwdDir2_Dir2BwdDir1))
                {
                    if (forient0 == StdRegions::eFwd)
                    {
                        prismorient.insert(i);
                    }
                }
                else
                {
                    if (forient0 == StdRegions::eBwd)
                    {
                        prismorient.insert(i);
                    }
                }
            }
        }
    }

    for (int i = 0; i < nel; ++i)
    {
        e = m_f->m_exp[0]->GetExp(i);
        if (m_config["tetonly"].as<bool>())
        {
            if (m_f->m_exp[0]->GetExp(i)->DetShapeType() !=
                LibUtilities::eTetrahedron)
            {
                continue;
            }
        }

        switch (e->DetShapeType())
        {
            case LibUtilities::eSegment:
            {
                int npoints0 = e->GetBasis(0)->GetNumPoints();

                newpoints =
                    LibUtilities::StdSegData::getNumberOfCoefficients(npoints0);
            }
            break;
            case LibUtilities::eTriangle:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np  = max(np0, np1);
                newpoints =
                    LibUtilities::StdTriData::getNumberOfCoefficients(np, np);
            }
            break;
            case LibUtilities::eQuadrilateral:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np  = max(np0, np1);

                newpoints =
                    LibUtilities::StdQuadData::getNumberOfCoefficients(np, np);
            }
            break;
            case LibUtilities::eTetrahedron:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np  = max(np0, max(np1, np2));

                newpoints = LibUtilities::StdTetData::getNumberOfCoefficients(
                    np, np, np);
            }
            break;
            case LibUtilities::ePrism:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np  = max(np0, max(np1, np2));

                newpoints = LibUtilities::StdPrismData::getNumberOfCoefficients(
                    np, np, np);
            }
            break;
            case LibUtilities::ePyramid:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np  = max(np0, max(np1, np2));

                newpoints = LibUtilities::StdPyrData::getNumberOfCoefficients(
                    np, np, np);
            }
            break;
            case LibUtilities::eHexahedron:
            {
                int np0 = e->GetBasis(0)->GetNumPoints();
                int np1 = e->GetBasis(1)->GetNumPoints();
                int np2 = e->GetBasis(2)->GetNumPoints();
                int np  = max(np0, max(np1, np2));

                newpoints = LibUtilities::StdHexData::getNumberOfCoefficients(
                    np, np, np);
            }
            break;
            default:
            {
                ASSERTL0(false, "Points not known");
            }
        }

        ppe.push_back(newpoints);
        newtotpoints += newpoints;

        if (e->DetShapeType() == LibUtilities::ePrism)
        {
            bool standard = true;

            if (prismorient.count(i))
            {
                standard = false; // reverse direction
            }

            e->GetSimplexEquiSpacedConnectivity(conn, standard);
        }
        else
        {

            if ((prevNcoeffs != e->GetNcoeffs()) ||
                (prevNpoints != e->GetTotPoints()))
            {
                prevNcoeffs = e->GetNcoeffs();
                prevNpoints = e->GetTotPoints();

                e->GetSimplexEquiSpacedConnectivity(conn);
            }
        }
        Array<OneD, int> newconn(conn.size());
        for (int j = 0; j < conn.size(); ++j)
        {
            newconn[j] = conn[j] + cnt;
        }

        ptsConn.push_back(newconn);
        cnt += newpoints;
    }

    nfields = m_f->m_variables.size();

    Array<OneD, Array<OneD, NekDouble> > pts(nfields + coordim);

    for (int i = 0; i < nfields + coordim; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(newtotpoints);
    }

    // Interpolate coordinates
    for (int i = 0; i < coordim; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(npts);
    }

    for (int i = coordim; i < 3; ++i)
    {
        coords[i] = NullNekDouble1DArray;
    }

    m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

    Array<OneD, NekDouble> tmp;

    for (int n = 0; n < coordim; ++n)
    {
        cnt      = 0;
        int cnt1 = 0;
        for (int i = 0; i < nel; ++i)
        {
            m_f->m_exp[0]->GetExp(i)->PhysInterpToSimplexEquiSpaced(
                coords[n] + cnt, tmp = pts[n] + cnt1);
            cnt1 += ppe[i];
            cnt += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
        }
    }

    for (int n = 0; n < m_f->m_variables.size(); ++n)
    {
        cnt      = 0;
        int cnt1 = 0;

        if (m_config["modalenergy"].as<bool>())
        {
            Array<OneD, const NekDouble> phys = m_f->m_exp[n]->GetPhys();
            for (int i = 0; i < nel; ++i)
            {
                GenOrthoModes(i, phys + cnt, tmp = pts[coordim + n] + cnt1);
                cnt1 += ppe[i];
                cnt += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
            }
        }
        else
        {
            Array<OneD, const NekDouble> phys = m_f->m_exp[n]->GetPhys();
            for (int i = 0; i < nel; ++i)
            {
                m_f->m_exp[0]->GetExp(i)->PhysInterpToSimplexEquiSpaced(
                    phys + cnt, tmp = pts[coordim + n] + cnt1);
                cnt1 += ppe[i];
                cnt += m_f->m_exp[0]->GetExp(i)->GetTotPoints();
            }
        }
    }

    m_f->m_fieldPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        coordim, m_f->m_variables, pts);
    if (shapedim == 1)
    {
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsSegBlock);
    }

    if (shapedim == 2)
    {
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTriBlock);
    }
    else if (shapedim == 3)
    {
        m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTetBlock);
    }
    m_f->m_fieldPts->SetConnectivity(ptsConn);
    if (homogeneous1D)
    {
        SetHomogeneousConnectivity();
    }

    // Clear m_exp
    m_f->m_exp = vector<MultiRegions::ExpListSharedPtr>();
}

void ProcessEquiSpacedOutput::SetHomogeneousConnectivity(void)
{
    LocalRegions::ExpansionSharedPtr e;
    int nel          = m_f->m_exp[0]->GetPlane(0)->GetExpSize();
    int nPlanes      = m_f->m_exp[0]->GetZIDs().size();
    int npts         = m_f->m_fieldPts->GetNpoints();
    int nptsPerPlane = npts / nPlanes;
    int coordim      = 3;

    if (m_f->m_exp[0]->GetHomogeneousBasis()->GetBasisType() ==
            LibUtilities::eFourier &&
        m_f->m_exp[0]->GetHomogeneousBasis()->GetPointsType() ==
            LibUtilities::eFourierEvenlySpaced)
    {
        // Write points with extra plane
        Array<OneD, Array<OneD, NekDouble> > pts;
        m_f->m_fieldPts->GetPts(pts);
        Array<OneD, Array<OneD, NekDouble> > newPts(pts.size());
        for (int i = 0; i < pts.size(); i++)
        {
            newPts[i] = Array<OneD, NekDouble>(npts + nptsPerPlane);
            // Copy old points
            Vmath::Vcopy(npts, pts[i], 1, newPts[i], 1);
            // Get new plane
            Array<OneD, NekDouble> extraPlane(nptsPerPlane, newPts[i] + npts);
            if (m_f->m_comm->GetColumnComm()->GetSize() == 1)
            {
                if ( i == (coordim-1))
                {
                    // z-coordinate
                    Vmath::Sadd(nptsPerPlane, m_f->m_exp[0]->GetHomoLen(),
                                pts[i], 1, extraPlane, 1);
                }
                else
                {
                    Vmath::Vcopy(nptsPerPlane, pts[i], 1, extraPlane, 1);
                }
            }
            else
            {
                // Determine to and from rank for communication
                int size     = m_f->m_comm->GetColumnComm()->GetSize();
                int rank     = m_f->m_comm->GetColumnComm()->GetRank();
                int fromRank = (rank + 1) % size;
                int toRank   = (rank == 0) ? size - 1 : rank - 1;
                // Communicate using SendRecv
                Array<OneD, NekDouble> send(nptsPerPlane, pts[i]);
                m_f->m_comm->GetColumnComm()->SendRecv(toRank, send, fromRank,
                                                       extraPlane);
                // Adjust z-coordinate of last process
                if (i == (coordim-1) && (rank == size - 1) )
                {
                    Vmath::Sadd(nptsPerPlane, m_f->m_exp[0]->GetHomoLen(),
                                extraPlane, 1, extraPlane, 1);
                }
            }
        }
        m_f->m_fieldPts->SetPts(newPts);
        // Now extend plane connectivity
        vector<Array<OneD, int> > oldConn;
        Array<OneD, int> conn;
        m_f->m_fieldPts->GetConnectivity(oldConn);
        vector<Array<OneD, int> > newConn = oldConn;
        int connPerPlane = oldConn.size() / nPlanes;
        for (int i = 0; i < connPerPlane; i++)
        {
            conn = Array<OneD, int>(oldConn[i].size());
            for (int j = 0; j < conn.size(); j++)
            {
                conn[j] = oldConn[i][j] + npts;
            }
            newConn.push_back(conn);
        }
        m_f->m_fieldPts->SetConnectivity(newConn);

        nPlanes++;
        npts += nptsPerPlane;
    }

    vector<Array<OneD, int> > oldConn;
    vector<Array<OneD, int> > newConn;
    Array<OneD, int>          conn;
    m_f->m_fieldPts->GetConnectivity(oldConn);

    // 2DH1D case
    if (m_f->m_fieldPts->GetPtsType() == LibUtilities::ePtsTriBlock)
    {
        for(int i = 0; i < nel; ++i)
        {
            int nLines = oldConn[i].size()/2;
            // Create array for new connectivity
            // (2 triangles between each plane for each line)
            conn = Array<OneD, int> (2*3*nLines*(nPlanes-1));
            int cnt = 0;
            for(int n = 0; n<nLines; n++)
            {
                // Define new connectivity with triangles
                for ( int p = 0; p<nPlanes-1; p++)
                {
                    conn[cnt++] = oldConn[i+  p  *nel][2*n+0];
                    conn[cnt++] = oldConn[i+  p  *nel][2*n+1];
                    conn[cnt++] = oldConn[i+(p+1)*nel][2*n+0];

                    conn[cnt++] = oldConn[i+  p  *nel][2*n+1];
                    conn[cnt++] = oldConn[i+(p+1)*nel][2*n+0];
                    conn[cnt++] = oldConn[i+(p+1)*nel][2*n+1];
                }
            }
            newConn.push_back(conn);
        }
    }
    else if(m_f->m_fieldPts->GetPtsType() == LibUtilities::ePtsTetBlock)
    {
        // Get maximum number of points per direction in the mesh
        int maxN = 0;
        for(int i = 0; i < nel; ++i)
        {
            e = m_f->m_exp[0]->GetPlane(0)->GetExp(i);

            int np0 = e->GetBasis(0)->GetNumPoints();
            int np1 = e->GetBasis(1)->GetNumPoints();
            int np = max(np0,np1);
            maxN = max(np, maxN);
        }

        // Create a global numbering for points in plane 0
        Array<OneD, int> vId(nptsPerPlane);
        int cnt1=0;
        int cnt2=0;
        for(int i = 0; i < nel; ++i)
        {
            e = m_f->m_exp[0]->GetPlane(0)->GetExp(i);
            int np0 = e->GetBasis(0)->GetNumPoints();
            int np1 = e->GetBasis(1)->GetNumPoints();
            int np = max(np0,np1);
            switch(e->DetShapeType())
            {
            case LibUtilities::eTriangle:
                {
                    int newpoints = LibUtilities::StdTriData::
                                        getNumberOfCoefficients(np,np);

                    // Vertex numbering
                    vId[cnt1]                 = e->GetGeom()->GetVid(0);
                    vId[cnt1+np-1]            = e->GetGeom()->GetVid(1);
                    vId[cnt1+newpoints-1]     = e->GetGeom()->GetVid(2);

                    // Edge numbering
                    StdRegions::Orientation             edgeOrient;
                    int edge1 = 0;
                    int edge2 = 0;
                    for (int n = 1; n < np-1; n++)
                    {
                        // Edge 0
                        edgeOrient          = e->GetGeom()->GetEorient(0);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            vId[cnt1+n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(0) + n;
                        }
                        else
                        {
                            vId[cnt1+np-1-n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(0) + n;
                        }

                        // Edge 1
                        edgeOrient          = e->GetGeom()->GetEorient(1);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            edge1 += np-n;
                            vId[cnt1+np-1+edge1] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(1) + n;
                        }
                        else
                        {
                            edge1 += n;
                            vId[cnt1+newpoints-1-edge1] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(1) + n;
                        }

                        // Edge 2
                        edgeOrient          = e->GetGeom()->GetEorient(2);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            edge2 += n+1;
                            vId[cnt1+newpoints-1-edge2] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(2) + n;
                        }
                        else
                        {
                            edge2 += np+1-n;
                            vId[cnt1+edge2] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(2) + n;
                        }
                    }

                    // Interior numbering
                    edge2 = 0;
                    for (int n = 1; n < np-1; n++)
                    {
                        edge2 += np+1-n;
                        for (int m = 1; m < np-n-1; m++)
                        {
                            vId[cnt1+edge2+m] = 4*nel + maxN*4*nel + cnt2;
                            cnt2++;
                        }
                    }
                    cnt1+= newpoints;
                }
                break;
            case LibUtilities::eQuadrilateral:
                {
                    int newpoints = LibUtilities::StdQuadData::
                                        getNumberOfCoefficients(np,np);
                    // Vertex numbering
                    vId[cnt1]           = e->GetGeom()->GetVid(0);
                    vId[cnt1+np-1]      = e->GetGeom()->GetVid(1);
                    vId[cnt1+np*np-1]   = e->GetGeom()->GetVid(2);
                    vId[cnt1+np*(np-1)] = e->GetGeom()->GetVid(3);

                    // Edge numbering
                    StdRegions::Orientation             edgeOrient;
                    for (int n = 1; n < np-1; n++)
                    {
                        // Edge 0
                        edgeOrient          = e->GetGeom()->GetEorient(0);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            vId[cnt1+n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(0) + n;
                        }
                        else
                        {
                            vId[cnt1+np-1-n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(0) + n;
                        }

                        // Edge 1
                        edgeOrient          = e->GetGeom()->GetEorient(1);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            vId[cnt1+np-1+n*np] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(1) + n;
                        }
                        else
                        {
                            vId[cnt1+np*np-1-n*np] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(1) + n;
                        }

                        // Edge 2
                        edgeOrient          = e->GetGeom()->GetEorient(2);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            vId[cnt1+np*np-1-n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(2) + n;
                        }
                        else
                        {
                            vId[cnt1+np*(np-1)+n] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(2) + n;
                        }

                        // Edge 3
                        edgeOrient          = e->GetGeom()->GetEorient(3);
                        if (edgeOrient==StdRegions::eForwards)
                        {
                            vId[cnt1+np*(np-1)-n*np] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(3) + n;
                        }
                        else
                        {
                            vId[cnt1+n*np] = 4*nel +
                                            maxN*e->GetGeom()->GetEid(3) + n;
                        }
                    }

                    // Interior numbering
                    for (int n = 1; n < np-1; n++)
                    {
                        for (int m = 1; m < np-1; m++)
                        {
                            vId[cnt1+m*np+n] = 4*nel + maxN*4*nel + cnt2;
                            cnt2++;
                        }
                    }
                    cnt1+= newpoints;
                }
                break;
            default:
                {
                    ASSERTL0(false,"Points not known");
                }
            }
        }

        // Crete new connectivity using homogeneous information
        for(int i = 0; i < nel; ++i)
        {
            int nTris = oldConn[i].size()/3;
            // Create array for new connectivity
            // (3 tetrahedra between each plane for each triangle)
            conn = Array<OneD, int> (4*3*nTris*(nPlanes-1));
            cnt2=0;
            for(int n = 0; n<nTris; n++)
            {
                // Get id of vertices in this triangle (on plane 0)
                int vId0 = vId[oldConn[i][3*n+0]];
                int vId1 = vId[oldConn[i][3*n+1]];
                int vId2 = vId[oldConn[i][3*n+2]];

                // Determine ordering based on global Id of pts
                Array<OneD, int> rot(3);
                if ( (vId0<vId1) && (vId0<vId2))
                {
                    rot[0] = 0;
                    if (vId1<vId2)
                    {
                        rot[1] = 1;
                        rot[2] = 2;
                    }
                    else
                    {
                        rot[1] = 2;
                        rot[2] = 1;
                    }
                }
                else if ( (vId1<vId0) && (vId1<vId2))
                {
                    rot[0] = 1;
                    if (vId0<vId2)
                    {
                        rot[1] = 0;
                        rot[2] = 2;
                    }
                    else
                    {
                        rot[1] = 2;
                        rot[2] = 0;
                    }
                }
                else
                {
                    rot[0] = 2;
                    if (vId0<vId1)
                    {
                        rot[1] = 0;
                        rot[2] = 1;
                    }
                    else
                    {
                        rot[1] = 1;
                        rot[2] = 0;
                    }
                }

                // Define new connectivity with tetrahedra
                for ( int p = 0; p<nPlanes-1; p++)
                {
                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[0]];
                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[1]];
                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[2]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[2]];

                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[0]];
                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[1]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[2]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[1]];

                    conn[cnt2++] = oldConn[i+  p  *nel][3*n+rot[0]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[1]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[2]];
                    conn[cnt2++] = oldConn[i+(p+1)*nel][3*n+rot[0]];
                }
            }
            newConn.push_back(conn);
        }
    }
    else
    {
        ASSERTL0( false, "Points type not supported.");
    }

    m_f->m_fieldPts->SetConnectivity(newConn);
}

void ProcessEquiSpacedOutput::GenOrthoModes(
    int n,
    const Array<OneD, const NekDouble> &phys,
    Array<OneD, NekDouble> &coeffs)
{
    LocalRegions::ExpansionSharedPtr e;
    e = m_f->m_exp[0]->GetExp(n);

    switch (e->DetShapeType())
    {
        case LibUtilities::eTriangle:
        {
            int np0 = e->GetBasis(0)->GetNumPoints();
            int np1 = e->GetBasis(1)->GetNumPoints();
            int np  = max(np0, np1);

            // to ensure points are correctly projected to np need
            // to increase the order slightly of coordinates
            LibUtilities::PointsKey pa(np + 1, e->GetPointsType(0));
            LibUtilities::PointsKey pb(np, e->GetPointsType(1));
            Array<OneD, NekDouble> tophys(np * (np + 1));

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, np, pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B, np, pb);
            StdRegions::StdTriExp OrthoExp(Ba, Bb);

            // interpolate points to new phys points!
            LibUtilities::Interp2D(e->GetBasis(0)->GetBasisKey(),
                                   e->GetBasis(1)->GetBasisKey(), phys, Ba, Bb,
                                   tophys);

            OrthoExp.FwdTrans(tophys, coeffs);
            break;
        }
        case LibUtilities::eQuadrilateral:
        {
            int np0 = e->GetBasis(0)->GetNumPoints();
            int np1 = e->GetBasis(1)->GetNumPoints();
            int np  = max(np0, np1);

            LibUtilities::PointsKey pa(np + 1, e->GetPointsType(0));
            LibUtilities::PointsKey pb(np + 1, e->GetPointsType(1));
            Array<OneD, NekDouble> tophys((np + 1) * (np + 1));

            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A, np, pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_A, np, pb);
            StdRegions::StdQuadExp OrthoExp(Ba, Bb);

            // interpolate points to new phys points!
            LibUtilities::Interp2D(e->GetBasis(0)->GetBasisKey(),
                                   e->GetBasis(1)->GetBasisKey(), phys, Ba, Bb,
                                   tophys);

            OrthoExp.FwdTrans(phys, coeffs);
            break;
        }
        default:
            ASSERTL0(false, "Shape needs setting up");
            break;
    }
}
}
}
