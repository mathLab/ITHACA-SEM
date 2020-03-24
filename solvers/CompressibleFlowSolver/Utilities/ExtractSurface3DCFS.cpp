///////////////////////////////////////////////////////////////////////////////
//
// File ExtractSurface3DCFS.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Extract 2D surfaces from 3D file and output relevant quantities
// for compressible flow solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <MultiRegions/DisContField3D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion.h>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField3D.h>
#include <SpatialDomains/MeshGraph.h>

#include <SolverUtils/SolverUtilsDeclspec.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    string fname = std::string(argv[2]);
    int fdot = fname.find_last_of('.');
    if (fdot != std::string::npos)
    {
        string ending = fname.substr(fdot);

        // If .chk or .fld we exchange the extension in the output file.
        // For all other files (e.g. .bse) we append the extension to avoid
        // conflicts.
        if (ending == ".chk" || ending == ".fld")
        {
            fname = fname.substr(0,fdot);
        }
    }

    fname = fname + ".txt";

    int cnt;
    int id1, id2;
    int i, j, n, e, b;
    Array<OneD, NekDouble> auxArray;

    int nBndEdgePts, nBndEdges, nBndRegions;

    if (argc < 3)
    {
        fprintf(stderr,
                "Usage: ExtractSurface3DCFS meshfile fieldFile\n");
        fprintf(stderr,
                "Extracts a surface from a 3D fld file"
                "(only for CompressibleFlowSolver and purely 3D .fld files)\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(3, argv);
    SpatialDomains::MeshGraphSharedPtr graphShPt =
        SpatialDomains::MeshGraph::Read(vSession);

    std::string                         m_ViscosityType;

    NekDouble                           m_gamma;
    NekDouble                           m_pInf;
    NekDouble                           m_rhoInf;
    NekDouble                           m_uInf;
    NekDouble                           m_vInf;
    NekDouble                           m_wInf;
    NekDouble                           m_gasConstant;
    NekDouble                           m_Twall;
    NekDouble                           m_mu;

    int m_spacedim = 3;
    int nDimensions = m_spacedim;
    int phys_offset;

    // Get gamma parameter from session file.
    ASSERTL0(vSession->DefinesParameter("Gamma"),
             "Compressible flow sessions must define a Gamma parameter.");
    vSession->LoadParameter("Gamma", m_gamma, 1.4);

    // Get E0 parameter from session file.
    ASSERTL0(vSession->DefinesParameter("pInf"),
             "Compressible flow sessions must define a pInf parameter.");
    vSession->LoadParameter("pInf", m_pInf, 101325);

    // Get rhoInf parameter from session file.
    ASSERTL0(vSession->DefinesParameter("rhoInf"),
             "Compressible flow sessions must define a rhoInf parameter.");
    vSession->LoadParameter("rhoInf", m_rhoInf, 1.225);

    // Get uInf parameter from session file.
    ASSERTL0(vSession->DefinesParameter("uInf"),
             "Compressible flow sessions must define a uInf parameter.");
    vSession->LoadParameter("uInf", m_uInf, 0.1);

    // Get vInf parameter from session file.
    if (m_spacedim == 2 || m_spacedim == 3)
    {
        ASSERTL0(vSession->DefinesParameter("vInf"),
                 "Compressible flow sessions must define a vInf parameter"
                 "for 2D/3D problems.");
        vSession->LoadParameter("vInf", m_vInf, 0.0);
    }

    // Get wInf parameter from session file.
    if (m_spacedim == 3)
    {
        ASSERTL0(vSession->DefinesParameter("wInf"),
                 "Compressible flow sessions must define a wInf parameter"
                 "for 3D problems.");
        vSession->LoadParameter("wInf", m_wInf, 0.0);
    }

    vSession->LoadParameter ("GasConstant",   m_gasConstant,   287.058);
    vSession->LoadParameter ("Twall",         m_Twall,         300.15);
    vSession->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
    vSession->LoadParameter ("mu",            m_mu,            1.78e-05);

    //--------------------------------------------------------------------------
    // Import field file
    string                                          fieldFile(argv[2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fieldDef;
    vector<vector<NekDouble> >                      fieldData;

    LibUtilities::Import(fieldFile, fieldDef, fieldData);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Set up Expansion information
    vector< vector<LibUtilities::PointsType> > pointsType;
    for (i = 0; i < fieldDef.size(); ++i)
    {
        vector<LibUtilities::PointsType> ptype;
        for (j = 0; j < 3; ++j)
        {
            ptype.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointsType.push_back(ptype);
    }
    graphShPt->SetExpansions(fieldDef, pointsType);

    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // Define Expansion
    int nfields = fieldDef[0]->m_fields.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields(nfields);

    for(i = 0; i < pFields.size(); i++)
    {
        pFields[i] = MemoryManager<MultiRegions
                                   ::DisContField3D>::AllocateSharedPtr(vSession, graphShPt,
                                                                        vSession->GetVariable(i));
    }

    MultiRegions::ExpList3DSharedPtr Exp3D;
    Exp3D = MemoryManager<MultiRegions::ExpList3D>
        ::AllocateSharedPtr(vSession, graphShPt);

    Exp[0] = Exp3D;

    for (i = 1; i < nfields; ++i)
    {
        Exp[i] = MemoryManager<MultiRegions::ExpList3D>
            ::AllocateSharedPtr(*Exp3D);
    }

    // Count of the point on the surface
    int nSurfacePts = 0;
    if (pFields[0]->GetBndCondExpansions().size())
    {
        nSurfacePts = 0;
        cnt = 0;
        nBndRegions = pFields[0]->GetBndCondExpansions().size();
        for (b = 0; b < nBndRegions; ++b)
        {
            nBndEdges = pFields[0]->GetBndCondExpansions()[b]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = pFields[0]->
                    GetBndCondExpansions()[b]->GetExp(e)->GetTotPoints();

                if (pFields[0]->GetBndConditions()[b]->
                    GetUserDefined() == "WallViscous" ||
                    pFields[0]->GetBndConditions()[b]->
                    GetUserDefined() == "WallAdiabatic" ||
                    pFields[0]->GetBndConditions()[b]->
                    GetUserDefined() == "Wall")
                {
                    nSurfacePts += nBndEdgePts;
                }
            }
        }
    }


    int nSolutionPts = pFields[0]->GetNpoints();
    int nTracePts    = pFields[0]->GetTrace()->GetTotPoints();
    int nElements    = pFields[0]->GetExpSize();

    Array<OneD, NekDouble> tmp(nSolutionPts, 0.0);

    Array<OneD, NekDouble> x(nSolutionPts);
    Array<OneD, NekDouble> y(nSolutionPts);
    Array<OneD, NekDouble> z(nSolutionPts);

    Array<OneD, NekDouble> traceX(nTracePts);
    Array<OneD, NekDouble> traceY(nTracePts);
    Array<OneD, NekDouble> traceZ(nTracePts);

    Array<OneD, NekDouble> surfaceX(nSurfacePts);
    Array<OneD, NekDouble> surfaceY(nSurfacePts);
    Array<OneD, NekDouble> surfaceZ(nSurfacePts);

    pFields[0]->GetCoords(x, y, z);

    pFields[0]->ExtractTracePhys(x, traceX);
    pFields[0]->ExtractTracePhys(y, traceY);
    pFields[0]->ExtractTracePhys(z, traceZ);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Copy data from field file
    Array<OneD, Array<OneD, NekDouble> > uFields(nfields);
    Array<OneD, Array<OneD, NekDouble> > traceFields(nfields);
    Array<OneD, Array<OneD, NekDouble> > surfaceFields(nfields);

    // Extract the physical values of the solution at the boundaries
    for (j = 0; j < nfields; ++j)
    {
        uFields[j]       = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        traceFields[j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        surfaceFields[j] = Array<OneD, NekDouble>(nSurfacePts, 0.0);


        for (i = 0; i < fieldData.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fieldDef[i], fieldData[i],
                                        fieldDef[i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(), Exp[j]->UpdatePhys());
        Vmath::Vcopy(nSolutionPts, Exp[j]->GetPhys(), 1, uFields[j], 1);
        pFields[0]->ExtractTracePhys(uFields[j], traceFields[j]);
    }


    //Fields to add in the output file
    int nfieldsAdded = 34;
    Array<OneD, Array<OneD, NekDouble> > traceFieldsAdded(nfieldsAdded);
    Array<OneD, Array<OneD, NekDouble> > surfaceFieldsAdded(nfieldsAdded);

    for (j = 0; j < nfieldsAdded; ++j)
    {
        traceFieldsAdded[j] = Array<OneD, NekDouble>(nTracePts, 0.0);
        surfaceFieldsAdded[j] = Array<OneD, NekDouble>(nSurfacePts, 0.0);
    }

    /******** Evaluation of normals and tangents on the trace *****************
     * nx -> traceFieldsAdded[0];
     * ny -> traceFieldsAdded[1];
     * nz -> traceFieldsAdded[2];
     * bx -> traceFieldsAdded[3];
     * by -> traceFieldsAdded[4];
     * bz -> traceFieldsAdded[5];
     * tx -> traceFieldsAdded[6];
     * ty -> traceFieldsAdded[7];
     * tz -> traceFieldsAdded[8];

     ***************************************************************************/

    Array<OneD, Array<OneD, NekDouble> > m_traceNormals (nDimensions);
    for(i = 0; i < nDimensions; ++i)
    {
        m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
    }
    pFields[0]->GetTrace()->GetNormals(m_traceNormals);

    Array<OneD, Array<OneD, NekDouble> > m_traceTangents (nDimensions);
    Array<OneD, Array<OneD, NekDouble> > m_traceBinormals (nDimensions);
    Array<OneD, Array<OneD, NekDouble> > h (nDimensions);
    Array<OneD, NekDouble > tmpNorm (nTracePts, 1.0);
    Array<OneD, NekDouble > NormH (nTracePts, 0.0);
    Array<OneD, NekDouble > tmpTrace (nTracePts, 0.0);


    for(i = 0; i < nDimensions; ++i)
    {
        m_traceTangents[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
        m_traceBinormals[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
        h[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
    }

    // Normals

    // nx
    Vmath::Vcopy(nTracePts,
                 &m_traceNormals[0][0], 1,
                 &traceFieldsAdded[0][0], 1);

    // ny
    Vmath::Vcopy(nTracePts,
                 &m_traceNormals[1][0], 1,
                 &traceFieldsAdded[1][0], 1);

    // nz
    Vmath::Vcopy(nTracePts,
                 &m_traceNormals[2][0], 1,
                 &traceFieldsAdded[2][0], 1);


    // Tangents and Binormals
    // h1
    Vmath::Vadd(nTracePts,
                &m_traceNormals[0][0], 1,
                &tmpNorm[0], 1,
                &h[0][0], 1);
    // h2
    Vmath::Vcopy(nTracePts,
                 &m_traceNormals[1][0], 1,
                 &h[1][0], 1);
    // h3
    Vmath::Vcopy(nTracePts,
                 &m_traceNormals[2][0], 1,
                 &h[2][0], 1);

    // Norm of h
    for (i = 0; i < m_spacedim; i++)
    {
        Vmath::Vvtvp (nTracePts, &h[i][0], 1, &h[i][0], 1,
                      &NormH[0],1, &NormH[0],1);
    }

    //b1
    Vmath::Vmul(nTracePts,
                &h[0][0], 1,
                &h[1][0], 1,
                &tmpTrace[0],1);

    Vmath::Vdiv(nTracePts,
                &tmpTrace[0],1,
                &NormH[0], 1,
                &tmpTrace[0],1);

    Vmath::Smul(nTracePts, -2.0,
                &tmpTrace[0], 1,
                &m_traceBinormals[0][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceBinormals[0][0], 1,
                 &traceFieldsAdded[3][0], 1);


    //b2
    Vmath::Vmul(nTracePts,
                &h[1][0], 1,
                &h[1][0], 1,
                &tmpTrace[0],1);

    Vmath::Vdiv(nTracePts,
                &tmpTrace[0],1,
                &NormH[0], 1,
                &tmpTrace[0],1);

    Vmath::Smul(nTracePts, -2.0,
                &tmpTrace[0], 1,
                &tmpTrace[0], 1);

    Vmath::Vadd(nTracePts,
                &tmpTrace[0], 1,
                &tmpNorm[0], 1,
                &m_traceBinormals[1][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceBinormals[1][0], 1,
                 &traceFieldsAdded[4][0], 1);


    //b3
    Vmath::Vmul(nTracePts,
                &h[1][0], 1,
                &h[2][0], 1,
                &tmpTrace[0],1);

    Vmath::Vdiv(nTracePts,
                &tmpTrace[0],1,
                &NormH[0], 1,
                &tmpTrace[0],1);

    Vmath::Smul(nTracePts, -2.0,
                &tmpTrace[0], 1,
                &m_traceBinormals[2][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceBinormals[2][0], 1,
                 &traceFieldsAdded[5][0], 1);


    //t1
    Vmath::Vmul(nTracePts,
                &h[0][0], 1,
                &h[2][0], 1,
                &tmpTrace[0],1);

    Vmath::Vdiv(nTracePts,
                &tmpTrace[0],1,
                &NormH[0], 1,
                &tmpTrace[0],1);

    Vmath::Smul(nTracePts, -2.0,
                &tmpTrace[0], 1,
                &m_traceTangents[0][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceTangents[0][0], 1,
                 &traceFieldsAdded[6][0], 1);

    //t2
    Vmath::Vcopy(nTracePts,
                 &m_traceBinormals[2][0], 1,
                 &m_traceTangents[1][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceTangents[1][0], 1,
                 &traceFieldsAdded[7][0], 1);


    //t3
    Vmath::Vmul(nTracePts,
                &h[2][0], 1,
                &h[2][0], 1,
                &tmpTrace[0],1);

    Vmath::Vdiv(nTracePts,
                &tmpTrace[0],1,
                &NormH[0], 1,
                &tmpTrace[0],1);

    Vmath::Smul(nTracePts, -2.0,
                &tmpTrace[0], 1,
                &tmpTrace[0], 1);

    Vmath::Vadd(nTracePts,
                &tmpTrace[0], 1,
                &tmpNorm[0], 1,
                &m_traceTangents[2][0], 1);

    Vmath::Vcopy(nTracePts,
                 &m_traceTangents[2][0], 1,
                 &traceFieldsAdded[8][0], 1);


    /******** Evaluation of the pressure ***************************************
     * P    = (E-1/2.*rho.*((rhou./rho).^2+(rhov./rho).^2))*(gamma - 1);
     * P -> traceFieldsAdded[9];
     ***************************************************************************/

    Array<OneD, NekDouble> pressure(nSolutionPts, 0.0);
    NekDouble gammaMinusOne    = m_gamma - 1.0;


    for (i = 0; i < m_spacedim; i++)
    {
        Vmath::Vmul(nSolutionPts,
                    &uFields[i + 1][0], 1,
                    &uFields[i + 1][0], 1,
                    &tmp[0],1);


        Vmath::Smul(nSolutionPts, 0.5,
                    &tmp[0], 1,
                    &tmp[0], 1);

        Vmath::Vadd(nSolutionPts,
                    &pressure[0], 1,
                    &tmp[0], 1,
                    &pressure[0], 1);
    }

    Vmath::Vdiv(nSolutionPts,
                &pressure[0], 1,
                &uFields[0][0], 1,
                &pressure[0],1);

    Vmath::Vsub(nSolutionPts,
                &uFields[nfields - 1][0], 1,
                &pressure[0], 1,
                &pressure[0],1);

    Vmath::Smul(nSolutionPts, gammaMinusOne,
                &pressure[0], 1,
                &pressure[0], 1);

    // Extract trace
    pFields[0]->ExtractTracePhys(pressure, traceFieldsAdded[9]);

    /******** Evaluation of the temperature ************************************
     * T = P/(R*rho);
     * T -> traceFieldsAdded[10];
     ***************************************************************************/

    Array<OneD, NekDouble> temperature(nSolutionPts, 0.0);

    Vmath::Vdiv(nSolutionPts,
                &pressure[0], 1,
                &uFields[0][0], 1,
                &temperature[0],1);

    NekDouble GasConstantInv =  1.0/m_gasConstant;
    Vmath::Smul(nSolutionPts, GasConstantInv,
                &temperature[0], 1,
                &temperature[0], 1);

    // Extract trace
    pFields[0]->ExtractTracePhys(temperature, traceFieldsAdded[10]);

    /*** Evaluation of the temperature gradient in the normal direction ********
     * DT_n -> traceFieldsAdded[11]
     ***************************************************************************/

    Array<OneD, Array<OneD, NekDouble> > Dtemperature(nDimensions);
    Array<OneD, Array<OneD, NekDouble> > traceDtemperature(nDimensions);

    for (i = 0; i < nDimensions; ++ i)
    {
        Dtemperature[i]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        traceDtemperature[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    for (i = 0; i < nDimensions; ++ i)
    {
        for (n = 0; n < nElements; n++)
        {
            phys_offset = pFields[0]->GetPhys_Offset(n);

            pFields[i]->GetExp(n)->PhysDeriv(
                                             i, temperature + phys_offset,
                                             auxArray = Dtemperature[i] + phys_offset);
        }
        // Extract trace
        pFields[0]->ExtractTracePhys(Dtemperature[i], traceDtemperature[i]);
    }

    for(i = 0; i < nDimensions; ++i)
    {
        Vmath::Vmul(nTracePts,
                    &m_traceNormals[i][0], 1,
                    &traceDtemperature[i][0], 1,
                    &tmp[0],1);

        Vmath::Vadd(nTracePts,
                    &traceFieldsAdded[11][0], 1,
                    &tmp[0], 1,
                    &traceFieldsAdded[11][0], 1);
    }

    /*** Evaluation of the pressure gradient ***********************************
     * DP_t -> traceFieldsAdded[12]   tangent direction
     * DP_b -> traceFieldsAdded[13]   binormal direction
     * DP_x -> traceFieldsAdded[14]
     * DP_y -> traceFieldsAdded[15]
     * DP_z -> traceFieldsAdded[16]
     ***************************************************************************/

    Array<OneD, Array<OneD, NekDouble> > Dpressure(nDimensions);
    Array<OneD, Array<OneD, NekDouble> > traceDpressure(nDimensions);

    for (i = 0; i < nDimensions; ++ i)
    {
        Dpressure[i]  = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        traceDpressure[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    for (i = 0; i < nDimensions; ++ i)
    {
        for (n = 0; n < nElements; n++)
        {
            phys_offset = pFields[0]->GetPhys_Offset(n);

            pFields[i]->GetExp(n)->PhysDeriv(
                                             i, pressure + phys_offset,
                                             auxArray = Dpressure[i] + phys_offset);
        }
        // Extract trace
        pFields[0]->ExtractTracePhys(Dpressure[i], traceDpressure[i]);
    }

    // Dp_t
    for(i = 0; i < nDimensions; ++i)
    {
        Vmath::Vmul(nTracePts,
                    &m_traceTangents[i][0], 1,
                    &traceDpressure[i][0], 1,
                    &tmp[0],1);

        Vmath::Vadd(nTracePts,
                    &traceFieldsAdded[12][0], 1,
                    &tmp[0], 1,
                    &traceFieldsAdded[12][0], 1);
    }

    // Dp_b
    for(i = 0; i < nDimensions; ++i)
    {
        Vmath::Vmul(nTracePts,
                    &m_traceBinormals[i][0], 1,
                    &traceDpressure[i][0], 1,
                    &tmp[0],1);

        Vmath::Vadd(nTracePts,
                    &traceFieldsAdded[13][0], 1,
                    &tmp[0], 1,
                    &traceFieldsAdded[13][0], 1);
    }


    // Dp_x
    Vmath::Vcopy(nTracePts,
                 &traceDpressure[0][0], 1,
                 &traceFieldsAdded[14][0], 1);

    // Dp_y
    Vmath::Vcopy(nTracePts,
                 &traceDpressure[1][0], 1,
                 &traceFieldsAdded[15][0], 1);

    // Dp_z
    Vmath::Vcopy(nTracePts,
                 &traceDpressure[2][0], 1,
                 &traceFieldsAdded[16][0], 1);


    /** Evaluation of the velocity gradient in the cartesian directions
     * Du_x:    traceFieldsAdded[17]
     * Du_y:    traceFieldsAdded[18]
     * Du_z:    traceFieldsAdded[19]
     * Dv_x:    traceFieldsAdded[20]
     * Dv_y:    traceFieldsAdded[21]
     * Dv_z:    traceFieldsAdded[22]
     * Dw_x:    traceFieldsAdded[23]
     * Dw_y:    traceFieldsAdded[24]
     * Dw_z:    traceFieldsAdded[25]
     **/

    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > Dvelocity(nDimensions);
    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceDvelocity(nDimensions);
    Array<OneD, Array<OneD, NekDouble> > velocity(nDimensions);

    for (i = 0; i < nDimensions; ++ i)
    {
        Dvelocity[i]      = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
        traceDvelocity[i] = Array<OneD, Array<OneD, NekDouble> >(nDimensions);
        velocity[i]       = Array<OneD, NekDouble>(nSolutionPts, 0.0);

        Vmath::Vdiv(nSolutionPts, uFields[i+1], 1, uFields[0], 1,
                    velocity[i], 1);

        for (j = 0; j < nDimensions; ++j)
        {
            Dvelocity[i][j]      = Array<OneD, NekDouble>(nSolutionPts, 0.0);
            traceDvelocity[i][j] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }
    }

    for (i = 0; i < nDimensions; ++i)
    {
        for (j = 0; j < nDimensions; ++j)
        {
            for (n = 0; n < nElements; n++)
            {
                phys_offset = pFields[0]->GetPhys_Offset(n);

                pFields[i]->GetExp(n)->PhysDeriv(
                                                 j, velocity[i] + phys_offset,
                                                 auxArray = Dvelocity[i][j] + phys_offset);
            }

            // Extract trace
            pFields[0]->ExtractTracePhys(Dvelocity[i][j], traceDvelocity[i][j]);
        }
    }

    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[0][0][0], 1,
                 &traceFieldsAdded[17][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[0][1][0], 1,
                 &traceFieldsAdded[18][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[0][2][0], 1,
                 &traceFieldsAdded[19][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[1][0][0], 1,
                 &traceFieldsAdded[20][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[1][1][0], 1,
                 &traceFieldsAdded[21][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[1][2][0], 1,
                 &traceFieldsAdded[22][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[2][0][0], 1,
                 &traceFieldsAdded[23][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[2][1][0], 1,
                 &traceFieldsAdded[24][0], 1);
    Vmath::Vcopy(nTracePts,
                 &traceDvelocity[2][2][0], 1,
                 &traceFieldsAdded[25][0], 1);


    /*** Evaluation of shear stresses ******************************************
     * tau_xx -> traceFieldsAdded[26]
     * tau_yy -> traceFieldsAdded[27]
     * tau_zz -> traceFieldsAdded[28]
     * tau_xy -> traceFieldsAdded[29]
     * tau_xz -> traceFieldsAdded[30]
     * tau_yz -> traceFieldsAdded[31]
     ***************************************************************************/

    // Stokes hypotesis
    const NekDouble lambda = -2.0/3.0;

    // Auxiliary variables
    Array<OneD, NekDouble > mu    (nSolutionPts, 0.0);
    Array<OneD, NekDouble > mu2   (nSolutionPts, 0.0);
    Array<OneD, NekDouble > divVel(nSolutionPts, 0.0);

    // Variable viscosity through the Sutherland's law
    if (m_ViscosityType == "Variable")
    {
        NekDouble mu_star = m_mu;
        NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
        NekDouble ratio;

        for (int i = 0; i < nSolutionPts; ++i)
        {
            ratio = temperature[i] / T_star;
            mu[i] = mu_star * ratio * sqrt(ratio) *
                (T_star + 110.0) / (temperature[i] + 110.0);
        }
    }
    else
    {
        Vmath::Fill(nSolutionPts, m_mu, &mu[0], 1);
    }

    // Computing diagonal terms of viscous stress tensor
    Array<OneD, Array<OneD, NekDouble> > temp(m_spacedim);
    Array<OneD, Array<OneD, NekDouble> > Sgg(m_spacedim);

    // mu2 = 2 * mu
    Vmath::Smul(nSolutionPts, 2.0, &mu[0], 1, &mu2[0], 1);

    // Velocity divergence
    Vmath::Vadd(nSolutionPts, &divVel[0], 1,
                &Dvelocity[0][0][0], 1, &divVel[0], 1);
    Vmath::Vadd(nSolutionPts, &divVel[0], 1,
                &Dvelocity[1][1][0], 1, &divVel[0], 1);

    // Velocity divergence scaled by lambda * mu
    Vmath::Smul(nSolutionPts, lambda, &divVel[0], 1, &divVel[0], 1);
    Vmath::Vmul(nSolutionPts, &mu[0], 1, &divVel[0], 1, &divVel[0], 1);

    // Diagonal terms of viscous stress tensor (Sxx, Syy)
    // Sjj = 2 * mu * du_j/dx_j - (2 / 3) * mu * sum_j(du_j/dx_j)
    for (j = 0; j < m_spacedim; ++j)
    {
        temp[j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);
        Sgg[j] = Array<OneD, NekDouble>(nSolutionPts, 0.0);

        Vmath::Vmul(nSolutionPts, &mu2[0], 1, &Dvelocity[j][j][0], 1,
                    &temp[j][0], 1);

        Vmath::Vadd(nSolutionPts, &temp[j][0], 1, &divVel[0], 1, &Sgg[j][0], 1);
    }

    // Extra diagonal terms of viscous stress tensor
    Array<OneD, NekDouble > Sxy(nSolutionPts, 0.0);
    Array<OneD, NekDouble > Sxz(nSolutionPts, 0.0);
    Array<OneD, NekDouble > Syz(nSolutionPts, 0.0);

    // Sxy = (du/dy + dv/dx)
    Vmath::Vadd(nSolutionPts, &Dvelocity[0][1][0], 1,
                &Dvelocity[1][0][0], 1, &Sxy[0], 1);

    // Sxz = (du/dz + dw/dx)
    Vmath::Vadd(nSolutionPts, &Dvelocity[0][2][0], 1,
                &Dvelocity[2][0][0], 1, &Sxz[0], 1);

    // Syz = (dv/dz + dw/dy)
    Vmath::Vadd(nSolutionPts, &Dvelocity[1][2][0], 1,
                &Dvelocity[2][1][0], 1, &Syz[0], 1);

    // Sxy = mu * (du/dy + dv/dx)
    Vmath::Vmul(nSolutionPts, &mu[0], 1, &Sxy[0], 1, &Sxy[0], 1);

    // Sxz = mu * (du/dy + dv/dx)
    Vmath::Vmul(nSolutionPts, &mu[0], 1, &Sxz[0], 1, &Sxz[0], 1);

    // Syz = mu * (du/dy + dv/dx)
    Vmath::Vmul(nSolutionPts, &mu[0], 1, &Syz[0], 1, &Syz[0], 1);



    pFields[0]->ExtractTracePhys(Sgg[0], traceFieldsAdded[26]);
    pFields[0]->ExtractTracePhys(Sgg[1], traceFieldsAdded[27]);
    pFields[0]->ExtractTracePhys(Sgg[2], traceFieldsAdded[28]);
    pFields[0]->ExtractTracePhys(Sxy, traceFieldsAdded[29]);
    pFields[0]->ExtractTracePhys(Sxz, traceFieldsAdded[30]);
    pFields[0]->ExtractTracePhys(Syz, traceFieldsAdded[31]);

    /*** Evaluation of dinamic viscosity ***************************************
     * mu -> traceFieldsAdded[32]
     ***************************************************************************/

    pFields[0]->ExtractTracePhys(mu, traceFieldsAdded[32]);

    /*** Evaluation of Mach number *********************************************
     * M -> traceFieldsAdded[33]
     ***************************************************************************/
    NekDouble gamma    = m_gamma;

    // Speed of sound
    Array<OneD, NekDouble> soundspeed(nSolutionPts, 0.0);

    Vmath::Vdiv (nSolutionPts, pressure, 1, uFields[0], 1, soundspeed, 1);
    Vmath::Smul (nSolutionPts, gamma, soundspeed, 1, soundspeed, 1);
    Vmath::Vsqrt(nSolutionPts, soundspeed, 1, soundspeed, 1);

    // Mach
    Array<OneD, NekDouble> mach(nSolutionPts, 0.0);

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vvtvp(nSolutionPts,
                     uFields[i + 1], 1,
                     uFields[i + 1], 1,
                     mach, 1, mach, 1);
    }

    Vmath::Vdiv(nSolutionPts,  mach, 1, uFields[0], 1, mach, 1);
    Vmath::Vdiv(nSolutionPts,  mach, 1, uFields[0], 1, mach, 1);
    Vmath::Vsqrt(nSolutionPts, mach, 1, mach, 1);
    Vmath::Vdiv(nSolutionPts,  mach, 1, soundspeed, 1, mach, 1);

    pFields[0]->ExtractTracePhys(mach, traceFieldsAdded[33]);



    /**************************************************************************/
    // Extract coordinates

    if (pFields[0]->GetBndCondExpansions().size())
    {
        id1 = 0;
        cnt = 0;
        nBndRegions = pFields[0]->GetBndCondExpansions().size();
        for (b = 0; b < nBndRegions; ++b)
        {
            nBndEdges = pFields[0]->GetBndCondExpansions()[b]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = pFields[0]->
                    GetBndCondExpansions()[b]->GetExp(e)->GetTotPoints();

                id2 = pFields[0]->GetTrace()->
                    GetPhys_Offset(pFields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));

                if (pFields[0]->GetBndConditions()[b]->
                        GetUserDefined() == "WallViscous" ||
                    pFields[0]->GetBndConditions()[b]->
                        GetUserDefined() == "WallAdiabatic" ||
                    pFields[0]->GetBndConditions()[b]->
                        GetUserDefined() == "Wall")
                {

                    Vmath::Vcopy(nBndEdgePts, &traceX[id2], 1,
                                 &surfaceX[id1], 1);

                    Vmath::Vcopy(nBndEdgePts, &traceY[id2], 1,
                                 &surfaceY[id1], 1);

                    Vmath::Vcopy(nBndEdgePts, &traceZ[id2], 1,
                                 &surfaceZ[id1], 1);

                    id1 += nBndEdgePts;
                }
            }
        }
    }

    // Extract fields
    if (pFields[0]->GetBndCondExpansions().size())
    {

        for (j = 0; j < nfields; ++j)
        {
            cout << "field " << j << endl;

            id1 = 0;
            cnt = 0;
            nBndRegions = pFields[j]->GetBndCondExpansions().size();
            for (b = 0; b < nBndRegions; ++b)
            {
                nBndEdges = pFields[j]->GetBndCondExpansions()[b]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = pFields[j]->
                        GetBndCondExpansions()[b]->GetExp(e)->GetTotPoints();

                    id2 = pFields[j]->GetTrace()->
                        GetPhys_Offset(pFields[j]->GetTraceMap()->
                                       GetBndCondIDToGlobalTraceID(cnt++));

                    if (pFields[j]->GetBndConditions()[b]->
                            GetUserDefined() == "WallViscous" ||
                        pFields[j]->GetBndConditions()[b]->
                            GetUserDefined() == "WallAdiabatic" ||
                        pFields[j]->GetBndConditions()[b]->
                            GetUserDefined() == "Wall")
                    {
                        Vmath::Vcopy(nBndEdgePts, &traceFields[j][id2], 1,
                                     &surfaceFields[j][id1], 1);

                        id1 += nBndEdgePts;
                    }
                }
            }
        }
    }

    // Extract fields added
    if (pFields[0]->GetBndCondExpansions().size())
    {
        for (j = 0; j < nfieldsAdded; ++j)
        {
            cout << "field added " << j << endl;

            id1 = 0;
            cnt = 0;
            nBndRegions = pFields[0]->GetBndCondExpansions().size();
            for (b = 0; b < nBndRegions; ++b)
            {
                nBndEdges = pFields[0]->GetBndCondExpansions()[b]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = pFields[0]->
                        GetBndCondExpansions()[b]->GetExp(e)->GetTotPoints();

                    id2 = pFields[0]->GetTrace()->
                        GetPhys_Offset(pFields[0]->GetTraceMap()->
                                       GetBndCondIDToGlobalTraceID(cnt++));

                    if (pFields[0]->GetBndConditions()[b]->
                            GetUserDefined() == "WallViscous" ||
                        pFields[0]->GetBndConditions()[b]->
                            GetUserDefined() == "WallAdiabatic" ||
                        pFields[0]->GetBndConditions()[b]->
                            GetUserDefined() == "Wall")
                    {
                        Vmath::Vcopy(nBndEdgePts, &traceFieldsAdded[j][id2], 1,
                                     &surfaceFieldsAdded[j][id1], 1);

                        id1 += nBndEdgePts;
                    }
                }
            }
        }
    }

    //==========================================================================
    //==========================================================================
    //==========================================================================

    // Print the surface coordinates and the surface solution in a .txt file
    ofstream outfile;
    outfile.open(fname.c_str());
    outfile <<  "%  x[m] " << " \t"
            << "y[m] " << " \t"
            << "z[m] " << " \t"
            << "nx[]  " << " \t"
            << "ny[]  " << " \t"
            << "nz[]  " << " \t"
            << "bx[]  " << " \t"
            << "by[]  " << " \t"
            << "bz[]  " << " \t"
            << "tx[]  " << " \t"
            << "ty[]  " << " \t"
            << "tz[]  " << " \t"
            << "rho[kg/m^3] " << " \t"
            << "rhou[kg/(m^2 s)] " << " \t"
            << "rhov[kg/(m^2 s)] " << " \t"
            << "rhow[kg/(m^2 s)] " << " \t"
            << "E[Pa] " << " \t"
            << "p[Pa] " << " \t"
            << "T[k]  " << " \t"
            << "dT/dn[k/m]  "  << " \t"
            << "dp/dT[Pa/m]  " << " \t"
            << "dp/dB[Pa/m]  " << " \t"
            << "dp/dx[Pa/m]  " << " \t"
            << "dp/dy[Pa/m]  " << " \t"
            << "dp/dz[Pa/m]  " << " \t"
            << "du/dx[s^-1]  " << " \t"
            << "du/dy[s^-1]  " << " \t"
            << "du/dz[s^-1]  " << " \t"
            << "dv/dx[s^-1]  " << " \t"
            << "dv/dy[s^-1]  " << " \t"
            << "dv/dz[s^-1]  " << " \t"
            << "dw/dx[s^-1]  " << " \t"
            << "dw/dy[s^-1]  " << " \t"
            << "dw/dz[s^-1]  " << " \t"
            << "tau_xx[Pa]   " << " \t"
            << "tau_yy[Pa]   " << " \t"
            << "tau_zz[Pa]   " << " \t"
            << "tau_xy[Pa]   " << " \t"
            << "tau_xz[Pa]   " << " \t"
            << "tau_yz[Pa]   " << " \t"
            << "mu[Pa s]     " << " \t"
            << "M[] " << " \t"
            << endl;
    for (i = 0; i < nSurfacePts; ++i)
    {
        outfile << scientific
                << setw (17)
                << setprecision(16)
                << surfaceX[i] << " \t "
                << surfaceY[i] << " \t "
                << surfaceZ[i] << " \t "
                << surfaceFieldsAdded[0][i] << " \t "
                << surfaceFieldsAdded[1][i] << " \t "
                << surfaceFieldsAdded[2][i] << " \t "
                << surfaceFieldsAdded[3][i] << " \t "
                << surfaceFieldsAdded[4][i] << " \t "
                << surfaceFieldsAdded[5][i] << " \t "
                << surfaceFieldsAdded[6][i] << " \t "
                << surfaceFieldsAdded[7][i] << " \t "
                << surfaceFieldsAdded[8][i] << " \t "
                << surfaceFields[0][i] << " \t "
                << surfaceFields[1][i] << " \t "
                << surfaceFields[2][i] << " \t "
                << surfaceFields[3][i] << " \t "
                << surfaceFields[4][i] << " \t "
                << surfaceFieldsAdded[9][i] << " \t "
                << surfaceFieldsAdded[10][i] << " \t "
                << surfaceFieldsAdded[11][i] << " \t "
                << surfaceFieldsAdded[12][i] << " \t "
                << surfaceFieldsAdded[13][i] << " \t "
                << surfaceFieldsAdded[14][i] << " \t "
                << surfaceFieldsAdded[15][i] << " \t "
                << surfaceFieldsAdded[16][i] << " \t "
                << surfaceFieldsAdded[17][i] << " \t "
                << surfaceFieldsAdded[18][i] << " \t "
                << surfaceFieldsAdded[19][i] << " \t "
                << surfaceFieldsAdded[20][i] << " \t "
                << surfaceFieldsAdded[21][i] << " \t "
                << surfaceFieldsAdded[22][i] << " \t "
                << surfaceFieldsAdded[23][i] << " \t "
                << surfaceFieldsAdded[24][i] << " \t "
                << surfaceFieldsAdded[25][i] << " \t "
                << surfaceFieldsAdded[26][i] << " \t "
                << surfaceFieldsAdded[27][i] << " \t "
                << surfaceFieldsAdded[28][i] << " \t "
                << surfaceFieldsAdded[29][i] << " \t "
                << surfaceFieldsAdded[30][i] << " \t "
                << surfaceFieldsAdded[31][i] << " \t "
                << surfaceFieldsAdded[32][i] << " \t "
                << surfaceFieldsAdded[33][i] << " \t "
                << endl;
    }
    outfile << endl << endl;
    outfile.close();

    return 0;
}
