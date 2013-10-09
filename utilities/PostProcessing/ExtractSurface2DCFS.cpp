#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <MultiRegions/DisContField2D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion.h>

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <SolverUtils/SolverUtilsDeclspec.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int cnt;
    int id1, id2;
    int i, j, e, b;
    
    int nBndEdgePts, nBndEdges, nBndRegions;
        
    if (argc < 3)
    {
        fprintf(stderr,
                "Usage: ExtractSurface2DCFS meshfile fieldFile\n");
        fprintf(stderr,
                "Extracts a surface from a 2D fld file" 
                "(only for CompressibleFlowSolver and purely 2D .fld files)\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(3, argv);

    //--------------------------------------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = 
        SpatialDomains::MeshGraph::Read(vSession);
    //--------------------------------------------------------------------------

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
        for (j = 0; j < 2; ++j)
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
        
    for(i = 0; i < pFields.num_elements(); i++)
    {
        pFields[i] = MemoryManager<MultiRegions
        ::DisContField2D>::AllocateSharedPtr(vSession, graphShPt, 
                                             vSession->GetVariable(i));
    }
    
    MultiRegions::ExpList2DSharedPtr Exp2D;
    Exp2D = MemoryManager<MultiRegions::ExpList2D>
        ::AllocateSharedPtr(vSession, graphShPt);
    
    Exp[0] = Exp2D;

    for (i = 1; i < nfields; ++i)
    {
        Exp[i] = MemoryManager<MultiRegions::ExpList2D>
            ::AllocateSharedPtr(*Exp2D);
    }
    
    int nSolutionPts = pFields[0]->GetNpoints();
    int nTracePts    = pFields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble> x(nSolutionPts);
    Array<OneD, NekDouble> y(nSolutionPts); 
    Array<OneD, NekDouble> z(nSolutionPts);
    
    Array<OneD, NekDouble> traceX(nTracePts);
    Array<OneD, NekDouble> traceY(nTracePts); 
    Array<OneD, NekDouble> traceZ(nTracePts);
    
    Array<OneD, NekDouble> surfaceX(nTracePts);
    Array<OneD, NekDouble> surfaceY(nTracePts); 
    Array<OneD, NekDouble> surfaceZ(nTracePts);

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
        surfaceFields[j] = Array<OneD, NekDouble>(nTracePts, 0.0);


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
    //--------------------------------------------------------------------------
    
    if (pFields[0]->GetBndCondExpansions().num_elements())
    {
        id1 = 0;
        cnt = 0;
        nBndRegions = pFields[0]->GetBndCondExpansions().num_elements();
        for (b = 0; b < nBndRegions; ++b)
        {
            nBndEdges = pFields[0]->GetBndCondExpansions()[b]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = pFields[0]->
                GetBndCondExpansions()[b]->GetExp(e)->GetNumPoints(0);
                    
                id2 = pFields[0]->GetTrace()->
                GetPhys_Offset(pFields[0]->GetTraceMap()->
                    GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                if (pFields[0]->GetBndConditions()[b]->
                    GetUserDefined() == SpatialDomains::eWallViscous || 
                    pFields[0]->GetBndConditions()[b]->
                    GetUserDefined() == SpatialDomains::eWall)
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
    
    if (pFields[0]->GetBndCondExpansions().num_elements())
    {
        for (j = 0; j < nfields; ++j)
        {
            id1 = 0;
            cnt = 0;
            nBndRegions = pFields[j]->GetBndCondExpansions().num_elements();
            for (b = 0; b < nBndRegions; ++b)
            {
                nBndEdges = pFields[j]->GetBndCondExpansions()[b]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = pFields[j]->
                    GetBndCondExpansions()[b]->GetExp(e)->GetNumPoints(0);
                                        
                    id2 = pFields[j]->GetTrace()->
                    GetPhys_Offset(pFields[j]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    if (pFields[j]->GetBndConditions()[b]->
                        GetUserDefined() == SpatialDomains::eWallViscous || 
                        pFields[j]->GetBndConditions()[b]->
                        GetUserDefined() == SpatialDomains::eWall)
                    {
                        Vmath::Vcopy(nBndEdgePts, &traceFields[j][id2], 1,
                                     &surfaceFields[j][id1], 1);
                                                
                        id1 += nBndEdgePts;
                    }
                }
            }
        }
    }
    
    // Print the surface coordinates and the surface solution in a .txt file
    ofstream outfile;
    outfile.open("surfaceQuantities.txt");
    for (i = 0; i < id1; ++i)
    {
        outfile << scientific 
        << setw (17) 
        << setprecision(16) 
        << surfaceX[i] << " \t " 
        << surfaceY[i] << " \t " 
        << surfaceZ[i] << " \t " 
        << surfaceFields[0][i] << " \t " 
        << surfaceFields[1][i] << " \t " 
        << surfaceFields[2][i] << " \t "
        << surfaceFields[3][i] << " \t "
        << endl;
    }
    outfile << endl << endl;
    outfile.close();
    
    /*
    //--------------------------------------------------------------------------
    // Copy data from field file
    
    // Extract the physical values of the solution at the boundaries
    for (j = 0; j < nfields; ++j)
    {
        cnt = 0;
        nBndRegions = Exp[j]->GetBndCondExpansions().num_elements();
        for (b = 0; b < nBndRegions; ++b)
        {
            nBndEdges = Exp[j]->GetBndCondExpansions()[b]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = Exp[j]->
                    GetBndCondExpansions()[b]->GetExp(e)->GetNumPoints(0);
                
                id1 = Exp[j]->GetBndCondExpansions()[b]->GetPhys_Offset(e);
                
                id2 = Exp[0]->GetTrace()->
                    GetPhys_Offset(Exp[0]->GetTraceMap()->
                        GetBndCondTraceToGlobalTraceMap(cnt++));
                
                for (i = 0; i < fieldData.size(); ++i)
                {
                    if (Exp[j]->GetBndConditions()[b]->
                            GetUserDefined() == SpatialDomains::eWallViscous || 
                        Exp[j]->GetBndConditions()[b]->
                            GetUserDefined() == SpatialDomains::eWall)
                    {
                        Exp[j]->ExtractDataToCoeffs(fieldDef[i], fieldData[i],
                                                    fieldDef[i]->m_fields[j],
                                                    Exp[j]->UpdateCoeffs());
                    }
                }
            }
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(), Exp[j]->UpdatePhys());
    }
    //--------------------------------------------------------------------------
*/

    /*    
    //----------------------------------------------
    // Probe data fields    
    Array<OneD, Array<OneD, NekDouble> > gloCoord();

    for (int i = 0; i < N; ++i)
    {
        gloCoord[0] = x0 + i*dx;
        gloCoord[1] = y0 + i*dy;
        gloCoord[2] = z0 + i*dz;
        cout << gloCoord[0] << "   " << gloCoord[1] << "   " << gloCoord[2];
        
        int ExpId = Exp[0]->GetExpIndex(gloCoord, NekConstants::kGeomFactorsTol);
        

        for (int j = 0; j < nfields; ++j)
        {
            Exp[j]->PutPhysInToElmtExp();
            cout << "   " << Exp[j]->GetExp(ExpId)->PhysEvaluate(gloCoord);
        }
        cout << endl;
    }
     */
    //----------------------------------------------
    return 0;
}

