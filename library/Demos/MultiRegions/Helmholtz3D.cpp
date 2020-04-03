///////////////////////////////////////////////////////////////////////////////
//
// File: Helmholtz3D.cpp
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
// Description: 3D Helmholtz solver demo.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField3D.h>
#include <SpatialDomains/MeshGraph.h>

using namespace std;
using namespace Nektar;

//#define TIMING
#ifdef TIMING
#include <time.h>
#define Timing(s) \
 fprintf(stdout,"%s Took %g seconds\n",s,(clock()-st)/cps);  \
 st = clock();
#else
#define Timing(s) \
 /* Nothing */
#endif

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    MultiRegions::ContField3DSharedPtr Exp, Fce;
    MultiRegions::ExpListSharedPtr DerExp1, DerExp2, DerExp3;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;
    StdRegions::VarCoeffMap varcoeffs;
#ifdef TIMING
    NekDouble st;
    NekDouble cps = (double)CLOCKS_PER_SEC;
#endif

    if(argc < 2)
    {
        fprintf(stderr,"Usage: Helmholtz3D  meshfile [solntype]\n");
        exit(1);
    }

    try
    {
        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateDefault(vSession);

        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph3D = 
            SpatialDomains::MeshGraph::Read(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        factors[StdRegions::eFactorLambda]             = 
            vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionMap &expansions = 
            graph3D->GetExpansions();
        LibUtilities::BasisKey              bkey0      = 
            expansions.begin()->second->m_basisKeyVector[0];
        
        if (vSession->GetComm()->GetRank() == 0)
        {
            cout << "Solving 3D Helmholtz:"  << endl;
            cout << "  - Communication: " 
                 << vSession->GetComm()->GetType() << " (" 
                 << vSession->GetComm()->GetSize() 
                 << " processes)" << endl;
            cout << "  - Solver type  : " 
                 << vSession->GetSolverInfo("GlobalSysSoln") << endl;
            cout << "  - Lambda       : " 
                 << factors[StdRegions::eFactorLambda] << endl;
            cout << "  - No. modes    : " 
                 << bkey0.GetNumModes() << endl;
            cout << endl;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField3D>
            ::AllocateSharedPtr(vSession, graph3D, vSession->GetVariable(0));
        //----------------------------------------------

        Timing("Read files and define exp ..");

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

        xc0 = Array<OneD, NekDouble>(nq, 0.0);
        xc1 = Array<OneD, NekDouble>(nq, 0.0);
        xc2 = Array<OneD, NekDouble>(nq, 0.0);

        switch(coordim)
        {
            case 3:
                Exp->GetCoords(xc0,xc1,xc2);
                break;
            default:
                ASSERTL0(false,"Coordim not valid");
                break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Set up variable coefficients if defined
        if (vSession->DefinesFunction("d00"))
        {
            Array<OneD, NekDouble> d00(nq,0.0);
            LibUtilities::EquationSharedPtr d00func =
                vSession->GetFunction("d00",0);
            d00func->Evaluate(xc0, xc1, xc2, d00);
            varcoeffs[StdRegions::eVarCoeffD00] = d00;
        }
        if (vSession->DefinesFunction("d11"))
        {
            Array<OneD, NekDouble> d11(nq,0.0);
            LibUtilities::EquationSharedPtr d11func =
                vSession->GetFunction("d11",0);
            d11func->Evaluate(xc0, xc1, xc2, d11);
            varcoeffs[StdRegions::eVarCoeffD11] = d11;
        }
        if (vSession->DefinesFunction("d22"))
        {
            Array<OneD, NekDouble> d22(nq,0.0);
            LibUtilities::EquationSharedPtr d22func =
                vSession->GetFunction("d22",0);
            d22func->Evaluate(xc0, xc1, xc2, d22);
            varcoeffs[StdRegions::eVarCoeffD22] = d22;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        fce = Array<OneD,NekDouble>(nq);
        LibUtilities::EquationSharedPtr ffunc =
            vSession->GetFunction("Forcing", 0);

        ffunc->Evaluate(xc0, xc1, xc2, fce);
        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------

        //---------------------------------------------- 
        //Helmholtz solution taking physical forcing after setting
        //initial condition to zero
        Vmath::Zero(Exp->GetNcoeffs(),Exp->UpdateCoeffs(),1);
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors,
                       varcoeffs);
        //----------------------------------------------
        Timing("Helmholtz Solve ..");

#ifdef TIMING
        for(i = 0; i < 20; ++i)
        {
            Vmath::Zero(Exp->GetNcoeffs(),Exp->UpdateCoeffs(),1);
            Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), factors);
        }

        Timing("20 Helmholtz Solves:... ");
#endif

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
        Exp->BwdTrans(
            Exp->GetCoeffs(), Exp->UpdatePhys());
        //----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol =
            vSession->GetFunction("ExactSolution", 0);

        //-----------------------------------------------
        // Write solution to file
        string out = vSession->GetSessionName() + ".fld";
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
            Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(xc0, xc1, xc2, fce);

            //----------------------------------------------

            //--------------------------------------------
            // Calculate L_inf error
            Fce->SetPhys(fce);
            Fce->SetPhysState(true);

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Exp->GetPhys(), Fce->GetPhys());
            NekDouble vL2Error   = Exp->L2(Exp->GetPhys(), Fce->GetPhys());
            NekDouble vH1Error   = Exp->H1(Exp->GetPhys(), Fce->GetPhys());
            if (vComm->GetRank() == 0)
            {
                cout << "L infinity error: " << vLinfError << endl;
                cout << "L 2 error:        " << vL2Error << endl;
                cout << "H 1 error:        " << vH1Error << endl;
                cout << "Time in ExpEval:  " << ex_sol->GetTime() << endl;
            }
            //--------------------------------------------
        }
        //----------------------------------------------
    }
    catch (const std::runtime_error&)
    {
        cout << "Caught an error" << endl;
        return 1;
    }

    vComm->Finalise();

    return 0;
}
