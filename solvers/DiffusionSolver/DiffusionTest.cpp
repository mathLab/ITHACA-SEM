///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolver.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Advection Diffusion Reaction framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    string sessionName;
    string fileName;
    SpatialDomains::MeshGraphSharedPtr graph;
    MultiRegions::ContField2DSharedPtr field;

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Filename of the session file
        fileName = session->GetFilename();

        // Save the basename of input file name for output details
        sessionName = session->GetSessionName();

        // Read the geometry and the expansion information
        graph = SpatialDomains::MeshGraph::Read(session);


        field = MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(session, graph, session->GetVariable(0));

        // Get coordinates of physical points
        unsigned int nq = field->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
        field->GetCoords(x0,x1,x2);

        // Evaluate initial condition
        LibUtilities::EquationSharedPtr ffunc
            = session->GetFunction("ExactSolution", "u");
        ffunc->Evaluate(x0,x1,x2,0.0,field->UpdatePhys());
        field->FwdTrans(field->GetPhys(), field->UpdateCoeffs());
        field->BwdTrans(field->GetCoeffs(), field->UpdatePhys());

        // Get parameters and set up arrays
        unsigned int nSteps = session->GetParameter("NumSteps");
        NekDouble delta_t = session->GetParameter("TimeStep");
        NekDouble epsilon = session->GetParameter("epsilon");
        Array<OneD, NekDouble> u_old(nq, 0.0);
        Array<OneD, NekDouble> u_tmp(nq, 0.0);
        NekDouble lambda = 1.0/delta_t/epsilon;
        cout << "Lambda: " << lambda << endl;
        StdRegions::ConstFactorMap factors;

        // Perform first step using backward euler
        if (nSteps > 0)
        {
            Vmath::Vcopy(nq, field->GetPhys(), 1, u_old, 1);

            factors[StdRegions::eFactorLambda] = lambda;
            Vmath::Smul(nq, -1.0/delta_t/epsilon, field->GetPhys(), 1, field->UpdatePhys(), 1);
cout << "before helmsolve: " << Vmath::Dot(nq, field->GetPhys(), field->GetPhys()) << endl;
            field->HelmSolve(field->GetPhys(), field->UpdateCoeffs(), NullFlagList, factors);
            field->BwdTrans(field->GetCoeffs(), field->UpdatePhys());
cout << "after helmsolve: " << Vmath::Dot(nq, field->GetPhys(), field->GetPhys()) << endl;

        }

        // Perform subsequent steps using BFDImplicitOrder2
        factors[StdRegions::eFactorLambda] = 3.0*lambda/2.0;
        for (unsigned int n = 1; n < nSteps; ++n)
        {
            cout << "eFactorLambda = " << factors[StdRegions::eFactorLambda] << endl;
            cout << "BDF2 Step " << n << endl;
cout << "u_old norm: " << Vmath::Dot(nq, u_old, 1, u_old, 1) << endl;
cout << "u_n norm  : " << Vmath::Dot(nq, field->GetPhys(), 1, field->GetPhys(), 1) << endl;
            Vmath::Svtvp(nq, -4.0, field->GetPhys(), 1, u_old, 1, u_tmp, 1);
            Vmath::Smul(nq, lambda/2.0, u_tmp, 1, u_tmp, 1);
            Vmath::Vcopy(nq, field->GetPhys(), 1, u_old, 1);
cout << "before helmsolve: " << Vmath::Dot(nq, u_tmp, u_tmp) << endl;
            field->HelmSolve(u_tmp, field->UpdateCoeffs(), NullFlagList, factors);
            field->BwdTrans(field->GetCoeffs(), field->UpdatePhys());
cout << "after helmsolve: " << Vmath::Dot(nq, field->GetPhys(), field->GetPhys()) << endl;
        }


        // Write solution to file
        string   out(session->GetSessionName() + ".fld");
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                    = field->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            field->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph->Write(out, FieldDef, FieldData);

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol = session->GetFunction("ExactSolution",0);
        Array<OneD, NekDouble> exact(nq);
        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(x0, x1, x2, (nSteps)*delta_t, exact);

            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = field->Linf(exact);
            NekDouble vL2Error   = field->L2(exact);
            NekDouble vH1Error   = field->H1(exact);
            cout << "L infinity error: " << vLinfError << endl;
            cout << "L 2 error:        " << vL2Error << endl;
            cout << "H 1 error:        " << vH1Error << endl;
            //--------------------------------------------
        }

        // Finalise session
        session->Finalise();
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
