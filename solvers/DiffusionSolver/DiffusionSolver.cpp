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
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ContField2D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    LibUtilities::FieldIOSharedPtr       fld;
    SpatialDomains::MeshGraphSharedPtr   graph;
    MultiRegions::ContField2DSharedPtr   field;
    LibUtilities::EquationSharedPtr      icond, ex_sol;
    StdRegions::ConstFactorMap           factors;

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Read the geometry and the expansion information
        graph = SpatialDomains::MeshGraph::Read(session);

        // Create Field I/O object.
        fld     = LibUtilities::FieldIO::CreateDefault(session);

        // Get some information about the session
        string       sessionName = session->GetSessionName();
        string       outFile     = sessionName + ".fld";
        unsigned int nSteps      = session->GetParameter("NumSteps");
        NekDouble    delta_t     = session->GetParameter("TimeStep");
        NekDouble    epsilon     = session->GetParameter("epsilon" );

        // Create field
        field = MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(session, graph, session->GetVariable(0));

        // Get coordinates of physical points
        unsigned int nq = field->GetNpoints();
        Array<OneD,NekDouble> x0(nq), x1(nq), x2(nq);
        field->GetCoords(x0,x1,x2);

        // Evaluate initial condition at these points
        icond = session->GetFunction("InitialConditions", "u");
        icond->Evaluate(x0, x1, x2, 0.0, field->UpdatePhys());

        // Compute lambda in the Helmholtz problem
        factors[StdRegions::eFactorLambda] = 1.0/delta_t/epsilon;

        // Zero field coefficients for initial guess for linear solver.
        Vmath::Zero(field->GetNcoeffs(), field->UpdateCoeffs(), 1);

        // Time integrate using backward Euler
        for (unsigned int n = 0; n < nSteps; ++n)
        {
            Vmath::Smul(nq, -1.0/delta_t/epsilon, field->GetPhys(),    1,
                                                  field->UpdatePhys(), 1);

            field->HelmSolve(field->GetPhys(), field->UpdateCoeffs(), factors);

            field->BwdTrans(field->GetCoeffs(), field->UpdatePhys());
        }

        // Write solution to file
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                        = field->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            field->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        fld->Write(outFile, FieldDef, FieldData);

        // Check for exact solution
        ex_sol = session->GetFunction("ExactSolution",0);
        if(ex_sol)
        {
            // Allocate storage
            Array<OneD, NekDouble> exact(nq);

            //----------------------------------------------
            // evaluate exact solution
            ex_sol->Evaluate(x0, x1, x2, (nSteps)*delta_t, exact);

            //--------------------------------------------
            // Calculate errors
            cout << "L inf error:      "
                 << field->Linf(field->GetPhys(), exact) << endl;
            cout << "L 2 error:        "
                 << field->L2(field->GetPhys(), exact) << endl;
            cout << "H 1 error:        "
                 << field->H1(field->GetPhys(), exact) << endl;
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
