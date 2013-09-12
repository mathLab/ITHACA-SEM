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

#include <cstdlib>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ContField2D.h>

using namespace Nektar;

class Operators
{
    public:
        Operators(int argc, char* argv[]);
        ~Operators();

        void TimeIntegrate();

        void DoCopy(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);
        void DoImplicitSolve(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                const NekDouble time,
                const NekDouble lambda);

    private:
        LibUtilities::SessionReaderSharedPtr session;
        string sessionName;
        string fileName;
        SpatialDomains::MeshGraphSharedPtr graph;
        MultiRegions::ContField2DSharedPtr field;

        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        LibUtilities::TimeIntegrationSchemeOperators ode;
        Array<OneD, Array<OneD, NekDouble> > fields;

        unsigned int nSteps;
        NekDouble delta_t;
        NekDouble epsilon;
        NekDouble lambda;

        void WriteSolution();
        void ExactSolution();
};


int main(int argc, char *argv[])
{
    Operators ops(argc, argv);
    ops.TimeIntegrate();
}

Operators::Operators(int argc, char* argv[])
{
//    LibUtilities::SessionReaderSharedPtr session;
//    string sessionName;
//    string fileName;
//    SpatialDomains::MeshGraphSharedPtr graph;
//    MultiRegions::ContField2DSharedPtr field;
//
//    Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
//    LibUtilities::TimeIntegrationSolutionSharedPtr u;
//    LibUtilities::TimeIntegrationSchemeOperators ode;
//    Array<OneD, Array<OneD, NekDouble> > fields(1);
    fields = Array<OneD, Array<OneD, NekDouble> >(1);
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
            ::AllocateSharedPtr(session, graph, session->GetVariable(0),
                                true, false);
        fields[0] = field->UpdatePhys();

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
        nSteps = session->GetParameter("NumSteps");
        delta_t = session->GetParameter("TimeStep");
        epsilon = session->GetParameter("epsilon");
        lambda = 1.0/delta_t/epsilon;
        cout << "Lambda: " << lambda << endl;
    }
    catch (const std::runtime_error& e)
    {
        exit(-1);
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
        exit(-1);
    }
}

Operators::~Operators()
{
    // Finalise session
    session->Finalise();
}
void Operators::TimeIntegrate()
{
    int nq = field->GetNpoints();
    Array<OneD, NekDouble> u_old(nq, 0.0);
    Array<OneD, NekDouble> u_tmp(nq, 0.0);
    // Set up time integration
    int numMultiSteps = 2;
    NekDouble curtime = 0.0;

    IntScheme = Array<OneD, LibUtilities::
        TimeIntegrationSchemeSharedPtr>(numMultiSteps);

    // Used in the first time step to initalize the scheme
    LibUtilities::TimeIntegrationSchemeKey IntKey0(
            LibUtilities::eBackwardEuler);

    // Used for all other time steps
    LibUtilities::TimeIntegrationSchemeKey IntKey1(
            LibUtilities::eBDFImplicitOrder2);
    IntScheme[0] = LibUtilities::
        TimeIntegrationSchemeManager()[IntKey0];
    IntScheme[1] = LibUtilities::
        TimeIntegrationSchemeManager()[IntKey1];

    ode.DefineImplicitSolve(&Operators::DoImplicitSolve, this);
cout << Vmath::Dot(nq, fields[0], 1, fields[0], 1) << endl;
    // Initialise the scheme for actual time integration scheme
    u = IntScheme[1]->InitializeScheme(delta_t, fields, curtime, ode);
cout << Vmath::Dot(nq, u->GetSolutionVector()[0][0], 1, u->GetSolutionVector()[0][0], 1) << endl;

    for (int n = 0; n < nSteps; ++n)
    {
        cout << "Step " << n << endl;
        cout << "u[0]" << Vmath::Dot(nq, u->GetSolutionVector()[0][0], 1, u->GetSolutionVector()[0][0], 1) << endl;
        cout << "u[1]" << Vmath::Dot(nq, u->GetSolutionVector()[1][0], 1, u->GetSolutionVector()[1][0], 1) << endl;
        fields = IntScheme[std::min(n, numMultiSteps-1)]->TimeIntegrate(
                            delta_t, u, ode);
        cout << "u[0]" << Vmath::Dot(nq, u->GetSolutionVector()[0][0], 1, u->GetSolutionVector()[0][0], 1) << endl;
        cout << "u[1]" << Vmath::Dot(nq, u->GetSolutionVector()[1][0], 1, u->GetSolutionVector()[1][0], 1) << endl;
        curtime += delta_t;
    }
    Vmath::Vcopy(nq, fields[0], 1, field->UpdatePhys(), 1);

    WriteSolution();
    ExactSolution();

}

void Operators::DoCopy(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
{
    unsigned int n = inarray[0].num_elements();
    for (int i = 0; i < inarray.num_elements(); ++i)
    {
        Vmath::Vcopy(n, inarray[i], 1, outarray[i], 1);
    }
}


void Operators::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble lambda)
{
    cout << "Do Implicit Solve: " << endl;
    cout << "   time = " << time << endl;
    cout << "   lambda = " << lambda << endl;
    int nvariables = 1;
    int nq = field->GetNpoints();
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0/lambda/epsilon;
    factors[StdRegions::eFactorTau] = 1.0;
Array<OneD, NekDouble> tmp(nq, 0.0);
Array<OneD, NekDouble> tmp2(field->GetNcoeffs(), 0.0);
    cout << "eFactorLambda = " << factors[StdRegions::eFactorLambda] << endl;
    for (int i = 0; i < nvariables; ++i)
    {
        // Multiply 1.0/timestep/lambda
        Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1, tmp, 1);
cout << "before helmsolve: " << Vmath::Dot(nq, tmp, tmp) << endl;
        // Solve a system of equations with Helmholtz solver
        field->HelmSolve(tmp, tmp2,NullFlagList,factors);
        field->BwdTrans(tmp2, tmp);
cout << "after helmsolve: " << Vmath::Dot(nq, tmp, tmp) << endl;
        field->SetPhysState(false);

        // The solution is Y[i]
        outarray[i] = tmp;
    }
}

void Operators::WriteSolution()
{
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

}

void Operators::ExactSolution()
{
    //----------------------------------------------
    // See if there is an exact solution, if so
    // evaluate and plot errors
    unsigned int nq = field->GetNpoints();
    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);
    field->GetCoords(x0,x1,x2);

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

}
