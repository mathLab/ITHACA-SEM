///////////////////////////////////////////////////////////////////////////////
//
// File DiffusionTestTI.cpp
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
// Description: Diffusion solver
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>

#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ContField2D.h>

using namespace std;
using namespace Nektar;

class Diffusion
{
    public:
        Diffusion( int argc, char* argv[] );
        ~Diffusion();

        void TimeIntegrate();

        void DoImplicitSolve( const Array<OneD, const Array<OneD, NekDouble> > & inarray,
                                    Array<OneD, Array<OneD, NekDouble> >       & outarray,
                              const NekDouble                                    time,
                              const NekDouble                                    lambda );

    private:
        LibUtilities::SessionReaderSharedPtr            session;
        LibUtilities::FieldIOSharedPtr                  fld;
        string                                          sessionName;
        SpatialDomains::MeshGraphSharedPtr              graph;
        MultiRegions::ContField2DSharedPtr              field;

        LibUtilities::TimeIntegrationSchemeSharedPtr    m_IntScheme;
        LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr  m_u;
        LibUtilities::TimeIntegrationSchemeOperators    ode;
        Array<OneD, Array<OneD, NekDouble> >            fields;

        string                                          m_scheme_name;
        unsigned int                                    nSteps;
        NekDouble                                       delta_t;
        NekDouble                                       epsilon;
        NekDouble                                       lambda;

        void WriteSolution();
        void ExactSolution();
};


Diffusion::Diffusion( int argc, char* argv[] )
{
    // Create session reader.
    session     = LibUtilities::SessionReader::CreateInstance( argc, argv );

    // Read the geometry and the expansion information
    graph       = SpatialDomains::MeshGraph::Read(session);

    // Create Field I/O object.
    fld         = LibUtilities::FieldIO::CreateDefault(session);

    // Get some information from the session
    sessionName   = session->GetSessionName();
    m_scheme_name = session->GetSolverInfo( "TimeIntegrationMethod" );
    nSteps        = session->GetParameter( "NumSteps" );
    delta_t       = session->GetParameter( "TimeStep" );
    epsilon       = session->GetParameter( "epsilon" );
    lambda        = 1.0 / delta_t / epsilon;

    // Set up the field
    field       = MemoryManager<MultiRegions::ContField2D>::
                    AllocateSharedPtr(session, graph, session->GetVariable(0));

    fields      = Array<OneD, Array<OneD, NekDouble> >(1);
    fields[0]   = field->UpdatePhys();

    // Get coordinates of physical points
    unsigned int nq = field->GetNpoints();
    Array<OneD,NekDouble> x0(nq), x1(nq), x2(nq);
    field->GetCoords(x0,x1,x2);

    // Evaluate initial condition
    LibUtilities::EquationSharedPtr icond
        = session->GetFunction("InitialConditions", "u");
    icond->Evaluate(x0,x1,x2,0.0,field->UpdatePhys());
}

Diffusion::~Diffusion()
{
    session->Finalise();
}

void Diffusion::TimeIntegrate()
{
    LibUtilities::TimeIntegrationSchemeFactory & fac =
        LibUtilities::GetTimeIntegrationSchemeFactory();

    m_IntScheme = fac.CreateInstance( m_scheme_name, "", 0,
				      std::vector<NekDouble>() );

    ode.DefineImplicitSolve( &Diffusion::DoImplicitSolve, this );

    // Initialise the scheme for actual time integration scheme
    m_u = m_IntScheme->InitializeScheme( delta_t, fields, 0.0, ode );

    // Zero field coefficients for initial guess for linear solver.
    Vmath::Zero(field->GetNcoeffs(), field->UpdateCoeffs(), 1);

    for (int n = 0; n < nSteps; ++n)
    {
      fields = m_IntScheme->TimeIntegrate( n, delta_t, m_u, ode );
    }
    Vmath::Vcopy(field->GetNpoints(), fields[0], 1, field->UpdatePhys(), 1);

    WriteSolution();
    ExactSolution();
}


void Diffusion::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble lambda)
{
    boost::ignore_unused(time);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0/lambda/epsilon;

    for (int i = 0; i < inarray.size(); ++i)
    {
        // Multiply RHS by 1.0/timestep/lambda
        Vmath::Smul(field->GetNpoints(), -factors[StdRegions::eFactorLambda],
                                         inarray [i], 1,
                                         outarray[i], 1);

        // Solve a system of equations with Helmholtz solver
        field->HelmSolve(outarray[i],
                         field->UpdateCoeffs(), factors);

        // Transform to physical space and store in solution vector
        field->BwdTrans (field->GetCoeffs(), outarray[i]);
    }
}

void Diffusion::WriteSolution()
{
    // Write solution to file
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                = field->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    for(int i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        field->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    fld->Write(session->GetSessionName() + ".fld", FieldDef, FieldData);

}


void Diffusion::ExactSolution()
{
    unsigned int nq = field->GetNpoints();
    Array<OneD,NekDouble> x0(nq), x1(nq), x2(nq);
    field->GetCoords(x0,x1,x2);

    LibUtilities::EquationSharedPtr ex_sol =
                                session->GetFunction("ExactSolution",0);

    if(ex_sol)
    {
        // evaluate exact solution
        Array<OneD, NekDouble> exact(nq);
        ex_sol->Evaluate(x0, x1, x2, (nSteps)*delta_t, exact);

        // Calculate errors
        cout << "L inf error:      "
             << field->Linf(field->GetPhys(), exact) << endl;
        cout << "L 2 error:        "
             << field->L2(field->GetPhys(), exact) << endl;
        cout << "H 1 error:        "
             << field->H1(field->GetPhys(), exact) << endl;
    }

}

int main(int argc, char *argv[])
{
    try
    {
        Diffusion ops(argc, argv);
        ops.TimeIntegrate();
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
