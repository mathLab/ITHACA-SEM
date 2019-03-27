///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationDemo.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Example of using time-integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

//--------------------------------------------------
// This file illustrates the use of the time-stepping framework.
//
// For more information, please consult the following reference:
//
// Vos, P.E.J., Eskilsson, C., Bolis, A., Chun, S., Kirby, R.M. and Sherwin, S.J.
// "A Generic Framework for Time-Stepping PDEs: general linear methods,
//  object-oriented implementation and application to fluid problems"
// International Journal of Computational Fluid Dynamics, to appear

// It solves the one-dimensional advection-diffusion problem, defined as
//
//  |    du     du     d^2 u
//  |    -- + V -- = D -----,
//  |    dt     dx     d x^2
//  |
//  |    subject to:
//  |    - periodic boundary conditions
//  |    - the initial condition
//  |        u(x,0) = sin(2*pi*k*x)
//  |
//  |    and with  x = [0,1]
//  |              t = [0,1]
//  |              U = 1
//  |              D = 0.05
//  |              k = 1
//  |
//
// using the finite difference method.
// The exact solution of the equation above is given by
//
//  u(x,t) = exp(-D * (2*pi*k)^2 * t) * sin(2*pi*k * (x - U*t) )
//
// The output is written out to the files
//
//   - OneDfinDiffAdvDiffSolverOutput.dat (containing the data)
//   - OneDfinDiffAdvDiffSolverOutput.p   (containing a gnuplot script)
//
// and can be visualised by gnuplot using the command
//
//    gnuplot OneDfinDiffAdvDiffSolverOutput.p
//
//--------------------------------------------------
#include <fstream>
#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegratorBase.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

namespace po = boost::program_options;

// We first implement a class that represents
// the 1D finite difference solver
class OneDfinDiffAdvDiffSolver
{
public:
    // constructor based upon the discretisation details
    OneDfinDiffAdvDiffSolver(int nPoints, int nTimeSteps):
        m_x0(0.0),
        m_xend(1.0),
        m_nPoints(nPoints),
        m_dx((m_xend-m_x0)/((double) m_nPoints-1.0)),
        m_t0(0.0),
        m_tend(1.0),
        m_nTimeSteps(nTimeSteps),
        m_dt((m_tend-m_t0)/(double) m_nTimeSteps),
        m_wavenumber(1.0),
        m_U(1.0),
        m_D(0.05)
    {
    }

    // -----------------------------------------------------------------
    // ---- These functions/methods below are the routines which will be used by
    // ---- the TimeIntegration framework (and are required for using it)...
    // ---- The implementation of these functions can be found at the end of the file
    void HelmSolve(const Array<OneD, const Array<OneD, double> >& inarray,
                         Array<OneD,       Array<OneD, double> >& outarray,
                   const NekDouble time,
                   const NekDouble lambda) const;

    void EvaluateAdvectionTerm(const Array<OneD, const  Array<OneD, double> >& inarray,
                                     Array<OneD,        Array<OneD, double> >& outarray,
                               const NekDouble time) const;

    void Project(const Array<OneD, const  Array<OneD, double> >& inarray,
                       Array<OneD,        Array<OneD, double> >& outarray,
                 const NekDouble time) const;
    // -----------------------------------------------------------------


    // -----------------------------------------------------------------
    // Below are some additional functions
    int GetNpoints() const;

    void EvaluateExactSolution(Array<OneD, Array<OneD, double> >& outarray,
                               const NekDouble time) const;

    double EvaluateL2Error(const Array<OneD, const  Array<OneD, double> >& approx,
                           const Array<OneD, const  Array<OneD, double> >& exact) const;

    void AppendOutput(ofstream& outfile,
                      const Array<OneD, const  Array<OneD, double> >& approx,
                      const Array<OneD, const  Array<OneD, double> >& exact) const;

    void GenerateGnuplotScript() const;

    double GetInitialTime() const;

    double GetTimeStep() const;
    // -----------------------------------------------------------------

private:
    // spatial discretisation
    double m_x0;         // the left boundary of the domain
    double m_xend;       // the right boundary of the domain
    int    m_nPoints;    // the number of grid-points used in the finite difference method
    double m_dx;         // the distance between 2 grid points

    // temporal discretisation
    double m_t0;         // the initial time
    double m_tend;       // the end time
    int    m_nTimeSteps; // the number of time-steps
    double m_dt;         // the size of a time-step

    // value of the coefficients
    double m_wavenumber; // wave number
    double m_U;          // advection speed
    double m_D;          // diffusion coefficient


    void solveTriDiagMatrix (int n, double alpha, double beta,
                             const Array<OneD, const double>& inarray,
                             Array<OneD,       double>& outarray) const;
};


int main(int argc, char *argv[])
{
    po::options_description desc("Usage:");
    desc.add_options()
        ("help, h",      "Display this help message.")
        ("Npoints, np",    po::value<int>(),
                          "Specify the number of grid points to be used.")
        ("Ntimesteps, nt", po::value<int>(),
                          "Specify the number of timesteps to be used.")
        ("TimeIntegrationMethod, m", po::value<int>(),
                    "TimeIntegrationMethod is a number in the range [1,8]\n"
                    "and defines the time-integration method to be used:\n"
                    "- 1: 1st order multi-step IMEX scheme (Euler Backwards/Euler Forwards)\n"
                    "- 2: 2nd order multi-step IMEX scheme\n"
                    "- 3: 3rd order multi-step IMEX scheme\n"
                    "- 4: 2nd order multi-stage DIRK IMEX scheme\n"
                    "- 5: 3nd order multi-stage DIRK IMEX scheme\n"
                    "- 6: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
                    "- 7: 2nd order Crank-Nicolson/Adams-Bashforth (CNAB)\n"
                    "- 8: 2nd order Modified Crank-Nicolson/Adams-Bashforth (MCNAB)" );
    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(desc).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    if( !vm.count("Npoints") || !vm.count("Ntimesteps") || !vm.count("TimeIntegrationMethod") || vm.count("help") )
    {
        cout << "Please specify points, timesteps and method." << endl << endl;
        cout << desc;
        return 1;
    }

    int nPoints = vm["Npoints"].as<int>();
    int nTimesteps = vm["Ntimesteps"].as<int>();
    int nMethod = vm["TimeIntegrationMethod"].as<int>();

    // Open a file for writing the solution
    ofstream outfile;
    outfile.open("OneDfinDiffAdvDiffSolverOutput.dat");



    // -----------------------------------------------------------------------------
    // THE IMPLEMENTATION BELOW SHOWS HOW THE TIME-STEPPING FRAMEWORK CAN BE
    // USED FOR TIME-INTEGRATION PDEs

    // 1. THE SPATIAL DISCRETISATION
    //    Create an object of the OneDfinDiffAdvDiffSolver class.
    //    This class can be thought of as representing the
    //    spatial (finite difference) discretisation.
    OneDfinDiffAdvDiffSolver* solver = new OneDfinDiffAdvDiffSolver( nPoints, nTimesteps );
    //    After this spatial discretisation, the PDE has actually been
    //    reduced (through the method-of-lines) to an ODE.
    //    In order to use the time-stepping framework, we need to give it the necessary
    //    information about this ODE.
    //    Therefore, we create an oject of the class TimeIntegrationSchemeOperators that
    //    contains a 'function pointer' (in fact a 'functor') to the
    //    - explicit term of the ODE (i.e. the advection term)
    //    - implicit solve routine (i.e. the Helmholtz solver)
    //    - projection operator (i.e. the identity operator in this case)
    LibUtilities::TimeIntegrationSchemeOperators ode;
    ode.DefineOdeRhs       ( &OneDfinDiffAdvDiffSolver::EvaluateAdvectionTerm, solver );
    ode.DefineImplicitSolve( &OneDfinDiffAdvDiffSolver::HelmSolve,             solver );
    ode.DefineProjection   ( &OneDfinDiffAdvDiffSolver::Project,               solver );

    // 2. THE TEMPORAL DISCRETISATION
    // 2.1 Read in which method should be used.  Create time integrator
    //
    LibUtilities::TimeIntegratorSharedPtr timeIntegrator;

    switch( nMethod )
    {
    case 1 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXOrder1" );     break;
    case 2 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXOrder2" );     break;
    case 3 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXOrder3" );     break;
    case 4 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXdirk_2_3_2" ); break;
    case 5 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXdirk_3_4_3" ); break;
    case 6 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "IMEXGear" );       break;
    case 7 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "CNAB" );           break;
    case 8 : timeIntegrator = LibUtilities::GetTimeIntegratorFactory().CreateInstance( "MCNAB" );          break;
    default :
      {
        cerr << "The third argument defines the time-integration method to be used:\n\n";
        cout << desc;
        exit( 1 );
      }
    }

    // 2.2 Initialise the arrays that contain the numerical and analytical solutions:
    //
    double t0 = solver->GetInitialTime();

    // Note, fidifsol (finite difference solution) is an array of arrays, though in this current
    // example, we only ever use one array (ie: fidifsol[0]).  This is because, in this demo, we are only
    // solving for one variable.  

    Array<OneD, Array<OneD, double> > fidifsol(1); // Array containing the numerical solution
    Array<OneD, Array<OneD, double> > exactsol(1); // Array containing the exact solution
    fidifsol[0] = Array<OneD, double>( nPoints );
    exactsol[0] = Array<OneD, double>( nPoints );
    solver->EvaluateExactSolution( fidifsol, t0 );  // Set the initial condition

    // I've commented out the following two lines as they duplicate what is done in the
    // first iteration of the below for loop (part 2.5).  
    // 
    //    solver->EvaluateExactSolution(exactsol,t0);  // Set the initial condition
    //    solver->AppendOutput(outfile,fidifsol,exactsol); // Write the initial condition to a file

    // 2.3 Initialize the time-integration scheme:
    //
    double dt = solver->GetTimeStep();

    LibUtilities::TimeIntegrationSolutionSharedPtr sol;
    sol = timeIntegrator->InitializeIntegrator( dt, fidifsol, t0, ode );

    // 2.4 Do the time-integration:
    //
    for( int i = 0; i < nTimesteps; i++ )
    {
        double t = t0 + i*dt;

        fidifsol = timeIntegrator->TimeIntegrate( i, dt, sol, ode ); // Time-integration for 1 time-step

        solver->EvaluateExactSolution( exactsol, t );        // Calculate the exact solution
        solver->AppendOutput( outfile, fidifsol, exactsol ); // Dump the output to a file
    }

    // Calculate the error and dump to screen
    cout << "L 2 error :" <<  solver->EvaluateL2Error(fidifsol,exactsol) << endl;

    // Some more writing out the results
    solver->GenerateGnuplotScript();
    outfile.close();

    delete solver;

    return 0;
}

void OneDfinDiffAdvDiffSolver::HelmSolve(const Array<OneD, const Array<OneD, double> >& inarray,
                                               Array<OneD,       Array<OneD, double> >& outarray,
                                         const NekDouble time,
                                         const NekDouble lambda) const
{
    // This function implements a 1D finite difference helmholtz solver.
    // The 1D Helmholtz equation leads to a cyclic triadiagonal matrix to be
    // solved.
    // In this function, this is solved as follows:
    // - Based upon the static condensation technique, we first solve
    //   for the periodic grid-points (i.e. the boundary nodes)
    // - Next, we solve for the interior grid-points. The associated tridiagonal
    //   system is solved based upon the Thomas algorithm
    double a = - m_D * lambda / (m_dx * m_dx); //off diagonal term
    double b  = 1.0 + 2.0 * lambda * m_D / (m_dx * m_dx); // diagonal term of triadiagonal matrix

    int nIntPoints = m_nPoints-2;

    Array<OneD, double> invD_f(nIntPoints);
    solveTriDiagMatrix(nIntPoints,a,b,inarray[0]+1,invD_f);

    Array<OneD, double> C(nIntPoints,0.0);
    Array<OneD, double> invD_C(nIntPoints,0.0);
    C[0]     = a;
    C[nIntPoints-1] = a;
    solveTriDiagMatrix(nIntPoints,a,b,C,invD_C);

    outarray[0][0] = ( inarray[0][0] - a*(invD_f[0]+invD_f[nIntPoints-1]) ) / (b - a*(invD_C[0]+invD_C[nIntPoints-1]) );
    outarray[0][m_nPoints-1] = outarray[0][0];

    Array<OneD, double> f(nIntPoints);
    for(int i = 0; i < nIntPoints; i++)
    {
        f[i] = inarray[0][i+1];
    }
    f[0] -= outarray[0][0]*a;
    f[nIntPoints-1] -= outarray[0][0]*a;

    Array<OneD, double> tmp;
    solveTriDiagMatrix(nIntPoints,a,b,f,tmp = outarray[0]+1); // Calls the Thomas algorithm
}

void OneDfinDiffAdvDiffSolver::EvaluateAdvectionTerm(const Array<OneD, const  Array<OneD, double> >& inarray,
                                                           Array<OneD,        Array<OneD, double> >& outarray,
                                                     const NekDouble time) const
{
    // The advection term can be evaluated using central or upwind differences
    if( true )
    {
        // Note: We are using a periodic boundary condition where the 1st point and last point are actually
        // the same point.  This is why the 1st (and last) point in the output array (index 0 and m_nPoints-1 respectively)
        // are NOT used for the central differences, and instead the 2nd point (index 1) and 2nd to last point
        // are used.
      
        // Central differences:
        outarray[0][0]           = - m_U * (inarray[0][1]-inarray[0][m_nPoints-2]) / (2.0 * m_dx);
        outarray[0][m_nPoints-1] = outarray[0][0];

        for(int i = 1; i < m_nPoints-1; i++)
        {
            outarray[0][i] = - m_U * (inarray[0][i+1]-inarray[0][i-1]) / (2.0 * m_dx);
        }
    }
    else
    {
        // upwind differences
        for(int i = 1; i < m_nPoints; i++)
        {
            outarray[0][i] = - m_U * (inarray[0][i]-inarray[0][i-1]) / (m_dx);
        }
        outarray[0][0] = outarray[0][m_nPoints-1];

    }

}

void OneDfinDiffAdvDiffSolver::Project(const Array<OneD, const  Array<OneD, double> >& inarray,
                                             Array<OneD,        Array<OneD, double> >& outarray,
                                       const NekDouble time) const
{
    // This is simply the identity operator for this case
    for(int i = 0; i < m_nPoints; i++)
    {
        outarray[0][i] = inarray[0][i];
    }
}

void OneDfinDiffAdvDiffSolver::solveTriDiagMatrix (int n, double a, double b,
                                                   const Array<OneD, const double>& inarray,
                                                         Array<OneD,       double>& outarray) const
{
    // Implementation of the Thomas algorithm for Tridiaginol systems
    Array<OneD, double> cprime(n);
    Array<OneD, double> dprime(n);

    cprime[0] = a / b;
    dprime[0] = inarray[0] / b;

    for(int i = 1; i < n; i++)
    {
        double id = 1.0 / (b - cprime[i-1] * a);
        cprime[i] = a * id;
        dprime[i] = (inarray[i] - dprime[i-1] * a) * id;
    }

    outarray[n-1] = dprime[n-1];
    for (int i = n - 2; i >= 0; i--)
    {
        outarray[i] = dprime[i] - cprime[i] * outarray[i+1];
    }

}
int OneDfinDiffAdvDiffSolver::GetNpoints() const
{
    return m_nPoints;
}

void OneDfinDiffAdvDiffSolver::EvaluateExactSolution(Array<OneD, Array<OneD, double> >& outarray,
                                                     const NekDouble time) const
{
    double x;
    for(int i = 0; i < m_nPoints; i++)
    {
        x = m_x0 + i*m_dx;
        outarray[0][i] = exp(-m_D * 2.0 * 2.0 * M_PI * M_PI * m_wavenumber*m_wavenumber*time) *
            sin( 2.0 * m_wavenumber * M_PI * (x - m_U*time) );
    }
}

// Note, this routine returns the Relative Error L2 Norm (as opposed to the absolute L2 norm)...
double OneDfinDiffAdvDiffSolver::EvaluateL2Error(const Array<OneD, const  Array<OneD, double> >& approx,
                                                 const Array<OneD, const  Array<OneD, double> >& exact) const
{
    double a = 0.0;
    double b = 0.0;

    for( int i = 0; i < m_nPoints; i++ )
    {
        a += ( approx[0][i] - exact[0][i] ) * ( approx[0][i] - exact[0][i] );
        b += exact[0][i] * exact[0][i];
    }

    // Note: Returning Relative Error L2 Norm.
    return sqrt( a / b );
}

void OneDfinDiffAdvDiffSolver::AppendOutput(ofstream& outfile,
                                            const Array<OneD, const  Array<OneD, double> >& approx,
                                            const Array<OneD, const  Array<OneD, double> >& exact) const
{
    for(int i = 0; i < m_nPoints; i++)
    {
        outfile << scientific
                << setw (17)
                << setprecision(10)
                << m_x0 + i*m_dx
                << "  "
                << approx[0][i]
                << "  "
                << exact[0][i]
                << endl;
    }
    outfile << endl << endl;
}

void OneDfinDiffAdvDiffSolver::GenerateGnuplotScript() const
{
    ofstream outfile;
    outfile.open("OneDfinDiffAdvDiffSolverOutput.p");

    outfile << "# Gnuplot script file" << endl;
    outfile << "set   autoscale" << endl;
    outfile << "unset log" << endl;
    outfile << "unset label" << endl;
    outfile << "set xtic auto" << endl;
    outfile << "set ytic auto" << endl;
    outfile << "set title \"Finite Difference Solution to the 1D advection-diffusion equation\"" << endl;
    outfile << "set xlabel \"x\"" << endl;
    outfile << "set ylabel \"u\"" << endl;
    outfile << "set xr [" << m_x0 << ":" << m_xend << "]" << endl;
    outfile << "set yr [-1.0:1.0]" << endl;

    double t;
    for(int i=0; i <= m_nTimeSteps; i++)
    {
        t = m_t0+i*m_dt;
        outfile << "plot    \"OneDfinDiffAdvDiffSolverOutput.dat\" ";
        outfile << "using 1:2 index ";
        outfile << i << " title 'Finite Difference Solution (t=" << t << ")' with linespoints , ";
        outfile << "\"OneDfinDiffAdvDiffSolverOutput.dat\" ";
        outfile << "using 1:3 index ";
        outfile << i << " title 'Exact Solution (t=" << t << ")' with linespoints" << endl;
        outfile << "pause " << 4.0/m_nTimeSteps << endl;
    }

    outfile.close();
}

double OneDfinDiffAdvDiffSolver::GetInitialTime() const
{
    return m_t0;
}

double OneDfinDiffAdvDiffSolver::GetTimeStep() const
{
    return m_dt;
}

