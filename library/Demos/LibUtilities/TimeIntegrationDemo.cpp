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
//    Peter E.J. Vos, Sehun Chun, Alessandro Bolis, Claes Eskilsson, Robert M.
//    Kirby and Spencer J. Sherwin,
//    "A Generic Framework for Time-Stepping PDEs: General Linear Methods,
//    Object-Oriented Implementations and Applications to Fluid Problems",
//    International Journal of Computational Fluid Dynamics, Vol. 25, Issue 3,
//    pages 107-125, 2011.
//
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
//   - OneDFiniteDiffAdvDiffSolver.dat (containing the data)
//   - OneDFiniteDiffAdvDiffSolver.p   (containing a gnuplot script)
//
// and can be visualised by gnuplot using the command
//
//    gnuplot OneDFiniteDiffAdvDiffSolverOutput.p
//
// -----------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

#include <time.h>       /* time */

using namespace Nektar;
using namespace Nektar::LibUtilities;

namespace po = boost::program_options;

// Base class for the solver.
class DemoSolver
{
public:
    // -----------------------------------------------------------------
    // Constructor based upon the discretisation details
    DemoSolver(int nVars, int nPoints, int nTimeSteps) :
      m_nVars(nVars),

      m_x0(0.0), m_xend(1.0), m_nPoints(nPoints),
      m_dx((m_xend - m_x0) / ((double)m_nPoints - 1.0)),

      m_t0(0.0), m_tend(1.0), m_nTimeSteps(nTimeSteps),
      m_dt((m_tend - m_t0) / (double)m_nTimeSteps)
    {
        m_minValue = +std::numeric_limits<double>::max();
        m_maxValue = -std::numeric_limits<double>::max();
    }

    virtual ~DemoSolver() {};

    // -----------------------------------------------------------------
    // Exact solution and project (identity)
    virtual void EvaluateExactSolution(
        Array<OneD, Array<OneD, double>> &outarray,
        const NekDouble time) const = 0;

    virtual void Project(const Array<OneD, const Array<OneD, double>> &inarray,
                               Array<OneD,       Array<OneD, double>> &outarray,
                         const NekDouble time) const;

    // -----------------------------------------------------------------
    // Misc functions for error and outputing
    void GetMinMaxValues(
        const Array<OneD, const Array<OneD, double>> &approx,
        const Array<OneD, const Array<OneD, double>> &exact,
        bool print );

    void EvaluateL2Error(
        const Array<OneD, const Array<OneD, double>> &approx,
        const Array<OneD, const Array<OneD, double>> &exact );

    void AppendOutput(
        std::ofstream &outfile, const int timeStepNumber, const NekDouble time,
        const Array<OneD, const Array<OneD, double>> &approx,
        const Array<OneD, const Array<OneD, double>> &exact) const;

    void GenerateGnuplotScript(const std::string &method) const;

    // -----------------------------------------------------------------
    // Access methods
    double GetInitialTime() const
    {
        return m_t0;
    }

    double GetFinialTime() const
    {
        return m_tend;
    }

    double GetDeltaT() const
    {
        return m_dt;
    }

    std::string GetFileName() const
    {
        return m_fileName;
    }

    void SetSchemeName( std::string name )
    {
        m_schemeName = name;
    }

    // -----------------------------------------------------------------

protected:

    std::string m_fileName;
    std::string m_schemeName;
    std::string m_title;

    // Number variables
    int m_nVars;

    // Spatial discretisation:
    double m_x0;   // the left boundary of the domain
    double m_xend; // the right boundary of the domain
    int m_nPoints; // the number of grid-points used in the finite difference
                   // method
    double m_dx;   // the distance between 2 grid points

    // Temporal discretisation:
    double m_t0;      // the initial time
    double m_tend;    // the end time
    int m_nTimeSteps; // the number of time-steps
    double m_dt;      // the size of a time-step

    // Min and max values in the solutions
    double m_minValue;
    double m_maxValue;

}; // end class DemoSolver

///////////////////////////////////////////////////////////////////////////////
// Class that represents the 1D finite difference solver
class OneDFiniteDiffAdvDiffSolver : public DemoSolver
{
public:
    // constructor based upon the discretisation details
    OneDFiniteDiffAdvDiffSolver(int nVars, int nPoints, int nTimeSteps) :
        DemoSolver(nVars, nPoints, nTimeSteps),
        m_wavenumber(1.0), m_U(1.0), m_D(0.05)
    {
        m_fileName = std::string("OneDFiniteDiffAdvDiffSolver");
        m_title = std::string("Finite Difference Solution to the 1D advection-diffusion equation");
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void HelmSolve(const Array<OneD, const Array<OneD, double>> &inarray,
                         Array<OneD,       Array<OneD, double>> &outarray,
                   const NekDouble time, const NekDouble lambda) const;

    void EvaluateAdvectionTerm(
        const Array<OneD, const Array<OneD, double>> &inarray,
              Array<OneD,       Array<OneD, double>> &outarray,
        const NekDouble time) const;

    // -----------------------------------------------------------------
    void EvaluateExactSolution(Array<OneD, Array<OneD, double>> &outarray,
                               const NekDouble time) const;

    // -----------------------------------------------------------------

private:
    // Value of the coefficients:
    double m_wavenumber; // wave number
    double m_U;          // advection speed
    double m_D;          // diffusion coefficient

    void solveTriDiagMatrix(int n, double alpha, double beta,
                            const Array<OneD, const double> &inarray,
                            Array<OneD, double> &outarray) const;

}; // end class OneDFiniteDiffAdvDiffSolver

///////////////////////////////////////////////////////////////////////////////
// Class that represents the 1D sinusoid solver
class OneDSinusoidSolver : public DemoSolver
{
public:
    // constructor based upon the discretisation details
    OneDSinusoidSolver(int nVars, int nPoints, int nTimeSteps) :
        DemoSolver(nVars, nPoints, nTimeSteps)
    {
        m_fileName = std::string("OneDSinusoidSolver");
        m_title = std::string("Solution to the 1D Sinusoid equation");

        // Frequencies and phases for the sinusoidal solver.
        m_freqs  = Array<OneD, double>(m_nVars, 0.0 );
        m_phases = Array<OneD, double>(m_nVars, 0.0 );

        // Assumption: the two-dimensional Lambda matrix is a diagonal
        // matrix thus values are non zero if and only i=j. As such,
        // the diagonal Lambda values are stored as two vectors so to
        // accomodate complex numbers lambda[0] real, lambda[1]
        // imaginary.
        m_A      = Array<TwoD, NekDouble>(m_nVars, m_nVars, 0.0);
        m_lambda = Array<OneD, std::complex<NekDouble>>(m_nVars, 0.0);

        // Initial values.
        m_z0 = Array<OneD, NekDouble>(m_nVars, 0.0);

        // Initialize a random seed using the time.
        srand( time(NULL) );

        for (int k = 0; k < m_nVars; k++)
        {
            m_freqs [k] = (double) rand() / (double) RAND_MAX;
            m_phases[k] = (double) rand() / (double) RAND_MAX;

            m_lambda[k] =
              std::complex< NekDouble>((double) rand() / (double) RAND_MAX, 0);

            m_z0[k] = (double) rand() / (double) RAND_MAX;
        }

        // Fixed values for now based on a known solution for 5 values.
        static double freqs[]  = { 0.1405, 1.0930, -0.3592, -1.3925, 0.3480 };
        static double phases[] = { -0.7197, 1.7775, -0.8679, 0.1710, 2.0424 };
        static double lambda[] = { 0.7405, -0.2949, -1.6821, -1.1045, -1.0519 };
        static double y0[]     = { 1.0541, -0.0169, 0.4792, 2.4027, 1.4006 };

        for (int k = 0; k < m_nVars; k++)
        {
            m_freqs [k] = freqs[k];
            m_phases[k] = phases[k];
            m_lambda[k] = std::complex<NekDouble>(lambda[k], 0);
            m_z0[k] = y0[k];
        }
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void EvaluateSinusoidTerm(
        const Array<OneD, const Array<OneD, double>> &inarray,
              Array<OneD,       Array<OneD, double>> &outarray,
        const NekDouble time) const;

    // -----------------------------------------------------------------
    void EvaluateExactSolution(Array<OneD, Array<OneD, double>> &outarray,
                               const NekDouble time) const;
    // -----------------------------------------------------------------

    Array<OneD, std::complex<NekDouble>> &GetLambda()
    {
        return m_lambda;
    }

private:
    // Value of the coefficients:
    Array<OneD, double> m_freqs;
    Array<OneD, double> m_phases;

    Array<TwoD, NekDouble> m_A;
    Array<OneD, std::complex<NekDouble>> m_lambda;

    Array<OneD, NekDouble> m_z0;

}; // end class OneDSinusoidSolver

///////////////////////////////////////////////////////////////////////////////
// Class that represents the 1D FDE solver
class OneDFDESolver : public DemoSolver
{
public:
    // constructor based upon the discretisation details
    OneDFDESolver(int nVars, int nPoints, int nTimeSteps) :
        DemoSolver(nVars, nPoints, nTimeSteps)
    {
        m_fileName = std::string("OneDFDESolver");
        m_title = std::string("Solution to the 1D constant equation");

        m_alpha = 0.3;

        // Initial values set to zero
        m_u0 = Array<OneD, NekDouble>(m_nVars, 0.0);
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void EvaluateFDETerm(
        const Array<OneD, const Array<OneD, double>> &inarray,
              Array<OneD,       Array<OneD, double>> &outarray,
        const NekDouble time) const;

    // -----------------------------------------------------------------
    void EvaluateExactSolution(Array<OneD, Array<OneD, double>> &outarray,
                               const NekDouble time) const;
    // -----------------------------------------------------------------

private:
    // Value of the coefficients:
    NekDouble m_alpha; // Fractional order

    Array<OneD, NekDouble> m_u0;

}; // end class OneDFDESolver

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    // So to be case insensitive for the --L2 or --l2 option.
    namespace po_style = boost::program_options::command_line_style;

    po::options_description desc("Usage:");

    desc.add_options()
      ("help,h", "Produce this help message.")
      ("verbose,v", "Print the solution values for each time step.")
      ("L2,l", "Print the L2 error for each time step.")
      ("dof,d", po::value<int>(), "Number of degrees of freedom (points or values).")
      ("timesteps,t", po::value<int>(), "Number of timesteps.")
      ("order,o", po::value<int>(), "Order of the scheme.")
      ("parameter,p", po::value<std::string>(), "Free parameters for the scheme.")
      ("method,m", po::value<int>(),
        "Number for the time-integration scheme:\n"
        "- 0: 1st order Forward Euler\n"
        "- 1: 1st order multi-step IMEX scheme\n"
        "     (Euler Backwards/Euler Forwards)\n"
        "- 2: 2nd order multi-step IMEX scheme\n"
        "- 3: 3rd order multi-step IMEX scheme\n"
        "- 4: 4th order multi-step IMEX scheme\n"
        "  \n"
        "- 5: 2nd order multi-stage DIRK IMEX scheme\n"
        "- 6: 3nd order multi-stage DIRK IMEX scheme\n"
        "- 7: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
        "- 8: 2nd order Crank-Nicolson/Adams-Bashforth (CNAB)\n"
        "- 9: 2nd order Modified Crank-Nicolson/Adams-Bashforth\n"
        "     (MCNAB)\n"
        "  \n"
        "- 10: Nth order multi-stage Runga-Kutta scheme\n"
        "- 11: Nth order multi-stage Runga-Kutta SSP scheme\n"
        "- 12: Nth order multi-stage Diagonally Implicit\n"
        "      Runga-Kutta scheme (DIRK)\n"
        "- 13: Nth order multi-step Adams-Bashforth scheme\n"
        "- 14: Nth order multi-step Adams-Moulton scheme\n"
        "- 15: Nth order multi-step BDFImplicit scheme\n"
        "- 16: Nth order multi-step IMEX scheme\n"
        "      (Euler Backwards/Euler Forwards)\n"
        "- 17: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
        "- 18: Nth order multi-stage IMEX DIRK scheme\n"
        "  \n"
        "- 20: Nth order multi-step Lawson-Euler exponential scheme\n"
        "- 21: Nth order multi-step Norsett-Euler exponential scheme\n"
        "  \n"
        // "- 25: Nth order multi-step Fractional In Time scheme\n"
       );

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc)
                  .style(po_style::unix_style|po_style::case_insensitive)
                  .run(), vm);
        po::notify(vm);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl
                  << desc;
        return 1;
    }

    if (!vm.count("dof") || !vm.count("timesteps") || !vm.count("method") ||
        ((vm["method"].as<int>() >= 10) && !vm.count("order")) ||
        ((vm["method"].as<int>() == 18) && !vm.count("parameter")) ||
        vm.count("help"))
    {
        std::cout << "Please specify the number of "
                  << "dof and timesteps along with the order and method."
                  << std::endl << std::endl
                  << desc;
        return 1;
    }

    bool L2      = vm.count("L2");
    bool verbose = vm.count("verbose");

    int nDoF    = vm["dof"].as<int>();
    int nTimeSteps = vm["timesteps"].as<int>();
    int nMethod    = vm["method"].as<int>();
    int nOrder     = 0;
    std::vector<NekDouble> freeParams;

    if((vm["method"].as<int>() >= 10) )
    {
        nOrder = vm["order"].as<int>();
    }

    // IMEX DIRK methods also require free parameters.
    if((vm["method"].as<int>() == 18) )
    {
        std::string sParameter = vm["parameter"].as<std::string>();

        // Parse the free parameters which are in a string.
        while( sParameter.size() )
        {
            size_t found = sParameter.find(" ");

            if( found == std::string::npos )
            {
                int fp = stoi( sParameter );

                freeParams.push_back(fp);

                sParameter = "";
            }
            else if( found > 0 )
            {
                int fp = stoi( sParameter.substr(0, found) );

                freeParams.push_back(fp);

                sParameter = sParameter.substr(found+1);
            }
        }
    }

    if( nDoF < 2 )
    {
        std::cout << "Please specify the number of "
                  << "dof to be greater than 1."
                  << std::endl << std::endl
                  << desc;
        return 1;
    }

    // -------------------------------------------------------------------------
    // The implementation below shows how the time-stepping framework can be
    // used for time-integration PDEs

    // 1. THE TEMPORAL DISCRETISATION
    //    Create time integrator requested.
    TimeIntegrationSchemeFactory &factory =
        LibUtilities::GetTimeIntegrationSchemeFactory();

    LibUtilities::TimeIntegrationSchemeSharedPtr tiScheme;

    switch (nMethod)
    {
        case 0:
            tiScheme = factory.CreateInstance("ForwardEuler", "", 1, freeParams);
            break;
        case 1:
            tiScheme = factory.CreateInstance("IMEXOrder1", "", 1, freeParams);
            break;
        case 2:
            tiScheme = factory.CreateInstance("IMEXOrder2", "", 2, freeParams);
            break;
        case 3:
            tiScheme = factory.CreateInstance("IMEXOrder3", "", 3, freeParams);
            break;
        case 4:
            tiScheme = factory.CreateInstance("IMEXOrder4", "", 4, freeParams);
            break;
        case 5:
            tiScheme = factory.CreateInstance("IMEXdirk_2_3_2", "", 2, freeParams);
            break;
        case 6:
            tiScheme = factory.CreateInstance("IMEXdirk_3_4_3", "", 3, freeParams);
            break;
        case 7:
            tiScheme = factory.CreateInstance("IMEX", "Gear", 2, freeParams);
            break;
        case 8:
            tiScheme = factory.CreateInstance("CNAB", "", 2, freeParams);
            break;
        case 9:
            tiScheme = factory.CreateInstance("MCNAB", "", 2, freeParams);
            break;
        case 10:
            tiScheme = factory.CreateInstance("RungeKutta", "", nOrder, freeParams);
            break;
        case 11:
            tiScheme = factory.CreateInstance("RungeKutta", "SSP", nOrder, freeParams);
            break;
        case 12:
            tiScheme = factory.CreateInstance("DIRK", "", nOrder, freeParams);
            break;
        case 13:
            tiScheme = factory.CreateInstance("AdamsBashforth", "", nOrder, freeParams);
            break;
        case 14:
            tiScheme = factory.CreateInstance("AdamsMoulton", "", nOrder, freeParams);
            break;
        case 15:
            tiScheme = factory.CreateInstance("BDFImplicit", "", nOrder, freeParams);
            break;
        case 16:
            tiScheme = factory.CreateInstance("IMEX", "", nOrder, freeParams);
            break;
        case 17:
            tiScheme = factory.CreateInstance("IMEX", "Gear", 2, freeParams);
            break;
        case 18:
            tiScheme = factory.CreateInstance("IMEX", "dirk", nOrder, freeParams);
            break;

        case 20:
            tiScheme = factory.CreateInstance("EulerExponential", "Lawson", nOrder, freeParams );
            break;
        case 21:
            tiScheme = factory.CreateInstance("EulerExponential", "Norsett", nOrder, freeParams);
            break;
        // case 25:
        //     tiScheme = factory.CreateInstance("FractionalInTime", "", nOrder, freeParams);
        //     break;
        default:
        {
            std::cerr << "The -m argument defines the "
                      << "time-integration method to be used:"
                      << std::endl << std::endl;
            std::cout << desc;
            exit(1);
        }
    }

    int nVariables;
    int nPoints;

    // 2. THE SPATIAL DISCRETISATION
    //    Create an object of the DemoSolver class.
    //    This class can be thought of as representing the spatial (finite
    //    difference) discretisation.
    LibUtilities::TimeIntegrationSchemeOperators ode;

    std::shared_ptr<DemoSolver> solverSharedPtr;

    if( tiScheme->GetIntegrationSchemeType() == eFractionalInTime )
    {
        if( nDoF > 2 )
        {
          std::cout << "Curently set up for 2 variables, setting nDoF to 2"
                    << std::endl;

          nDoF = 2;
        }

        nVariables = nDoF;
        nPoints    = 1;

        OneDFDESolver *tmpSolver =
          new OneDFDESolver(nVariables, nPoints, nTimeSteps);

        ode.DefineOdeRhs(&OneDFDESolver::EvaluateFDETerm, tmpSolver);

        solverSharedPtr = std::shared_ptr<DemoSolver> (tmpSolver);
    }
    else if( tiScheme->GetIntegrationSchemeType() == eExponential )
    {
        if( nDoF > 5 )
        {
          std::cout << "Curently set up for 5 variables, setting nDoF to 5"
                    << std::endl;

          nDoF = 5;
        }

        nVariables = nDoF;
        nPoints    = 1;

        OneDSinusoidSolver *tmpSolver =
          new OneDSinusoidSolver(nVariables, nPoints, nTimeSteps);

        ode.DefineOdeRhs(&OneDSinusoidSolver::EvaluateSinusoidTerm, tmpSolver);

        solverSharedPtr = std::shared_ptr<DemoSolver> (tmpSolver);

        // For exponential integrators, the coefficents for each
        // variable needs to be set.
        tiScheme->SetExponentialCoefficients(tmpSolver->GetLambda());
    }
    else
    {
        nVariables = 1;
        nPoints    = nDoF;

        OneDFiniteDiffAdvDiffSolver *tmpSolver =
          new OneDFiniteDiffAdvDiffSolver(nVariables, nPoints, nTimeSteps);

        // After this spatial discretisation, the PDE has actually
        // been reduced (through the method-of-lines) to an ODE. In
        // order to use the time-stepping framework, we need to give
        // it the necessary information about this ODE. Therefore, we
        // create an oject of the class TimeIntegrationSchemeOperators
        // that contains a
        // 'function pointer' (in fact a 'functor') to the
        // - explicit term of the ODE (i.e. the advection term)
        // - implicit solve routine (i.e. the Helmholtz solver)
        // - projection operator (i.e. the identity operator in this case)
        ode.DefineOdeRhs(&OneDFiniteDiffAdvDiffSolver::EvaluateAdvectionTerm,
                         tmpSolver);
        ode.DefineImplicitSolve(&OneDFiniteDiffAdvDiffSolver::HelmSolve,
                                tmpSolver);

        solverSharedPtr = std::shared_ptr<DemoSolver> (tmpSolver);
    }

    ode.DefineProjection(&DemoSolver::Project, solverSharedPtr);

    // 3. Initialise the arrays that contain the approximate and exact
    // solutions.

    // Note, approxSol and exactSol are an array of arrays. The rows
    // contain the variables while the colums contains the points.

    // For exponential solvers the number of variables > 1 and the
    // number of points is 1, for all other solvers the number of
    // variables = 1 and the number of points > 1.

    // Array containing the approximate solution
    Array<OneD, Array<OneD, double>> approxSol(nVariables);
    // Array containing the exact solution
    Array<OneD, Array<OneD, double>>  exactSol(nVariables);

    for( int k=0; k<nVariables; ++k )
    {
        approxSol[k] = Array<OneD, double>(nPoints);
        exactSol[k]  = Array<OneD, double>(nPoints);
    }

    // Initialize both solutions using the exact solution.
    double t0 = solverSharedPtr->GetInitialTime();

    solverSharedPtr->EvaluateExactSolution(approxSol, t0);
    solverSharedPtr->EvaluateExactSolution(exactSol,  t0);

    // 3.1 Initialize the time-integration scheme.
    double dt = solverSharedPtr->GetDeltaT();

    LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr
        solutionPtr;

    // The Fractional in Time needs the end time whereas the GLMs need
    // the initial time.
    if( tiScheme->GetIntegrationSchemeType() == eFractionalInTime )
    {
        double t1 = solverSharedPtr->GetFinialTime();
        solutionPtr = tiScheme->InitializeScheme(dt, approxSol, t1, ode);
    }
    else
    {
        solutionPtr = tiScheme->InitializeScheme(dt, approxSol, t0, ode);
    }

    // 4. Open a file for writing the solution
    std::ofstream outfile;
    outfile.open( solverSharedPtr->GetFileName() + ".dat");

    // Save the scheme name for outputting.
    solverSharedPtr->SetSchemeName( tiScheme->GetName() );

    // Write the initial conditions to the output file
    solverSharedPtr->AppendOutput( outfile, 0, 0, approxSol, exactSol);

    int timeStep = 0;
    double time = t0;

    // 5. Do the time integration.
    while( timeStep < nTimeSteps )
    {
        // Time integration for one time step
        approxSol = tiScheme->TimeIntegrate(timeStep, dt, solutionPtr, ode);

        // Advance the time for the exact solution.
        time += dt;

        // Calculate the exact solution
        solverSharedPtr->EvaluateExactSolution(exactSol, time);

        // At this point the time step is finished so increment the time step.
        ++timeStep;

        // Save the solutions the output file
        solverSharedPtr->AppendOutput(outfile, timeStep, time, approxSol, exactSol);

        // 6. Calculate the min / max values. If the last argument is
        // true the values will be dumped to screen.

        // Write the time step and time
        if( verbose || L2 )
          std::cout << "Time step: " << timeStep << "  "
                    << "Time: " << time << std::endl;

        solverSharedPtr->GetMinMaxValues(approxSol, exactSol, verbose);

        if( L2 )
          solverSharedPtr->EvaluateL2Error(approxSol, exactSol);
    }

    // 7. Calculate the error for the last time step and dump to screen
    std::cout << tiScheme->GetFullName() << std::endl;
    std::cout << "Time step: " << timeStep << "  "
              << "Time: " << time << std::endl;

    solverSharedPtr->EvaluateL2Error(approxSol, exactSol);

    // 7. Some more writing out the results
    outfile.close();

    solverSharedPtr->GenerateGnuplotScript(tiScheme->GetFullName());

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void DemoSolver::Project(
    const Array<OneD, const Array<OneD, double>> &inarray,
          Array<OneD,       Array<OneD, double>> &outarray,
    const NekDouble time) const
{
    boost::ignore_unused(time);

    // This is simply the identity operator.
    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            outarray[k][i] = inarray[k][i];
        }
    }
}

// Calculate the Relative Error L2 Norm (as opposed to the absolute L2
// norm).
void DemoSolver::GetMinMaxValues(
    const Array<OneD, const Array<OneD, double>> &approx,
    const Array<OneD, const Array<OneD, double>> &exact,
    bool print )
{
    // Get the min and max value and write the approximate solution
    if( print )
        std::cout << "  approximate ";

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            // Do not get the min/max values becasue the approximate
            // solution can blow up.
            // if( m_minValue > approx[k][i] )
            //     m_minValue = approx[k][i];

            // if( m_maxValue < approx[k][i] )
            //     m_maxValue = approx[k][i];

          if( print )
              std::cout << approx[k][i] << "  ";
        }
    }

    if( print )
        std::cout << std::endl;

    // Get the min and max value and write the exact solution
    if( print )
        std::cout << "  exact       ";

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            if( m_minValue > exact[k][i] )
                m_minValue = exact[k][i];

            if( m_maxValue < exact[k][i] )
                m_maxValue = exact[k][i];

            if( print )
                std::cout << exact[k][i] << "  ";
        }
    }

    if( print )
        std::cout << std::endl;
}

// Calculate the Relative Error L2 Norm (as opposed to the absolute L2
// norm).
void DemoSolver::EvaluateL2Error(
    const Array<OneD, const Array<OneD, double>> &approx,
    const Array<OneD, const Array<OneD, double>> &exact )
{
    // Calcualate the sum of squares for the L2 Norm.
    double a = 0.0;
    double b = 0.0;

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            a += (approx[k][i] - exact[k][i]) * (approx[k][i] - exact[k][i]);
            b += exact[k][i] * exact[k][i];
        }
    }

    ASSERTL1( b>DBL_EPSILON,
              "Exact solution is near zero. L2 Norn is invalid" );

    // Calculate the relative error L2 Norm.
    double norm = sqrt(a / b);

    std::cout << "L 2 error :" << norm << std::endl;
}

void DemoSolver::AppendOutput(
    std::ofstream &outfile, int timeStepNumber, const NekDouble time,
    const Array<OneD, const Array<OneD, double>> &approx,
    const Array<OneD, const Array<OneD, double>> &exact) const
{
    if (timeStepNumber == 0)
    {
        // Save some data provenance and other useful info in output file...
        outfile << "# Data in this file consists of " << m_nTimeSteps
                << " time steps (there will be" << std::endl
                << "# a blank line between each data set for each time step)."
                << std::endl
                << "#" << std::endl
                << "# Delta T: " << m_dt
                << std::endl
                << "# Method:  " << m_schemeName
                << std::endl
                << "#" << std::endl
                << "# There are 3 columns with the following headers:"
                << std::endl
                << "#" << std::endl;

        if( m_nVars > 1 )
        {
            outfile << "#   Variable";
        }
        else if( m_nPoints > 1 )
        {
            outfile << "#   Location";
        }

        outfile << "  |  Approximate Solution  |  Exact Solution"
                << std::endl
                << "#"
                << std::endl << std::endl
                << "# Initial condition (at time " << time << "):"
                << std::endl;
    }
    else
    {
        outfile << "# Time step: " << timeStepNumber << ", time: " << time
                << std::endl;    }

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            if( m_nVars > 1 )
            {
                outfile << std::scientific
                        << std::setw(17) << std::setprecision(10)
                        << k
                        << "  " << approx[k][i] << "  " << exact[k][i]
                        << std::endl;
            }
            else if( m_nPoints > 1 )
            {
                outfile << std::scientific
                        << std::setw(17) << std::setprecision(10)
                        << m_x0 + i * m_dx
                        << "  " << approx[k][i] << "  " << exact[k][i]
                        << std::endl;
            }
        }
    }

    // Gnuplot uses two blank lines between each set of data.
    outfile << std::endl << std::endl;
}

void DemoSolver::GenerateGnuplotScript(const std::string &method) const
{
    std::ofstream outfile;
    outfile.open(m_fileName + ".p");
    outfile << "# Gnuplot script file" << std::endl
            << "set   autoscale" << std::endl
            << "unset log" << std::endl
            << "unset label" << std::endl
            << "set xtic auto" << std::endl
            << "set ytic auto" << std::endl
      // FIXME
      // << "set title 'Finite Difference Solution to the 1D "
      // << "advection-diffusion equation "
            << "set title 'Approximate vs exact solution using method "
            << method << "'" << std::endl
            << "set xlabel 'x'" << std::endl
            << "set ylabel 'u'" << std::endl;

    if( m_nVars > 1 )
    {
        outfile << "set xr [" << 0          << ":" << m_nVars-1  << "]"
                << std::endl
                << "set yr [" << m_minValue << ":" << m_maxValue << "]"
                << std::endl;
    }
    else if( m_nPoints > 1 )
    {
        outfile << "set xr [" << m_x0       << ":" << m_xend     << "]"
                << std::endl
                << "set yr [" << m_minValue << ":" << m_maxValue << "]"
                << std::endl;
    }

    for (int i = 0; i <= m_nTimeSteps; i++)
    {
        double t = m_t0 + (i * m_dt);

        outfile << "plot '" << m_fileName << ".dat' using 1:2 index "
                << i << " title 'Approximate Solution (t=" << t
                << ")' with linespoints , "

                << "'" << m_fileName << ".dat' using 1:3 index "
                << i << " title 'Exact Solution (t=" << t
                << ")' with linespoints " << std::endl

                << "pause " << 4.0 / m_nTimeSteps << std::endl;
    }

    // Keep window open until the user clicks
    outfile << "pause mouse any" << std::endl;

    outfile.close();
}


///////////////////////////////////////////////////////////////////////////////
void OneDFiniteDiffAdvDiffSolver::HelmSolve(
    const Array<OneD, const Array<OneD, double>> &inarray,
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time,
    const NekDouble lambda) const
{
    boost::ignore_unused(time);

    // This function implements a 1D finite difference helmholtz solver.
    // The 1D Helmholtz equation leads to a cyclic triadiagonal matrix to be
    // solved.
    // In this function, this is solved as follows:
    // - Based upon the static condensation technique, we first solve
    //   for the periodic grid-points (i.e. the boundary nodes)
    // - Next, we solve for the interior grid-points. The associated tridiagonal
    //   system is solved based upon the Thomas algorithm
    double a = -m_D * lambda / (m_dx * m_dx); // off diagonal term
    double b = 1.0 +
               2.0 * lambda * m_D /
                   (m_dx * m_dx); // diagonal term of triadiagonal matrix

    int nIntPoints = m_nPoints - 2;

    for (int k = 0; k < m_nVars; k++)
    {
        Array<OneD, double> invD_f(nIntPoints);
        solveTriDiagMatrix(nIntPoints, a, b, inarray[k] + 1, invD_f);

        Array<OneD, double> C(nIntPoints, 0.0);
        Array<OneD, double> invD_C(nIntPoints, 0.0);
        C[0]              = a;
        C[nIntPoints - 1] = a;
        solveTriDiagMatrix(nIntPoints, a, b, C, invD_C);

        outarray[k][0] =
          (inarray[k][0] - a * (invD_f[0] + invD_f[nIntPoints - 1])) /
          (b - a * (invD_C[0] + invD_C[nIntPoints - 1]));
        outarray[k][m_nPoints - 1] = outarray[k][0];

        Array<OneD, double> f(nIntPoints);
        for (int i = 0; i < nIntPoints; i++)
          {
            f[i] = inarray[k][i + 1];
          }
        f[0] -= outarray[k][0] * a;
        f[nIntPoints - 1] -= outarray[k][0] * a;

        Array<OneD, double> tmp;
        solveTriDiagMatrix(nIntPoints, a, b, f,
                           tmp = outarray[k] + 1); // Calls the Thomas algorithm
    }
}

void OneDFiniteDiffAdvDiffSolver::EvaluateAdvectionTerm(
    const Array<OneD, const Array<OneD, double>> &inarray,
          Array<OneD,       Array<OneD, double>> &outarray,
    const NekDouble time) const
{
    boost::ignore_unused(time);

    for (int k = 0; k < m_nVars; k++)
    {
        // The advection term can be evaluated using central or upwind
        // differences
        if (true)
        {
            // Note: We are using a periodic boundary condition where
            // the 1st point and last point are actually the same
            // point.  This is why the 1st (and last) point in the
            // output array (index 0 and m_nPoints-1 respectively) are
            // NOT used for the central differences, and instead the
            // 2nd point (index 1) and 2nd to last point are used.

            // Central differences:
            outarray[k][0] =
              -m_U * (inarray[k][1] - inarray[k][m_nPoints - 2]) / (2.0 * m_dx);
            outarray[k][m_nPoints - 1] = outarray[k][0];

            for (int i = 1; i < m_nPoints - 1; i++)
            {
                outarray[k][i] =
                  -m_U * (inarray[k][i + 1] - inarray[k][i - 1]) / (2.0 * m_dx);
            }
        }
        else
        {
            // upwind differences
            for (int i = 1; i < m_nPoints; i++)
            {
                outarray[k][i] =
                  -m_U * (inarray[k][i] - inarray[k][i - 1]) / (m_dx);
            }

            outarray[k][0] = outarray[k][m_nPoints - 1];
        }
    }
}

void OneDFiniteDiffAdvDiffSolver::solveTriDiagMatrix(
    int n, double a, double b,
    const Array<OneD, const double> &inarray,
          Array<OneD,       double> &outarray) const
{
    // Implementation of the Thomas algorithm for Tridiaginol systems
    Array<OneD, double> cprime(n);
    Array<OneD, double> dprime(n);

    cprime[0] = a / b;
    dprime[0] = inarray[0] / b;

    for (int i = 1; i < n; i++)
    {
        double id = 1.0 / (b - cprime[i - 1] * a);
        cprime[i] = a * id;
        dprime[i] = (inarray[i] - dprime[i - 1] * a) * id;
    }

    outarray[n - 1] = dprime[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        outarray[i] = dprime[i] - cprime[i] * outarray[i + 1];
    }
}

void OneDFiniteDiffAdvDiffSolver::EvaluateExactSolution(
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            double x       = m_x0 + i * m_dx;
            double wn = 2.0 * M_PI * m_wavenumber;
            outarray[k][i] =
              exp(-m_D * wn * wn * time) * sin(wn * (x - m_U * time));

            outarray[k][i] = exp(-m_D * 2.0 * 2.0 * M_PI * M_PI * m_wavenumber *
                                 m_wavenumber * time) *
              sin(2.0 * m_wavenumber * M_PI * (x - m_U * time));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void OneDSinusoidSolver::EvaluateSinusoidTerm(
    const Array<OneD, const Array<OneD, double>> &inarray,
          Array<OneD,       Array<OneD, double>> &outarray,
    const NekDouble time) const
{
    boost::ignore_unused(inarray);

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            outarray[k][i] = sin(m_freqs[k]*time + m_phases[k]);
        }
    }
}

void OneDSinusoidSolver::EvaluateExactSolution(
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    if( time == GetInitialTime() )
    {
        for (int k = 0; k < m_nVars; k++)
        {
            for (int i = 0; i < m_nPoints; i++)
            {
                outarray[k][i] = m_z0[k];
            }
        }
    }
    else
    {
        // Computes the exact solution at time t to the ODE
        //
        //   y' = A*y + sin(freq*t + phi),
        //   y(0) = y0
        //
        // where y, freq, and phi are vectors and A is a matrix. The
        // matrix A has the eigenvalue decomposition
        //
        //   A = V * diag(lambda) * inv(V)

        // ARS Note - the above is not yet true.

        // Compute right-hand side of the normal equations
        for (int k = 0; k < m_nVars; k++)
        {
            for (int i = 0; i < m_nPoints; i++)
            {
                double v = cos( m_phases[k] );
                double w = sin( m_phases[k] );

                double sinft = sin(m_freqs[k]*time);
                double cosft = cos(m_freqs[k]*time);

                double lambdaFreq2 =
                  m_freqs[k]*m_freqs[k] + m_lambda[k].real()*m_lambda[k].real();

                outarray[k][i] =
                  // exp(lambda*T) term
                  (std::exp(m_lambda[k].real() * time) *
                   (m_z0[k] + (m_lambda[k].real() * w + m_freqs[k] * v) / lambdaFreq2)) +

                  // sin(f*T) term
                  ((( m_freqs[k] * w - m_lambda[k].real() * v) / lambdaFreq2) * sinft) +

                  // cos(f*T) term
                  (((-m_lambda[k].real() * w - m_freqs[k] * v) / lambdaFreq2) * cosft);
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
void OneDFDESolver::EvaluateFDETerm(
    const Array<OneD, const Array<OneD, double>> &inarray,
          Array<OneD,       Array<OneD, double>> &outarray,
    const NekDouble time) const
{
    boost::ignore_unused(inarray, time);

    for (int i = 0; i < m_nVars; i++)
    {
        for (int j = 0; j < m_nPoints; j++)
        {
            if( i == 0 )
                outarray[i][j] =       tgamma( m_alpha + 1.0 );
            else if( i == 1 )
                outarray[i][j] = 0.5 * tgamma( m_alpha + 1.0 );
        }
    }
}

void OneDFDESolver::EvaluateExactSolution(
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    if( time == GetInitialTime() )
    {
        for (int i = 0; i < m_nVars; i++)
        {
            for (int j = 0; j < m_nPoints; j++)
            {
                outarray[i][j] = m_u0[i];
            }
        }
    }
    else
    {
        // Computes the exact solution at time t to the ODE
        //
        //   y' = time^alpha
        //   y(0) = y0
        //

        // Compute right-hand side of the normal equations
        for (int i = 0; i < m_nVars; i++)
        {
            for (int j = 0; j < m_nPoints; j++)
            {
                if( i == 0 )
                    outarray[i][j] =       std::pow( time, m_alpha );
                else if( i == 1 )
                    outarray[i][j] = 0.5 * std::pow( time, m_alpha );
            }
        }
    }
}
