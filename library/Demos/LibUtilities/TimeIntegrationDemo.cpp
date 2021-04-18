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
//  |              V = 1
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
//   - OneDFiniteDiffAdvDiffSolverError.dat (containing the error data)
//   - OneDFiniteDiffAdvDiffSolverError.p   (containing a gnuplot script)
//
//   - OneDFiniteDiffAdvDiffSolverL2Norm.dat (containing the l2 Norm data)
//   - OneDFiniteDiffAdvDiffSolverL2Norm.p   (containing a gnuplot script)
//
// and can be visualised by gnuplot using the command
//
//    gnuplot OneDFiniteDiffAdvDiffSolverOutput.p
//    gnuplot OneDFiniteDiffAdvDiffSolverOutputError.p
//    gnuplot OneDFiniteDiffAdvDiffSolverOutputL2Norm.p
//
// -----------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/program_options.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/TimeIntegration/EulerExponentialTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>

#include <time.h> /* time */

using namespace Nektar;
using namespace Nektar::LibUtilities;

namespace po = boost::program_options;

// Base class for the solver.
class DemoSolver
{
public:
    // -----------------------------------------------------------------
    // Constructor based upon the discretisation details
    DemoSolver(int nVars, int nPoints, int nTimeSteps, bool test)
        : m_nVars(nVars),

          m_x0(0.0), m_xend(1.0), m_nPoints(nPoints),
          m_dx((m_xend - m_x0) / ((double)m_nPoints - 1.0)),

          m_t0(0.0), m_tend(1.0), m_nTimeSteps(nTimeSteps),
          m_dt((m_tend - m_t0) / (double)m_nTimeSteps)
    {
        boost::ignore_unused(test);

        m_minValue = +std::numeric_limits<double>::max();
        m_maxValue = -std::numeric_limits<double>::max();

        m_minError = +std::numeric_limits<double>::max();
        m_maxError = -std::numeric_limits<double>::max();

        m_minL2Norm = +std::numeric_limits<double>::max();
        m_maxL2Norm = -std::numeric_limits<double>::max();
    }

    virtual ~DemoSolver(){};

    // -----------------------------------------------------------------
    // Exact solution and project (identity)
    virtual void EvaluateExactSolution(
        Array<OneD, Array<OneD, double>> &outarray,
        const NekDouble time) const = 0;

    virtual void Project(const Array<OneD, const Array<OneD, double>> &inarray,
                         Array<OneD, Array<OneD, double>> &outarray,
                         const NekDouble time) const;

    // -----------------------------------------------------------------
    // Misc functions for error and outputing
    void GetMinMaxValues(const Array<OneD, const Array<OneD, double>> &exact,
                         const Array<OneD, const Array<OneD, double>> &approx,
                         bool print);

    double EvaluateL2Error(const Array<OneD, const Array<OneD, double>> &exact,
                           const Array<OneD, const Array<OneD, double>> &approx,
                            bool print);

    void AppendOutput(std::ofstream &outfile, std::ofstream &errorfile,
                      std::ofstream &L2normfile, const int timeStepNumber,
                      const NekDouble time,
                      const Array<OneD, const Array<OneD, double>> &exact,
                      const Array<OneD, const Array<OneD, double>> &approx,
                      const double L2Norm) const;

    void GenerateGnuplotScript(const std::string &method) const;

    // -----------------------------------------------------------------
    // Access methods
    double GetInitialTime() const
    {
        return m_t0;
    }

    double GetFinalTime() const
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

    void SetSchemeName(std::string name)
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

    // Min and max errors in the solutions
    double m_minError;
    double m_maxError;

    // Min and max L2Norms in the solutions
    double m_minL2Norm;
    double m_maxL2Norm;

}; // end class DemoSolver

///////////////////////////////////////////////////////////////////////////////
// Class that represents the 1D finite difference solver
class OneDFiniteDiffAdvDiffSolver : public DemoSolver
{
public:
    // constructor based upon the discretisation details
    OneDFiniteDiffAdvDiffSolver(int nVars, int nPoints, int nTimeSteps,
                                bool test)
        : DemoSolver(nVars, nPoints, nTimeSteps, test), m_wavenumber(1.0),
          m_V(1.0), m_D(0.05)
    {
        m_fileName = std::string("OneDFiniteDiffAdvDiffSolver");
        m_title    = std::string("Finite Difference Solution to the 1D "
                              "advection-diffusion equation");
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void HelmSolve(const Array<OneD, const Array<OneD, double>> &inarray,
                   Array<OneD, Array<OneD, double>> &outarray,
                   const NekDouble time, const NekDouble lambda) const;

    void EvaluateAdvectionTerm(
        const Array<OneD, const Array<OneD, double>> &inarray,
        Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const;

    // -----------------------------------------------------------------
    void EvaluateExactSolution(Array<OneD, Array<OneD, double>> &outarray,
                               const NekDouble time) const;

    // -----------------------------------------------------------------

private:
    // Value of the coefficients:
    double m_wavenumber; // wave number
    double m_V;          // advection speed
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
    OneDSinusoidSolver(int nVars, int nPoints, int nTimeSteps, bool test)
        : DemoSolver(nVars, nPoints, nTimeSteps, test)
    {
        m_fileName = std::string("OneDSinusoidSolver");
        m_title    = std::string("Solution to the 1D Sinusoid equation");

        // Frequencies and phases for the sinusoidal solver.
        m_freqs  = Array<OneD, double>(m_nVars, 0.0);
        m_phases = Array<OneD, double>(m_nVars, 0.0);

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
        if (test)
            srand(0);
        else
            srand(time(NULL));

        // Randomly generate the jacobian in a way that essentially
        // ensures diagonalizability with real eigenvalues and
        // invertibility. Real eigenvalues aren't necessary, but keep
        // things real-valued for ease.
        DNekMat jac(m_nVars, m_nVars, 0.0, eFULL);
        DNekMat jacDiag(m_nVars, m_nVars, 0.0, eDIAGONAL);
        DNekMat jacInvert(m_nVars, m_nVars, 0.0, eFULL);

        DNekMat metric(m_nVars, m_nVars, 0.0, eFULL);

        for (int k = 0; k < m_nVars; ++k)
        {
            jacDiag(k, k) = (double)rand() / (double)RAND_MAX;

            for (int l = 0; l < m_nVars; ++l)
            {
                jac(l, k) = (double)rand() / (double)RAND_MAX;
            }
        }

        metric = jac * jacDiag;
        jac.Invert();
        metric = metric * jac;

        // Compute eigenvalues/eigenvectors of the metric tensor using
        // ideal mapping.
        char jobvl = 'N', jobvr = 'V';
        int worklen = 8 * m_nVars, info;

        DNekMat eval(m_nVars, m_nVars, 0.0, eDIAGONAL); // Eigen Values
        DNekMat evec(m_nVars, m_nVars, 0.0, eFULL);     // Eigen Vectors
        Array<OneD, NekDouble> vl(m_nVars * m_nVars);
        Array<OneD, NekDouble> work(worklen);
        Array<OneD, NekDouble> wi(m_nVars);

        Lapack::Dgeev(jobvl, jobvr, m_nVars, metric.GetRawPtr(), m_nVars,
                      &(eval.GetPtr())[0], &wi[0], &vl[0], m_nVars,
                      &(evec.GetPtr())[0], m_nVars, &work[0], worklen, info);

        for (int k = 0; k < m_nVars; k++)
        {
            m_freqs[k]  = (double)rand() / (double)RAND_MAX;
            m_phases[k] = (double)rand() / (double)RAND_MAX;

            m_lambda[k] = eval(k, k);

            m_z0[k] = (double)rand() / (double)RAND_MAX;
        }
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void EvaluateSinusoidTerm(
        const Array<OneD, const Array<OneD, double>> &inarray,
        Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const;

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
  OneDFDESolver(int nVars, int nPoints, int nTimeSteps, bool test,
                NekDouble alpha)
        : DemoSolver(nVars, nPoints, nTimeSteps, test)
    {
        m_fileName = std::string("OneDFDESolver");
        m_title    = std::string("Solution to the 1D constant equation");

        m_alpha = alpha;

        if (test)
            srand(0);
        else
            srand(time(NULL));

        // Initial values set to zero
        m_u0 = Array<OneD, Array<OneD, NekDouble>>(m_nVars);

        for (int i = 0; i < m_nVars; i++)
        {
            m_u0[i] = Array<OneD, NekDouble>(m_nPoints, 0.0);

            for (int j = 0; j < m_nPoints; j++)
            {
                m_u0[i][j] = (double)rand() / (double)RAND_MAX;
            }
        }
    }

    // -----------------------------------------------------------------
    // These functions/methods below are the routines which will be
    // used by the TimeIntegration framework (and are required for
    // using it).
    void EvaluateFDETerm(const Array<OneD, const Array<OneD, double>> &inarray,
                         Array<OneD, Array<OneD, double>> &outarray,
                         const NekDouble time) const;

    // -----------------------------------------------------------------
    void EvaluateExactSolution(Array<OneD, Array<OneD, double>> &outarray,
                               const NekDouble time) const;
    // -----------------------------------------------------------------

private:
    // Value of the coefficients:
    NekDouble m_alpha; // Fractional order

    Array<OneD, Array<OneD, NekDouble>> m_u0;

}; // end class OneDFDESolver

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    // So to be case insensitive for the --L2 or --l2 option.
    namespace po_style = boost::program_options::command_line_style;

    po::options_description desc("Usage:");

    desc.add_options()("help,h", "Produce this help message.")(
        "test", "Run in regession test mode.")(
        "butcher,b", "Print the Butcher Tableau for each phase.")(
        "solution,s", "Print the solution values for each time step.")(
        "L2,l", "Print the L2 error for each time step.")(
        "dof,d", po::value<int>(),
        "Number of degrees of freedom (points or values).")(
        "timesteps,t", po::value<int>(), "Number of timesteps.")(
        "order,o", po::value<int>(), "Order of the scheme.")(
        "parameter,p", po::value<std::vector<NekDouble>>()->multitoken(),
        "Free parameters for the scheme.")(
        "variant,v", po::value<std::string>(),
        "Method variant."
        "- Forward:  1st order Forward  Euler\n"
        "- Backward: 1st order Backward Euler\n"
        "  \n"
        "- SSP: 1st-3rd order multi-stage Runga-Kutta SSP scheme\n"
        "  \n"
        "- Gear: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
        "- DIRK: 1st-4th order multi-stage IMEX DIRK scheme\n"
        "  \n"
        "- Lawson:  1st order multi-step Lawson-Euler  exponential scheme\n"
        "- Norsett: 1st-4th order multi-step Norsett-Euler exponential "
        "scheme\n")(
        "method,m", po::value<std::string>(),
        "Name of the time-integration scheme:\n"
        "- Euler -variant Forward:  1st order Forward  Euler\n"
        "- Euler -variant Backward: 1st order Backward Euler\n"
        "- CNAB: 2nd order Crank-Nicolson/Adams-Bashforth (CNAB)\n"
        "- MCNAB: 2nd order Modified Crank-Nicolson/Adams-Bashforth\n"
        "     (MCNAB)\n"
        "  \n"
        "- RungeKutta: 1st-5th order multi-stage Runga-Kutta scheme\n"
        "- RungeKutta -variant SSP: 1st-4th order multi-stage Runga-Kutta SSP "
        "scheme\n"
        "- DIRK:1st-4th order multi-stage Diagonally Implicit\n"
        "      Runga-Kutta scheme (DIRK)\n"
        "- AdamsBashforth: 1st-4th order multi-step Adams-Bashforth scheme\n"
        "- AdamsMoulton: 1st-4th order multi-step Adams-Moulton scheme\n"
        "- BDFImplicit: 1st-4th order multi-step BDFImplicit scheme\n"
        "- IMEX: 1st-4th order multi-step IMEX scheme\n"
        "- IMEX -variant Gear: 2nd order IMEX Gear (Extrapolated Gear/SBDF-2)\n"
        "- IMEX -variant DIRK: 1st-3rd order multi-stage IMEX DIRK scheme\n"
        "  \n"
        "- EulerExponential -variant Lawson: 1st order multi-step Lawson-Euler "
        "exponential scheme\n"
        "- EulerExponential -variant Norsett:1st-4th order multi-step "
        "Norsett-Euler exponential scheme\n"
        "  \n"
        "- FractionalInTime: 1st-4th order multi-step Fractional In Time "
        "scheme\n");

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .style(po_style::unix_style | po_style::case_insensitive)
                      .run(),
                  vm);
        po::notify(vm);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl << desc;
        return 1;
    }

    // These three parameters are manditory.
    if (!vm.count("dof") || !vm.count("timesteps") || !vm.count("method"))
    {
        std::cout << std::endl
                  << "Please specify the number of dof and timesteps "
                  << "along with the order and method.";

        std::cout << std::endl << std::endl << desc;

        return 1;
    }

    int nDoF            = vm["dof"].as<int>();
    int nTimeSteps      = vm["timesteps"].as<int>();
    std::string sMethod = vm["method"].as<std::string>();

    if (nDoF < 2)
    {
        std::cout << "Please specify the number of "
                  << "dof to be greater than 1." << std::endl
                  << std::endl
                  << desc;
        return 1;
    }

    // Maybe needed parameters
    std::string sVariant =
        vm.count("variant") ? vm["variant"].as<std::string>() : "";
    int nOrder = vm.count("order") ? vm["order"].as<int>() : 0;

    // Check the varaibles and parse the free parameters which are in
    // a string.
    std::vector<NekDouble> freeParams;
    if (vm.count("parameter"))
    {
        freeParams = vm["parameter"].as<std::vector<NekDouble>>();
    }

    // Optional parameters
    bool L2      = vm.count("L2");
    bool printS  = vm.count("solution");
    bool printBT = vm.count("butcher");
    bool test    = vm.count("test");

    // Bullet proofing.

    // These methods need the order only
    if (((sMethod == "RungeKutta"     || sMethod == "DIRK" ||
          sMethod == "AdamsBashforth" || sMethod == "AdamsMoulton" ||
          sMethod == "BDFImplicit"    || sMethod == "EulerExponential") &&
         nOrder == 0) ||

        // No variant but the order
        (sMethod == "IMEX" && sVariant == "" && nOrder == 0) ||

        // Needs an order and parameters
        (sMethod == "IMEX" && sVariant == "DIRK" &&
         (nOrder == 0 || freeParams.size() == 0)) ||

        // Needs an order and possibly parameters
        (sMethod == "FractionalInTime" &&
         (nOrder == 0 ||
          (freeParams.size() != 0 &&
           freeParams.size() != 1 &&
           freeParams.size() != 2 &&
           freeParams.size() != 6))) ||

        // Needs a variant
        (sMethod == "Euler" && sVariant == "") ||

        // Needs a variant and order
        (sMethod == "EulerExponential" && sVariant == "" && nOrder == 0) ||

        vm.count("help"))
    {
        std::cout << std::endl
                  << "Please specify the number of dof and timesteps "
                  << "along with the ";

        if (sMethod == "IMEX" && sVariant == "" && nOrder == 0)
        {
            std::cout << "method and order.";
        }
        else if (sMethod == "IMEX" && sVariant == "DIRK" &&
                 (nOrder == 0 || freeParams.size() == 0))
        {
            std::cout << "method, order, and free parameters "
		      << "<implicit stages, explicit stages>.";
        }

        else if(sMethod == "FractionalInTime" &&
                (nOrder == 0 ||
                 (freeParams.size() != 0 &&
                  freeParams.size() != 1 &&
                  freeParams.size() != 2 &&
                  freeParams.size() != 6)))
        {
            std::cout << "method, order, and optionally the free parameters "
                      << "[alpha] | [alpha, base] | [alpha, base, nQuadPts,sigma,mu0,nu].";
        }
        else if (sMethod == "Euler")
        {
            std::cout << "method and variant.";
        }
        else if (sMethod == "EulerExponential")
        {
            std::cout << "method, variant, and order.";
        }
        else
        {
            std::cout << "method and order.";
        }

        std::cout << std::endl << std::endl << desc;

        return vm.count("help") ? 0 : 1;
    }

    // -------------------------------------------------------------------------
    // The implementation below shows how the time-stepping framework can be
    // used for time-integration PDEs

    // 1. THE TEMPORAL DISCRETISATION
    //    Create time integrator requested.
    TimeIntegrationSchemeFactory &factory =
        LibUtilities::GetTimeIntegrationSchemeFactory();

    LibUtilities::TimeIntegrationSchemeSharedPtr tiScheme =
        factory.CreateInstance(sMethod, sVariant, nOrder, freeParams);

    int nVariables;
    int nPoints;

    // 2. THE SPATIAL DISCRETISATION
    //    Create an object of the DemoSolver class.
    //    This class can be thought of as representing the spatial (finite
    //    difference) discretisation.
    LibUtilities::TimeIntegrationSchemeOperators ode;

    std::shared_ptr<DemoSolver> solverSharedPtr;

    if (tiScheme->GetIntegrationSchemeType() == eFractionalInTime)
    {
        nVariables = nDoF;
        nPoints    = 1;

        NekDouble alpha = 0.3; // Fractional order default

        if( freeParams.size() > 0 )
        {
          alpha = freeParams[0];
        }

        OneDFDESolver *tmpSolver =
          new OneDFDESolver(nVariables, nPoints, nTimeSteps, test, alpha);

        ode.DefineOdeRhs(&OneDFDESolver::EvaluateFDETerm, tmpSolver);

        solverSharedPtr = std::shared_ptr<DemoSolver>(tmpSolver);
    }
    else if (tiScheme->GetIntegrationSchemeType() == eExponential)
    {
        nVariables = nDoF;
        nPoints    = 1;

        OneDSinusoidSolver *tmpSolver =
            new OneDSinusoidSolver(nVariables, nPoints, nTimeSteps, test);

        ode.DefineOdeRhs(&OneDSinusoidSolver::EvaluateSinusoidTerm, tmpSolver);

        solverSharedPtr = std::shared_ptr<DemoSolver>(tmpSolver);

        // For exponential integrators, the coefficents for each
        // variable needs to be set.
        ((EulerExponentialTimeIntegrationScheme *)(&(*tiScheme)))
            ->SetExponentialCoefficients(tmpSolver->GetLambda());
    }
    else
    {
        nVariables = 1;
        nPoints    = nDoF;

        OneDFiniteDiffAdvDiffSolver *tmpSolver =
            new OneDFiniteDiffAdvDiffSolver(nVariables, nPoints, nTimeSteps,
                                            test);

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

        solverSharedPtr = std::shared_ptr<DemoSolver>(tmpSolver);
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
    Array<OneD, Array<OneD, double>> exactSol(nVariables);

    for (int k = 0; k < nVariables; ++k)
    {
        approxSol[k] = Array<OneD, double>(nPoints);
        exactSol[k]  = Array<OneD, double>(nPoints);
    }

    // Initialize both solutions using the exact solution.
    double t0 = solverSharedPtr->GetInitialTime();

    solverSharedPtr->EvaluateExactSolution(approxSol, t0);
    solverSharedPtr->EvaluateExactSolution(exactSol, t0);

    // 3.1 Initialize the time-integration scheme.
    double dt = solverSharedPtr->GetDeltaT();

    // The Fractional in Time needs the end time whereas the GLMs need
    // the initial time.
    if (tiScheme->GetIntegrationSchemeType() == eFractionalInTime)
    {
        double t1 = solverSharedPtr->GetFinalTime();
        tiScheme->InitializeScheme(dt, approxSol, t1, ode);
    }
    else
    {
        tiScheme->InitializeScheme(dt, approxSol, t0, ode);
    }

    // 4. Open a file for writing the solution
    std::ofstream outfile, errorfile, L2Normfile;
    outfile.open(solverSharedPtr->GetFileName() + ".dat");
    errorfile.open(solverSharedPtr->GetFileName() + "Error.dat");
    L2Normfile.open(solverSharedPtr->GetFileName() + "L2Norm.dat");

    // Save the scheme name for outputting.
    solverSharedPtr->SetSchemeName(tiScheme->GetName());

    // Time step and time values
    int timeStep = 0;
    double time  = t0;

    LibUtilities::Timer timer;
    NekDouble intTime = 0.0;

    // Write the time step and time
    if (printS || L2)
        std::cout << "Time step: " << timeStep << "  "
                  << "Time: " << time << std::endl;

    solverSharedPtr->GetMinMaxValues(exactSol, approxSol, printS);

    double L2Norm = solverSharedPtr->EvaluateL2Error(exactSol, approxSol, L2);

    // Write the initial conditions, error, and L2 Norm to the output file.
    solverSharedPtr->AppendOutput(outfile, errorfile, L2Normfile, 0, 0,
                                  exactSol, approxSol, L2Norm);

    // 5. Do the time integration.
    while (timeStep < nTimeSteps)
    {
        // Time integration for one time step
        timer.Start();
        approxSol = tiScheme->TimeIntegrate(timeStep, dt, ode);
        timer.Stop();
        intTime += timer.TimePerTest(1);

        // Advance the time for the exact solution.
        time += dt;

        // Calculate the exact solution
        solverSharedPtr->EvaluateExactSolution(exactSol, time);

        // At this point the time step is finished so increment the time step.
        ++timeStep;

        // 6. Calculate the min / max values. If the last argument is
        // true the values will be dumped to screen.

        // Write the time step and time
        if (printS || L2)
        {
            std::cout << "Time step: " << timeStep << "  "
                      << "Time: " << time << std::endl;
        }

        solverSharedPtr->GetMinMaxValues(exactSol, approxSol, printS);

        L2Norm = solverSharedPtr->EvaluateL2Error(exactSol, approxSol, L2);

        // Save the solutions, error, and L2 Norm to the output file
        solverSharedPtr->AppendOutput(outfile, errorfile, L2Normfile, timeStep,
                                      time, exactSol, approxSol, L2Norm);
    }

    // Printing preable so the user knows what was done.
    if (printBT)
    {
        tiScheme->printFull(std::cout);
    }
    else
    {
        std::cout << tiScheme << std::endl;
    }

    std::cout << "Number of time steps performed: " << timeStep << std::endl
              << "Time increment: " << dt << std::endl
              << "Total time: " << time << std::endl
              << "CPU time : " << std::setw(8) << std::left << intTime
              << std::endl;

    // 7. Calculate the error for the last time step and dump to screen
    solverSharedPtr->EvaluateL2Error(exactSol, approxSol, true);

    // 8. Some more writing out the results
    outfile.close();
    errorfile.close();

    solverSharedPtr->GenerateGnuplotScript(tiScheme->GetFullName());

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void DemoSolver::Project(const Array<OneD, const Array<OneD, double>> &inarray,
                         Array<OneD, Array<OneD, double>> &outarray,
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
    const Array<OneD, const Array<OneD, double>> &exact,
    const Array<OneD, const Array<OneD, double>> &approx, bool print)
{
    // Get the min and max value and write the exact solution
    if (print)
    {
        std::cout << "  exact       ";
    }

    if (print && m_nVars > 1 && m_nPoints > 1)
    {
        std::cout << std::endl;
    }

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            // Get the min/max values only from the exact solution
            // becasue the approximate solution can blow up which
            // casues the axis to be infinite.
            if (m_minValue > exact[k][i])
            {
                m_minValue = exact[k][i];
            }

            if (m_maxValue < exact[k][i])
            {
                m_maxValue = exact[k][i];
            }

            if (print)
            {
                std::cout << exact[k][i] << "  ";
            }
        }

        if (print && m_nVars > 1 && m_nPoints > 1)
        {
            std::cout << std::endl;
        }
    }

    if (print)
    {
        std::cout << std::endl;
    }

    // Get the min and max value and write the approximate solution
    if (print)
    {
        std::cout << "  approximate ";
    }

    if (print && m_nVars > 1 && m_nPoints > 1)
    {
        std::cout << std::endl;
    }

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            // Do not get the min/max values from the approximate
            // solution becasue it can blow up which casues the axis
            // to be infinite.

            if (print)
            {
                std::cout << approx[k][i] << "  ";
            }
        }

        if (print && m_nVars > 1 && m_nPoints > 1)
        {
            std::cout << std::endl;
        }
    }

    if (print)
    {
        std::cout << std::endl;
    }
}

// Calculate the Relative Error L2 Norm (as opposed to the absolute L2
// norm).
double DemoSolver::EvaluateL2Error(
    const Array<OneD, const Array<OneD, double>> &exact,
    const Array<OneD, const Array<OneD, double>> &approx, bool print)
{
    // Calcualate the sum of squares for the L2 Norm.
    double a = 0.0;
    double b = 0.0;

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            double diff = approx[k][i] - exact[k][i];

            if (m_minError > diff)
            {
                m_minError = diff;
            }

            if (m_maxError < diff)
            {
                m_maxError = diff;
            }

            a += diff * diff;
            b += exact[k][i] * exact[k][i];
        }
    }

    // Calculate the relative error L2 Norm.
    double norm = sqrt(a / b);

    if (m_minL2Norm > norm)
    {
        m_minL2Norm = norm;
    }

    if (m_maxL2Norm < norm)
    {
        m_maxL2Norm = norm;
    }

    if (print)
    {
        ASSERTL1(b > DBL_EPSILON,
                 "Exact solution is near zero. L2 Norm is invalid");

        std::cout << "L 2 error :" << norm << std::endl;
    }

    return norm;
}

void DemoSolver::AppendOutput(
    std::ofstream &outfile, std::ofstream &errorfile, std::ofstream &L2Normfile,
    int timeStepNumber, const NekDouble time,
    const Array<OneD, const Array<OneD, double>> &exact,
    const Array<OneD, const Array<OneD, double>> &approx,
    const double L2Norm) const
{
    if (timeStepNumber == 0)
    {
        // Save some data provenance and other useful info in output file...
        outfile << "# Data in this file consists of " << m_nTimeSteps
                << " time steps (there will be" << std::endl
                << "# a blank line between each data set for each time step)."
                << std::endl
                << "#" << std::endl
                << "# Delta T: " << m_dt << std::endl
                << "# Method:  " << m_schemeName << std::endl
                << "#" << std::endl
                << "# There are 3 columns with the following headers:"
                << std::endl
                << "#" << std::endl;

        // Save some data provenance and other useful info in output file...
        errorfile << "# Data in this file consists of " << m_nTimeSteps
                  << " time steps (there will be" << std::endl
                  << "# a blank line between each data set for each time step)."
                  << std::endl
                  << "#" << std::endl
                  << "# Delta T: " << m_dt << std::endl
                  << "# Method:  " << m_schemeName << std::endl
                  << "#" << std::endl
                  << "# There are 2 columns with the following headers:"
                  << std::endl
                  << "#" << std::endl;

        // Save some data provenance and other useful info in output file...
        L2Normfile << "# Data in this file consists of " << m_nTimeSteps
                   << " time steps." << std::endl
                   << "#" << std::endl
                   << "# Delta T: " << m_dt << std::endl
                   << "# Method:  " << m_schemeName << std::endl
                   << "#" << std::endl
                   << "# There are 2 columns with the following headers:"
                   << std::endl
                   << "#" << std::endl
                   << "#     Time      |      L2 Norm" << std::endl
                   << "#" << std::endl;

        if (m_nVars > 1)
        {
            outfile << "#   Variable";
            errorfile << "#   Variable";
        }
        else if (m_nPoints > 1)
        {
            outfile << "#   Location";
            errorfile << "#   Location";
        }

        outfile << "  |  Exact Solution  |  Approximate Solution" << std::endl
                << "#" << std::endl
                << std::endl
                << "# Initial condition (at time " << time << "):" << std::endl;

        errorfile << "  |  Error (approx-exact)" << std::endl
                  << "#" << std::endl
                  << std::endl
                  << "# Initial condition (at time " << time
                  << "):" << std::endl;
    }
    else
    {
        outfile << "# Time step: " << timeStepNumber << ", time: " << time
                << std::endl;
        errorfile << "# Time step: " << timeStepNumber << ", time: " << time
                  << std::endl;
    }

    L2Normfile << std::scientific << std::setw(17) << std::setprecision(10)
               << time << "  " << L2Norm << std::endl;

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            if (m_nVars > 1)
            {
                outfile << std::scientific << std::setw(17)
                        << std::setprecision(10) << k << "  " << exact[k][i]
                        << "  " << approx[k][i] << std::endl;

                errorfile << std::scientific << std::setw(17)
                          << std::setprecision(10) << k << "  "
                          << approx[k][i] - exact[k][i] << std::endl;
            }
            else if (m_nPoints > 1)
            {
                outfile << std::scientific << std::setw(17)
                        << std::setprecision(10) << m_x0 + i * m_dx << "  "
                        << exact[k][i] << "  " << approx[k][i] << std::endl;

                errorfile << std::scientific << std::setw(17)
                          << std::setprecision(10) << m_x0 + i * m_dx << "  "
                          << approx[k][i] - exact[k][i] << std::endl;
            }
        }
    }

    // Gnuplot uses two blank lines between each set of data.
    outfile << std::endl << std::endl;
    errorfile << std::endl << std::endl;
}

void DemoSolver::GenerateGnuplotScript(const std::string &method) const
{
    std::vector<std::string> fileSuffix{"", "Error", "L2Norm"};

    for (int j = 0; j < 3; ++j)
    {
        std::ofstream outfile;
        outfile.open(m_fileName + fileSuffix[j] + ".p");
        outfile << "# Gnuplot script file" << std::endl
                << "set   autoscale" << std::endl
                << "unset log" << std::endl
                << "unset label" << std::endl
                << "set xtic auto" << std::endl
                << "set ytic auto" << std::endl;
        // FIXME
        // << "set title 'Finite Difference Solution to the 1D "
        // << "advection-diffusion equation "
        if (j == 0)
        {
            outfile << "set title 'Approximate vs exact solution using method ";

            if (m_nVars > 1)
                outfile << method << "'" << std::endl
                        << "set xlabel 'variable'" << std::endl
                        << "set ylabel 'u'" << std::endl;
            else
                outfile << method << "'" << std::endl
                        << "set xlabel 'x'" << std::endl
                        << "set ylabel 'u'" << std::endl;
        }
        else if (j == 0)
        {
            outfile << "set title 'Error using method ";

            if (m_nVars > 1)
                outfile << method << "'" << std::endl
                        << "set xlabel 'variable'" << std::endl
                        << "set ylabel 'u'" << std::endl;
            else
                outfile << method << "'" << std::endl
                        << "set xlabel 'x'" << std::endl
                        << "set ylabel 'u'" << std::endl;
        }
        else if (j == 1)
        {
            outfile << "set title 'L2Norm using method ";

            outfile << method << "'" << std::endl
                    << "set xlabel 'time'" << std::endl
                    << "set ylabel 'L2Norm'" << std::endl;
        }

        double minValue;
        double maxValue;

        if (j == 0)
        {
            minValue = m_minValue;
            maxValue = m_maxValue;
        }
        else if (j == 1)
        {
            minValue = m_minError;
            maxValue = m_maxError;
        }
        else if (j == 2)
        {
            minValue = m_minL2Norm;
            maxValue = m_maxL2Norm;
        }

        if (j == 0 || j == 1)
        {
            if (m_nVars > 1)
            {
                outfile << "set xr [" << 0 << ":" << m_nVars - 1 << "]"
                        << std::endl
                        << "set yr [" << minValue << ":" << maxValue << "]"
                        << std::endl;
            }
            else if (m_nPoints > 1)
            {
                outfile << "set xr [" << m_x0 << ":" << m_xend << "]"
                        << std::endl
                        << "set yr [" << minValue << ":" << maxValue << "]"
                        << std::endl;
            }
        }
        else
        {
            outfile << "set xr [" << m_t0 << ":" << m_tend << "]" << std::endl
                    << "set yr [" << minValue << ":" << maxValue << "]"
                    << std::endl;

            outfile << "plot '" << m_fileName << fileSuffix[j]
                    << ".dat' using 1:2 index " << 0
                    << " title 'L2 Norm' with linespoints " << std::endl;
        }

        for (int i = 0; i <= m_nTimeSteps; i++)
        {
            double t = m_t0 + (i * m_dt);

            if (j == 0)
            {
                outfile << "plot '" << m_fileName << fileSuffix[j]
                        << ".dat' using 1:2 index " << i
                        << " title 'Exact Solution (t=" << t
                        << ")' with linespoints "
                        << ", '" << m_fileName << fileSuffix[j]
                        << ".dat' using 1:3 index " << i
                        << " title 'Approximate Solution (t=" << t
                        << ")' with linespoints " << std::endl
                        << "pause " << 4.0 / m_nTimeSteps << std::endl;
            }
            else if (j == 1)
            {
                outfile << "plot '" << m_fileName << fileSuffix[j]
                        << ".dat' using 1:2 index " << i
                        << " title 'Error (t=" << t << ")' with linespoints "
                        << std::endl
                        << "pause " << 4.0 / m_nTimeSteps << std::endl;
            }
        }

        // Keep window open until the user clicks
        outfile << "pause mouse any" << std::endl;

        outfile.close();
    }
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
    double b = 1.0 + 2.0 * lambda * m_D /
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
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
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
            outarray[k][0] = -m_V *
                             (inarray[k][1] - inarray[k][m_nPoints - 2]) /
                             (2.0 * m_dx);
            outarray[k][m_nPoints - 1] = outarray[k][0];

            for (int i = 1; i < m_nPoints - 1; i++)
            {
                outarray[k][i] = -m_V *
                                 (inarray[k][i + 1] - inarray[k][i - 1]) /
                                 (2.0 * m_dx);
            }
        }
        else
        {
            // upwind differences
            for (int i = 1; i < m_nPoints; i++)
            {
                outarray[k][i] =
                    -m_V * (inarray[k][i] - inarray[k][i - 1]) / (m_dx);
            }

            outarray[k][0] = outarray[k][m_nPoints - 1];
        }
    }
}

void OneDFiniteDiffAdvDiffSolver::solveTriDiagMatrix(
    int n, double a, double b, const Array<OneD, const double> &inarray,
    Array<OneD, double> &outarray) const
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
            double x  = m_x0 + i * m_dx;
            double wn = 2.0 * M_PI * m_wavenumber;
            outarray[k][i] =
                exp(-m_D * wn * wn * time) * sin(wn * (x - m_V * time));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void OneDSinusoidSolver::EvaluateSinusoidTerm(
    const Array<OneD, const Array<OneD, double>> &inarray,
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    boost::ignore_unused(inarray);

    for (int k = 0; k < m_nVars; k++)
    {
        for (int i = 0; i < m_nPoints; i++)
        {
            outarray[k][i] = sin(m_freqs[k] * time + m_phases[k]);
        }
    }
}

void OneDSinusoidSolver::EvaluateExactSolution(
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    if (time == GetInitialTime())
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

        // Compute right-hand side of the normal equations
        for (int k = 0; k < m_nVars; k++)
        {
            for (int i = 0; i < m_nPoints; i++)
            {
                double v = cos(m_phases[k]);
                double w = sin(m_phases[k]);

                double sinft = sin(m_freqs[k] * time);
                double cosft = cos(m_freqs[k] * time);

                double lambdaFreq2 = m_freqs[k] * m_freqs[k] +
                                     m_lambda[k].real() * m_lambda[k].real();

                outarray[k][i] =
                    // exp(lambda*T) term
                    (std::exp(m_lambda[k].real() * time) *
                     (m_z0[k] + (m_lambda[k].real() * w + m_freqs[k] * v) /
                                    lambdaFreq2)) +

                    // sin(f*T) term
                    (((m_freqs[k] * w - m_lambda[k].real() * v) / lambdaFreq2) *
                     sinft) +

                    // cos(f*T) term
                    (((-m_lambda[k].real() * w - m_freqs[k] * v) /
                      lambdaFreq2) *
                     cosft);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void OneDFDESolver::EvaluateFDETerm(
    const Array<OneD, const Array<OneD, double>> &inarray,
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    boost::ignore_unused(inarray, time);

    for (int i = 0; i < m_nVars; i++)
    {
        NekDouble w = (NekDouble)(m_nVars - i) / (NekDouble)(m_nVars);

        for (int j = 0; j < m_nPoints; j++)
        {
            outarray[i][j] = w * tgamma(m_alpha + 1.0);
        }
    }
}

void OneDFDESolver::EvaluateExactSolution(
    Array<OneD, Array<OneD, double>> &outarray, const NekDouble time) const
{
    // Computes the exact solution at time t to the ODE
    //
    //   y' = time^alpha
    //   y(0) = y0
    //

    // Compute right-hand side of the normal equations
    for (int i = 0; i < m_nVars; i++)
    {
        NekDouble w = (NekDouble)(m_nVars - i) / (NekDouble)(m_nVars);

        for (int j = 0; j < m_nPoints; j++)
        {
            outarray[i][j] = m_u0[i][j] + w * std::pow(time, m_alpha);
        }
    }
}
