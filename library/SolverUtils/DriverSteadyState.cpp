///////////////////////////////////////////////////////////////////////////////
//
// File DriverSteadyState.cpp
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
// Description: Incompressible Navier Stokes solver
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/DriverSteadyState.h>
#include <SolverUtils/AdvectionSystem.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

string DriverSteadyState::className
        = GetDriverFactory().RegisterCreatorFunction(
                "SteadyState", DriverSteadyState::create);
string DriverSteadyState::driverLookupId
        = LibUtilities::SessionReader::RegisterEnumValue(
                "Driver","SteadyState",0);

/**
 *
 */
DriverSteadyState::DriverSteadyState(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverModifiedArnoldi(pSession, pGraph)
{
}


/**
 *
 */
DriverSteadyState:: ~DriverSteadyState()
{
}


/**
 *
 */
void DriverSteadyState::v_InitObject(ostream &out)
{
    DriverModifiedArnoldi::v_InitObject(out);
}


void DriverSteadyState::v_Execute(ostream &out)

{
    // With a loop over "DoSolve", this Driver implements the
    // "encapsulated" Selective Frequency Damping method(SFD)
    // to find the steady state of a flow above the critical Reynolds
    // number.
    m_equ[m_nequ - 1]->PrintSummary(out);

    m_session->LoadParameter("IO_InfoSteps", m_infosteps, 1000);
    m_session->LoadParameter("IO_CheckSteps", m_checksteps, 100000);
    m_session->LoadParameter("ControlCoeff",m_X, 1);
    m_session->LoadParameter("FilterWidth", m_Delta, 2);

    // To evaluate optimum SFD parameters if growth rate provided in the
    // xml file
    m_session->LoadParameter("GrowthRateEV", GrowthRateEV, 0.0);

    // To evaluate optimum SFD parameters if frequency provided in the xml
    // file
    m_session->LoadParameter("FrequencyEV", FrequencyEV, 0.0);

    // Determine when SFD method is converged
    m_session->LoadParameter("TOL", TOL, 1.0e-08);

    // Used only for the Adaptive SFD method
    m_session->LoadParameter("AdaptiveTOL", AdaptiveTOL, 1.0e-02);

    // Used only for the Adaptive SFD method
    m_session->LoadParameter("AdaptiveTime", AdaptiveTime, 50.0*m_Delta);

    if (m_comm->GetRank() == 0)
    {
        PrintSummarySFD();
    }

    timer.Start();

    // Definition of shared pointer (used only for the Adaptive SFD method)
    AdvectionSystemSharedPtr A = m_equ[0]->as<AdvectionSystem>();

    // Condition necessary to run SFD for the compressible case
    NumVar_SFD = m_equ[m_nequ - 1]->UpdateFields()[0]->GetCoordim(0);
    if (m_session->GetSolverInfo("EqType") == "EulerCFE" ||
        m_session->GetSolverInfo("EqType") == "NavierStokesCFE")
    {
        // Number of variables for the compressible equations
        NumVar_SFD += 2;
    }
    if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        if (m_session->GetSolverInfo("HOMOGENEOUS") == "1D")
        {
            NumVar_SFD += 1;
        }
    }

    // We store the time step
    m_dt = m_equ[m_nequ - 1]->GetTimeStep();

    // Evaluate optimum SFD parameters if dominent EV given by xml file
    if (GrowthRateEV != 0.0 && FrequencyEV != 0.0)
    {
        cout << "Besed on the dominant EV given in the xml file,"
             << "a 1D model is used to evaluate the optumum parameters"
             << "of the SFD method:" << endl;
        complex<NekDouble> EV = polar(exp(GrowthRateEV), FrequencyEV);
        GradientDescentMethod(EV, m_X, m_Delta);
    }


    // We set up the elements of the operator of the encapsulated
    // formulation of the selective frequencive damping method
    SetSFDOperator(m_X, m_Delta);

    // m_steps is set to 1. Then "m_equ[m_nequ - 1]->DoSolve()" will run
    // for only one time step
    m_equ[m_nequ - 1]->SetSteps(1);
    ofstream m_file("ConvergenceHistory.txt", ios::out | ios::trunc);

    Array<OneD, Array<OneD, NekDouble> > q0(NumVar_SFD);
    Array<OneD, Array<OneD, NekDouble> > q1(NumVar_SFD);
    Array<OneD, Array<OneD, NekDouble> > qBar0(NumVar_SFD);
    Array<OneD, Array<OneD, NekDouble> > qBar1(NumVar_SFD);

    for(int i = 0; i < NumVar_SFD; ++i)
    {
        q0[i] = Array<OneD, NekDouble> (m_equ[m_nequ-1]->GetTotPoints(),
                                        0.0); //q0 is initialised
        qBar0[i] = Array<OneD, NekDouble> (m_equ[m_nequ-1]->GetTotPoints(),
                                        0.0);
        m_equ[m_nequ - 1]->CopyFromPhysField(i, qBar0[i]);
    }

    ///Definition of variables used in this algorithm
    m_stepCounter               = 0;
    m_Check                     = 0;
    m_Check_BaseFlow            = 1;
    m_NonConvergingStepsCounter = 0;
    Diff_q_qBar                 = 1.0;
    Diff_q1_q0                  = 1.0;
    cpuTime                     = 0.0;
    elapsed                     = 0.0;
    totalTime                   = 0.0;
    FlowPartiallyConverged      = false;

    while (max(Diff_q_qBar, Diff_q1_q0) > TOL)
    {
        ///Call the Navier-Stokes solver for one time step
        m_equ[m_nequ - 1]->DoSolve();

        for(int i = 0; i < NumVar_SFD; ++i)
        {
            ///Copy the current flow field into q0
            m_equ[m_nequ - 1]->CopyFromPhysField(i, q0[i]);

            ///Apply the linear operator to the outcome of the solver
            ComputeSFD(i, q0, qBar0, q1, qBar1);

            ///Update qBar
            qBar0[i] = qBar1[i];

            ///Copy the output of the SFD method into the current flow field
            m_equ[m_nequ - 1]->CopyToPhysField(i, q1[i]);
        }

        if(m_infosteps && !((m_stepCounter+1)%m_infosteps))
        {
            ConvergenceHistory(qBar1, q0, Diff_q_qBar, Diff_q1_q0);

            ///Loop for the adaptive SFD method
            if (m_EvolutionOperator    == eAdaptiveSFD &&
                FlowPartiallyConverged == false)
            {

                m_NonConvergingStepsCounter++;

                if (Diff_q_qBar < AdaptiveTOL)
                {
                    if (m_comm->GetRank() == 0)
                    {
                    cout << "\n\t The SFD method is converging: we compute "
                         << "stability analysis using the 'partially "
                         << "converged' steady state as base flow:\n" << endl;
                    }

                    m_equ[m_nequ - 1]->Checkpoint_BaseFlow(m_Check_BaseFlow);
                    m_Check_BaseFlow++;

                    A->GetAdvObject()->SetBaseFlow(q0,m_equ[0]->UpdateFields());
                    DriverModifiedArnoldi::v_Execute(out);

                    if (m_comm->GetRank() == 0)
                    {
                        // Compute the update of the SFD parameters only on
                        // one processor
                        ComputeOptimization();
                    }
                    else
                    {
                        // On all the other processors, the parameters are set
                        // to 0
                        m_X = 0;
                        m_Delta = 0;
                    }
                    // The we give to all the processors the value of X and
                    // Delta of the first processor
                    m_comm->AllReduce(m_X,     Nektar::LibUtilities::ReduceSum);
                    m_comm->AllReduce(m_Delta, Nektar::LibUtilities::ReduceSum);

                    FlowPartiallyConverged = true;
                }
                else if (m_NonConvergingStepsCounter * m_dt * m_infosteps
                            >= AdaptiveTime)
                {
                    if (m_comm->GetRank() == 0)
                    {
                        cout << "\n\t We compute stability analysis using"
                             << " the current flow field as base flow:\n"
                             << endl;
                    }

                    m_equ[m_nequ - 1]->Checkpoint_BaseFlow(m_Check_BaseFlow);
                    m_Check_BaseFlow++;

                    A->GetAdvObject()->SetBaseFlow(q0,m_equ[0]->UpdateFields());
                    DriverModifiedArnoldi::v_Execute(out);

                    if (m_comm->GetRank() == 0)
                    {
                        // Compute the update of the SFD parameters only on
                        // one processor
                        ComputeOptimization();
                    }
                    else
                    {
                        // On all the other processors, the parameters are set
                        // to 0
                        m_X = 0;
                        m_Delta = 0;
                    }
                    // The we give to all the processors the value of X and
                    // Delta of the first processor
                    m_comm->AllReduce(m_X, Nektar::LibUtilities::ReduceSum);
                    m_comm->AllReduce(m_Delta, Nektar::LibUtilities::ReduceSum);

                    m_NonConvergingStepsCounter = 0;
                }
            }
        }

        if(m_checksteps && m_stepCounter&&(!((m_stepCounter+1)%m_checksteps)))
        {
            m_Check++;
            m_equ[m_nequ - 1]->Checkpoint_Output(m_Check);
        }
        m_stepCounter++;
    }

    m_file.close();

    ///We save the final solution into a .fld file
    m_equ[m_nequ - 1]->Output();

    for(int j = 0; j < m_equ[m_nequ - 1]->GetNvariables(); ++j)
    {
        NekDouble vL2Error = m_equ[m_nequ - 1]->L2Error(j,false);
        NekDouble vLinfError = m_equ[m_nequ - 1]->LinfError(j);
        if (m_comm->GetRank() == 0)
        {
            out << "L 2 error (variable " << m_equ[m_nequ - 1]->GetVariable(j)
                << ") : " << vL2Error << endl;
            out << "L inf error (variable " << m_equ[m_nequ - 1]->GetVariable(j)
                << ") : " << vLinfError << endl;
        }
    }
}


/**
 * This routine defines the encapsulated SFD operator with first-order
 * splitting and exact resolution of the second subproblem.
 *(See http://scitation.aip.org/content/aip/journal/pof2/26/3/10.1063/1.4867482 for details)
 */
void DriverSteadyState::SetSFDOperator(const NekDouble X_input,
                       const NekDouble Delta_input)
{
    NekDouble X = X_input*m_dt;
    NekDouble Delta = Delta_input/m_dt;
    NekDouble coeff = 1.0/(1.0 + X*Delta);
    M11 = coeff*(1.0 + X*Delta*exp(-(X + 1.0/Delta)));
    M12 = coeff*(X*Delta*(1.0 - exp(-(X + 1.0/Delta))));
    M21 = coeff*(1.0 - exp(-(X + 1.0/Delta)));
    M22 = coeff*(X*Delta + exp(-(X + 1.0/Delta)));
}


void DriverSteadyState::ComputeSFD(const int i,
                   const Array<OneD, const Array<OneD, NekDouble> > &q0,
                const Array<OneD, const Array<OneD, NekDouble> > &qBar0,
                Array<OneD, Array<OneD, NekDouble> > &q1,
                Array<OneD, Array<OneD, NekDouble> > &qBar1)
{
    q1[i]    = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(),0.0);
    qBar1[i] = Array<OneD, NekDouble> (m_equ[m_nequ - 1]->GetTotPoints(),0.0);

    ///Encapsulated SFD method
    Vmath::Svtvp(q1[i].size(), M11, q0[i],    1, q1[i], 1, q1[i], 1 );
    Vmath::Svtvp(q1[i].size(), M12, qBar0[i], 1, q1[i], 1, q1[i], 1 );

    Vmath::Svtvp(qBar1[i].size(), M21, q0[i],    1, qBar1[i], 1,
                                                            qBar1[i], 1 );
    Vmath::Svtvp(qBar1[i].size(), M22, qBar0[i], 1, qBar1[i], 1,
                                                            qBar1[i], 1 );
}


void DriverSteadyState::ComputeOptimization()
{
    NekDouble growthEV(0.0);
    NekDouble frequencyEV(0.0);

    // m_kdim is the dimension of Krylov subspace (defined in the xml file and
    // used in DriverArnoldi.cpp)
    ReadEVfile(m_kdim, growthEV, frequencyEV);

    cout << "\n\tgrowthEV = " << growthEV << endl;
    cout << "\tfrequencyEV = " << frequencyEV << endl;

    complex<NekDouble> ApproxEV = polar(exp(growthEV), frequencyEV);

    NekDouble X_new = m_X;
    NekDouble Delta_new = m_Delta;

    GradientDescentMethod(ApproxEV, X_new, Delta_new);

    m_X = X_new;
    m_Delta = Delta_new;

    SetSFDOperator(m_X, m_Delta);
}

/**
 * This routine implements a gradient descent method to find the parameters X
 * end Delta which give the minimum eigenlavue of the SFD problem applied to
 * the scalar case u(n+1) = \alpha*u(n).
 */
void DriverSteadyState::GradientDescentMethod(
        const complex<NekDouble> &alpha,
                       NekDouble &X_output,
                       NekDouble &Delta_output)
{
    cout << "\n\tWe enter the Gradient Descent Method [...]" << endl;
    bool OptParmFound = false;
    bool Descending = true;
    NekDouble X_input = X_output;
    NekDouble Delta_input = Delta_output;

    NekDouble X_init = X_output;
    NekDouble Delta_init = Delta_output;
    int stepCounter(0);

    NekDouble F0(0.0);
    NekDouble F0x(0.0);
    NekDouble F0y(0.0);
    NekDouble F1(0.0);
    NekDouble dx = 0.00000001;
    NekDouble dirx(0.0);
    NekDouble diry(0.0);
    NekDouble s(0.0);
    NekDouble CurrentMin = 1.0;

    while (OptParmFound == false)
    {
        Descending = true;
        EvalEV_ScalarSFD(X_input, Delta_input, alpha, F0);
        EvalEV_ScalarSFD(X_input + dx, Delta_input, alpha, F0x);
        EvalEV_ScalarSFD(X_input, Delta_input + dx, alpha, F0y);
        dirx =  (F0 - F0x);
        diry =  (F0 - F0y);
        s = abs(0.000001/dirx);
        X_output = X_input + s*dirx;
        Delta_output = Delta_input + s*diry;
        F1 = F0;

        while (Descending == true)
        {
            CurrentMin = F1;
            X_input = X_output;
            Delta_input = Delta_output;
            EvalEV_ScalarSFD(X_output, Delta_output, alpha, F1);

            if (F1 > CurrentMin)
            {
                Descending = false;
            }
            else
            {
                s = s + s*0.01;
                X_output = X_input + s*dirx;
                Delta_output = Delta_input + s*diry;
            }

            if(stepCounter > 9999999)
            {
                //We are stuck in this loop..
                //Then we restart it with different initail conditions
                Descending = false;
                X_input = X_init;
                Delta_init = Delta_init + Delta_init*0.1;
                Delta_input = Delta_init;
                stepCounter = 0;
            }
            stepCounter++;
        }

        if (abs(F0-F1) < dx)
        {
            cout << "\tThe Gradient Descent Method has converged!" << endl;
            EvalEV_ScalarSFD(X_output, Delta_output, alpha, F1);
            cout << "\n\tThe updated parameters are: X_tilde = " << X_output
                 << " and Delta_tilde = " << Delta_output << endl;
            OptParmFound = true;
        }

    }
}


/**
 * This routine evaluates the maximum eigenvalue of the SFD system when applied
 * to the 1D model u(n+1) = alpha*u(n)
 */
void DriverSteadyState::EvalEV_ScalarSFD(
        const NekDouble &X_input,
        const NekDouble &Delta_input,
        const complex<NekDouble> &alpha,
              NekDouble &MaxEV)
{
    NekDouble A11 = ( 1.0 + X_input * Delta_input *
            exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);
    NekDouble A12 = ( X_input*Delta_input - X_input * Delta_input *
            exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);
    NekDouble A21 = ( 1.0 - 1.0 *
            exp(-(X_input + 1.0/Delta_input)) )/(1 + X_input*Delta_input);
    NekDouble A22 = ( X_input*Delta_input + 1.0 *
            exp(-(X_input + 1.0/Delta_input)) )/(1.0 + X_input*Delta_input);

    complex<NekDouble> B11 = alpha;
    NekDouble B12 = 0.0;
    NekDouble B21 = 0.0;
    NekDouble B22 = 1.0;

    complex<NekDouble> a = A11*B11 + A12*B21;
    NekDouble b = A11*B12 + A12*B22;
    complex<NekDouble> c = A21*B11 + A22*B21;
    NekDouble d = A21*B12 + A22*B22;

    complex<NekDouble> delt = sqrt((a-d)*(a-d) + 4.0*b*c);
    complex<NekDouble> lambda_1 = 0.5*(a+d + delt);
    complex<NekDouble> lambda_2 = 0.5*(a+d - delt);

    NekDouble NORM_1 = abs(lambda_1);
    NekDouble NORM_2 = abs(lambda_2);

    MaxEV = max(NORM_1, NORM_2);
}


void DriverSteadyState::ReadEVfile(
        int &KrylovSubspaceDim,
        NekDouble &growthEV,
        NekDouble &frequencyEV)
{
    // This routine reads the .evl file written by the Arnoldi algorithm
    // (written in September 2014)
    std::string   line;
    int           NumLinesInFile = 0;
    std::string   EVfileName = m_session->GetSessionName() +  ".evl";
    std::ifstream EVfile(EVfileName.c_str());
    ASSERTL0(EVfile.good(), "Cannot open .evl file.");

    if(EVfile)
    {
        // This block counts the total number of lines of the .evl file
        // We keep going util we reach the end of the file
        while(getline(EVfile, line))
        {
            ++NumLinesInFile;
        }

        // It may happen that the Stability method that have produced the .elv
        // file converges in less than m_kdim iterations. In this case,
        // KrylovSubspaceDim has to be changed here
        if(NumLinesInFile < KrylovSubspaceDim*2.0
                                + KrylovSubspaceDim*(KrylovSubspaceDim+1.0)/2.0)
        {
            for(int i = 1; i <= KrylovSubspaceDim; ++i)
            {
                if(NumLinesInFile == i*2.0 + i*(i+1.0)/2.0)
                {
                    KrylovSubspaceDim = i;
                }
            }
        }

        // go back to the beginning of the file
        EVfile.clear();
        EVfile.seekg(0, ios::beg);

        // We now want to go to the line where the most unstable eigenlavue was
        // written
        for(int i = 0; i < (NumLinesInFile - KrylovSubspaceDim); ++i)
        {
            std::getline(EVfile, line);
            cout << "Discard line: " << line << endl;
        }

        std::vector<std::string> tokens;
        std::getline(EVfile, line);
        boost::algorithm::split(tokens, line,
                boost::is_any_of("\t "),boost::token_compress_on);

        ASSERTL0(tokens.size() >= 5,
                 "Unexpected formatting of .evl file while reading line:\n"
                 + line);
        growthEV = boost::lexical_cast<NekDouble>(tokens[4]);
        frequencyEV = boost::lexical_cast<NekDouble>(tokens[5]);
    }
    else
    {
        cout << "An error occurred when opening the .evl file" << endl;
    }
    EVfile.close();
}


/**
 * This routine evaluates |q-qBar|_inf (and |q1-q0|_inf) and writes the values
 * in "ConvergenceHistory.txt"
 */
void DriverSteadyState::ConvergenceHistory(
        const Array<OneD, const Array<OneD, NekDouble> > &qBar1,
        const Array<OneD, const Array<OneD, NekDouble> > &q0,
              NekDouble &MaxNormDiff_q_qBar,
              NekDouble &MaxNormDiff_q1_q0)
{
    Array<OneD, NekDouble > NormDiff_q_qBar(NumVar_SFD, 1.0);
    Array<OneD, NekDouble > NormDiff_q1_q0(NumVar_SFD, 1.0);
    MaxNormDiff_q_qBar=0.0;
    MaxNormDiff_q1_q0=0.0;

    for(int i = 0; i < NumVar_SFD; ++i)
    {
        NormDiff_q_qBar[i] = m_equ[m_nequ - 1]->LinfError(i, qBar1[i]);
        NormDiff_q1_q0[i] = m_equ[m_nequ - 1]->LinfError(i, q0[i]);

        if (MaxNormDiff_q_qBar < NormDiff_q_qBar[i])
        {
            MaxNormDiff_q_qBar = NormDiff_q_qBar[i];
        }
        if (MaxNormDiff_q1_q0 < NormDiff_q1_q0[i])
        {
            MaxNormDiff_q1_q0 = NormDiff_q1_q0[i];
        }
    }

    timer.Stop();
    elapsed  = timer.TimePerTest(1);
    cpuTime += elapsed;
    totalTime += elapsed;

    if (m_comm->GetRank() == 0)
    {
        cout << "SFD - Step: " <<  left <<  m_stepCounter+1
             << ";\tTime: " << left << m_equ[m_nequ - 1]->GetFinalTime()
             << ";\tCPU time = " << left << cpuTime << " s"
             << ";\tTot time = " << left << totalTime << " s"
             << ";\tX = " << left << m_X
             << ";\tDelta = " << left << m_Delta
             << ";\t|q-qBar|inf = " << left << MaxNormDiff_q_qBar << endl;
        std::ofstream m_file( "ConvergenceHistory.txt", std::ios_base::app);
        m_file << m_stepCounter+1 << "\t"
               << m_equ[m_nequ - 1]->GetFinalTime() << "\t"
               << totalTime << "\t"
               << MaxNormDiff_q_qBar << "\t"
               << MaxNormDiff_q1_q0 << endl;
        m_file.close();
    }

    cpuTime = 0.0;
    timer.Start();
}


void DriverSteadyState::PrintSummarySFD()
{
    cout << "\n====================================="
            "=================================="
         << endl;
    cout << "Parameters for the SFD method:" << endl;
    cout << "\tControl Coefficient: X = " << m_X << endl;
    cout << "\tFilter Width: Delta = " << m_Delta << endl;
    cout << "The simulation is stopped when |q-qBar|inf < " << TOL << endl;
    if (m_EvolutionOperator == eAdaptiveSFD)
    {
        cout << "\nWe use the adaptive SFD method:" << endl;
        cout << "  The parameters are updated every " << AdaptiveTime
             << " time units;" << endl;
        cout << "  until |q-qBar|inf becomes smaller than " << AdaptiveTOL
             << endl;
    }
    cout << "====================================="
            "==================================\n" << endl;
}

}
}

