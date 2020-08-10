///////////////////////////////////////////////////////////////////////////////
//
// File FilterAeroForcesSPM.cpp
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
// Description: Output values of aerodynamic forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Filters/FilterAeroForcesSPM.h>

using namespace std;

namespace Nektar
{
std::string FilterAeroForcesSPM::className =
        SolverUtils::GetFilterFactory().RegisterCreatorFunction(
                "AeroForcesSPM", FilterAeroForcesSPM::create);

/**
 *
 */
FilterAeroForcesSPM::FilterAeroForcesSPM(
    const LibUtilities::SessionReaderSharedPtr       &pSession,
    const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    if (!(m_outputFile.length() >= 4
          && m_outputFile.substr(m_outputFile.length() - 4) == ".fce"))
    {
        m_outputFile += ".fce";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // Time after which we need to calculate the forces
    it = pParams.find("StartTime");
    if (it == pParams.end())
    {
        m_startTime = 0;
    }
    else
    {
        LibUtilities::Equation equ(
            m_session->GetInterpreter(), it->second);
        m_startTime = equ.Evaluate();
    }
}


/**
 *
 */
FilterAeroForcesSPM::~FilterAeroForcesSPM()
{

}

/**
 *
 */
void FilterAeroForcesSPM::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Save space dimension
    m_spaceDim = pFields[0]->GetGraph()->GetMeshDimension();

    // Fill in the directions vector
    m_dirNames.push_back("X");
    if (m_spaceDim > 1)
    {
        m_dirNames.push_back("Y");
    }
    if (m_spaceDim > 2)
    {
        m_dirNames.push_back("Z");
    }

    // Allocate the aero-forces vector
    m_Forces = Array<OneD, NekDouble>(m_spaceDim);

    // Write header
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    if (vComm->GetRank() == 0)
    {
        // Open output stream
        bool adaptive;
        m_session->MatchSolverInfo("Driver", "Adaptive",
                                    adaptive, false);
        if (adaptive)
        {
            m_outputStream.open(m_outputFile.c_str(), ofstream::app);
        }
        else
        {
            m_outputStream.open(m_outputFile.c_str());
        }
        m_outputStream << "# Forces acting on bodies\n";
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        for (int i = 0; i < m_spaceDim; ++i)
        {
            m_outputStream.width(14);
            m_outputStream <<  "F_" << m_dirNames[i];
        }

        m_outputStream << endl;
    }

    m_index = 0;
}

/**
 *
 */
void FilterAeroForcesSPM::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency  || (time < m_startTime))
    {
        return;
    }

    // Communicators to exchange results
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    //Write Results
    if (vComm->GetRank() == 0)
    {
        // Write time
        m_outputStream.width(8);
        m_outputStream << setprecision(6) << time;
        // Write forces
        for (int i = 0; i < m_spaceDim; ++i)
        {
            m_outputStream.width(15);
            m_outputStream << setprecision(8)
                           << m_Forces[i];
        }
        m_outputStream.width(10);
        m_outputStream << endl;
    }
}


/**
 *
 */
void FilterAeroForcesSPM::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}


/**
 *
 */
bool FilterAeroForcesSPM::v_IsTimeDependent()
{
    return true;
}

/**
 * @brief Determine the total force on the body defined by \f$\Phi\f$
 * (note that if the shape function represents more than one
 * body, this function calculates the value of the final force after adding
 * up the values for each body). This value must be scaled with the
 * density to get the real force vector.
 *
 * @param pIntVel
 * @param pUpPrev
 * @param pPhi
 * @param time
 * @param dt
 */
void FilterAeroForcesSPM::CalculateForces(
        const Array<OneD, Array<OneD, NekDouble> > &pIntVel,
        const Array<OneD, Array<OneD, NekDouble> > &pUpPrev,
        const MultiRegions::ExpListSharedPtr &pPhi,
        NekDouble time,
        NekDouble dt)
{
    // Only output every m_outputFrequency.
    if ((m_index % m_outputFrequency)  || (time < m_startTime))
    {
        return;
    }

    int nq = pIntVel[0].size();
    Array<OneD, NekDouble> tmp(nq);

    // Aerodynamic forces are computed according eq. 18a in Luo et al. (2009).
    // Smoothed profile method for particulate flows: Error analysis and
    // simulations. Journal of Computational Physics, 228(5)
    for (int i = 0; i < m_spaceDim; ++i)
    {
        // "Scalar" field to be integrated
        Vmath::Vsub(nq, pIntVel[i], 1, pUpPrev[i], 1,
                    tmp, 1);
        Vmath::Vmul(nq, pPhi->GetPhys(), 1, tmp, 1, tmp, 1);

        // Integration of force throughout the domain
        m_Forces[i] = pPhi->Integral(tmp) / dt;
    }
}

}
