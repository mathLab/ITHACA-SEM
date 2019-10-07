///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy1D.cpp
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
// Description: Outputs orthogonal expansion of 1D elements.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <SolverUtils/Filters/FilterEnergy1D.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterEnergy1D::className = GetFilterFactory().
    RegisterCreatorFunction("Energy1D", FilterEnergy1D::create);

/**
 * @brief Set up filter with output file and frequency parameters.
 *
 * @param pSession  Current session.
 * @param pParams   Map of parameters defined in XML file.
 */
FilterEnergy1D::FilterEnergy1D(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>      &pEquation,
    const ParamMap &pParams) :
    Filter(pSession, pEquation),
    m_index(0)
{
    ASSERTL0(pSession->GetComm()->GetSize() == 1,
             "The 1D energy filter currently only works in serial.");

    std::string outName;

    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        outName = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        outName = it->second;
    }
    outName += ".eny";
    m_out.open(outName.c_str());

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
}

/**
 * @brief Destructor.
 */
FilterEnergy1D::~FilterEnergy1D()
{

}

/**
 * @brief Initialize filter.
 */
void FilterEnergy1D::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(time);
    ASSERTL0(pFields[0]->GetExp(0)->GetNumBases() == 1,
             "The Energy 1D filter is only valid in 1D.");
}

/**
 * @brief Update filter output with the current timestep's orthogonal
 * coefficients.
 */
void FilterEnergy1D::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    int nElmt = pFields[0]->GetExpSize();

    // Loop over all elements
    m_out << "##" << endl;
    m_out << "## Time = " << time << endl;
    m_out << "##" << endl;

    for (int i = 0; i < nElmt; ++i)
    {
        // Figure out number of modes in this expansion.
        LocalRegions::ExpansionSharedPtr exp = pFields[0]->GetExp(i);
        int nModes = exp->GetBasis(0)->GetNumModes();

        // Set uo basis key for orthogonal basis
        LibUtilities::BasisType btype = LibUtilities::eOrtho_A;
        LibUtilities::BasisKey bkeyOrth(
            btype, nModes, exp->GetBasis(0)->GetPointsKey());

        // Get basis key for existing expansion
        LibUtilities::BasisKey bkey(
            exp->GetBasis(0)->GetBasisType(),
            exp->GetBasis(0)->GetNumModes(),
            exp->GetBasis(0)->GetPointsKey());

        // Find coeffs for this element in the list of all coefficients
        Array<OneD, NekDouble> coeffs =
            pFields[0]->GetCoeffs() + pFields[0]->GetCoeff_Offset(i);

        // Storage for orthogonal coefficients
        Array<OneD, NekDouble> coeffsOrth(nModes);

        // Project from coeffs -> orthogonal coeffs
        LibUtilities::InterpCoeff1D(bkey, coeffs, bkeyOrth, coeffsOrth);

        // Write coeffs to file
        m_out << "# Element " << i << " (ID "
              << exp->GetGeom()->GetGlobalID() << ")" << endl;
        for (int j = 0; j < nModes; ++j)
        {
            m_out << coeffsOrth[j] << endl;
        }
    }
    m_out << endl;
}

void FilterEnergy1D::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);
    m_out.close();
}

bool FilterEnergy1D::v_IsTimeDependent()
{
    return true;
}
}
}
