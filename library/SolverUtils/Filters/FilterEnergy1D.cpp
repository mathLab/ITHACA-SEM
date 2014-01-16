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
// Description: Outputs solution fields during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <SolverUtils/Filters/FilterEnergy1D.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterEnergy1D::className = GetFilterFactory().
            RegisterCreatorFunction("Energy1D", FilterEnergy1D::create);

        FilterEnergy1D::FilterEnergy1D(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            m_out.open("energy-out.txt");
        }

        FilterEnergy1D::~FilterEnergy1D()
        {

        }

        void FilterEnergy1D::v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
        }

        void FilterEnergy1D::v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            int nElmt = pFields[0]->GetExpSize();

            // Loop over all elements
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
                Array<OneD, NekDouble> coeffs = pFields[0]->GetCoeffs() + pFields[0]->GetCoeff_Offset(i);

                // Storage for orthogonal coefficients
                Array<OneD, NekDouble> coeffsOrth(nModes);

                // Project from coeffs -> orthogonal coeffs
                LibUtilities::InterpCoeff1D(bkey, coeffs, bkeyOrth, coeffsOrth);

                // Write coeffs to file
                for (int j = 0; j < nModes; ++j)
                {
                    m_out << coeffsOrth[j] << endl;
                }
            }
        }

        void FilterEnergy1D::v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            m_out.close();
        }

        bool FilterEnergy1D::v_IsTimeDependent()
        {
            return true;
        }
    }
}
