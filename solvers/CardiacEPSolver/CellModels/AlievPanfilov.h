///////////////////////////////////////////////////////////////////////////////
//
// File AlievPanfilov.h
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
// Description: Aliev-Panfilov phenomological cell model.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_ALIEVPANFILOV
#define NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_ALIEVPANFILOV

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    /// Aliev Panfilov model.
    class CellModelAlievPanfilov : public CellModel
    {
    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<CellModelAlievPanfilov>::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        CellModelAlievPanfilov(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField);

        virtual ~CellModelAlievPanfilov() {}

    protected:
        virtual void v_Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        virtual void v_GenerateSummary(SummaryList& s);

        virtual void v_SetInitialConditions();

    private:
        /// Trigger parameter a.
        NekDouble m_a;
        /// Scaling parameter k.
        NekDouble m_k;
        /// Restitution parameter \f$\mu_1\f$.
        NekDouble m_mu1;
        /// Restitution parameter \f$\mu_2\f$.
        NekDouble m_mu2;
        /// Restitution parameter \f$\epsilon\f$.
        NekDouble m_eps;

        /// Temporary space for storing \f$u^2\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uu;
        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp1;
        /// Workspace for computing reaction term.
        Array<OneD, NekDouble> m_tmp2;
    };
}
#endif /* ALIEVPANFILOV_H_ */
