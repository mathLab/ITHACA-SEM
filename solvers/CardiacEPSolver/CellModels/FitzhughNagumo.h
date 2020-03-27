///////////////////////////////////////////////////////////////////////////////
//
// File FitzhughNagumo.h
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
// Description: Fitzhugh-Nagumo phenomological cell model.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_FITZHUGHNAGUMO
#define NEKTAR_SOLVERS_ADRSOLVER_CELLMODELS_FITZHUGHNAGUMO

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    /// FitzHugh-Nagumo model.
    class CellModelFitzHughNagumo : public CellModel
    {
    public:
        /// Creates an instance of this class
        static CellModelSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField)
        {
            return MemoryManager<CellModelFitzHughNagumo>::AllocateSharedPtr(pSession, pField);
        }

        /// Name of class
        static std::string className;

        CellModelFitzHughNagumo(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const MultiRegions::ExpListSharedPtr& pField);

        virtual ~CellModelFitzHughNagumo() {}

    protected:
        virtual void v_Update(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

        virtual void v_GenerateSummary(SummaryList& s);

        virtual void v_SetInitialConditions();

    private:
        NekDouble              m_beta;
        NekDouble              m_epsilon;

        /// Temporary space for storing \f$u^3\f$ when computing reaction term.
        Array<OneD, NekDouble> m_uuu;
    };
}
#endif /* FITZHUGHNAGUMO_H_ */
