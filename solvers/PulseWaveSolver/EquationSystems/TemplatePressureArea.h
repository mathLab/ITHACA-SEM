///////////////////////////////////////////////////////////////////////////////
//
// File LymphaticPressureArea.h
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
// Description: TemplatePressureArea header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_TEMPLATEPRESSUREAREA_H
#define NEKTAR_TEMPLATEPRESSUREAREA_H

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>
#include <string>

namespace Nektar
{

class TemplatePressureArea;

typedef std::shared_ptr<TemplatePressureArea> TemplatePressureAreaSharedPtr;

// A global linear system.
class TemplatePressureArea : public PulseWavePressureArea
{
    public:
        static PulseWavePressureAreaSharedPtr
        create(Array<OneD, MultiRegions::ExpListSharedPtr> &pVessel,
               const LibUtilities::SessionReaderSharedPtr &pSession)
        {
            return MemoryManager<BetaPressureArea>::AllocateSharedPtr(pVessel,
                                                                      pSession);
        }

        static std::string className;

        TemplatePressureArea(Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
                           const LibUtilities::SessionReaderSharedPtr pSession);

        virtual ~TemplatePressureArea();

    protected:
        virtual void v_GetPressure(NekDouble &P, const NekDouble &beta,
                const NekDouble &A, const NekDouble &A0, const NekDouble &dAUdx,
                      const NekDouble &gamma = 0, const NekDouble &alpha = 0.5);

        virtual void v_GetC(NekDouble &c, const NekDouble &beta,
         const NekDouble &A, const NekDouble &A0, const NekDouble &alpha = 0.5);

        virtual void v_GetW1(NekDouble &W1, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        virtual void v_GetW2(NekDouble &W2, const NekDouble &u,
                 const NekDouble &beta, const NekDouble &A, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        virtual void v_GetAFromChars(NekDouble &A, const NekDouble &W1,
                const NekDouble &W2, const NekDouble &beta, const NekDouble &A0,
                                                  const NekDouble &alpha = 0.5);

        virtual void v_GetUFromChars(NekDouble &u, const NekDouble &W1,
                                                           const NekDouble &W2);

        virtual void v_GetCharIntegral(NekDouble &I, const NekDouble &beta,
         const NekDouble &A, const NekDouble &A0, const NekDouble &alpha = 0.5);

        virtual void v_GetJacobianInverse(NekMatrix<NekDouble> &invJ,
             const Array<OneD, NekDouble> &Au, const Array<OneD, NekDouble> &uu,
           const Array<OneD, NekDouble> &beta, const Array<OneD, NekDouble> &A0,
                  const Array<OneD, NekDouble> &alpha, const std::string &type);

    private:
};

} // namespace Nektar

#endif
