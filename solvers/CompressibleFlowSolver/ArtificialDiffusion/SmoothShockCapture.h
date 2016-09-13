///////////////////////////////////////////////////////////////////////////////
//
// File: SmoothShockCapture.h
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
// Description: Smooth artificial diffusion for shock capture
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_SMOOTH
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_SMOOTH

#include "ArtificialDiffusion.h"


namespace Nektar
{

/**
 * @brief Smooth artificial diffusion for shock capture for compressible flow
 * problems.
*/
class SmoothShockCapture : public ArtificialDiffusion
{
    public:

        friend class MemoryManager<SmoothShockCapture>;

        /// Creates an instance of this class
        static ArtificialDiffusionSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim)
        {
            ArtificialDiffusionSharedPtr p = MemoryManager<SmoothShockCapture>::
                                AllocateSharedPtr(pSession, pFields, spacedim);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:

        virtual void v_DoArtificialDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD,       Array<OneD, NekDouble> > &outarray);

        virtual void v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu);

    private:
        SmoothShockCapture(const LibUtilities::SessionReaderSharedPtr& pSession,
               const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
               const int spacedim);

        virtual ~SmoothShockCapture(void){};

        void GetForcingTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > outarrayForcing);
};

}

#endif
