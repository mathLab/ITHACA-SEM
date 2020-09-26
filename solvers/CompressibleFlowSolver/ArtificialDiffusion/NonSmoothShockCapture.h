///////////////////////////////////////////////////////////////////////////////
//
// File: NonSmoothShockCapture.h
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
// Description: NonSmooth artificial diffusion for shock capture
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_NONSMOOTH
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_ARTIFICIALDIFFUSION_NONSMOOTH

#include "ArtificialDiffusion.h"


namespace Nektar
{

/**
 * @brief Non Smooth artificial diffusion for shock capture for compressible
 * flow problems.
*/
class NonSmoothShockCapture : public ArtificialDiffusion
{
    public:

        friend class MemoryManager<NonSmoothShockCapture>;

        /// Creates an instance of this class
        static ArtificialDiffusionSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const int spacedim)
        {
            ArtificialDiffusionSharedPtr p =
                                MemoryManager<NonSmoothShockCapture>::
                                AllocateSharedPtr(pSession, pFields, spacedim);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:

        virtual void v_GetArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble  >             &mu);

    private:
        NonSmoothShockCapture(
               const LibUtilities::SessionReaderSharedPtr& pSession,
               const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
               const int spacedim);

        virtual ~NonSmoothShockCapture(void){};

        /// Parameters
        int             m_offset;
};

}

#endif
