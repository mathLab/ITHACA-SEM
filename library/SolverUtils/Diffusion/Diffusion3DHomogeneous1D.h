///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDG3DHomogeneous1D.h
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
// Description: LDG diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSION3DHOMOGENEOUS1D
#define NEKTAR_SOLVERUTILS_DIFFUSION3DHOMOGENEOUS1D

#include <SolverUtils/Advection/Advection3DHomogeneous1D.h>
#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    namespace SolverUtils
    {
        class Diffusion3DHomogeneous1D : public Diffusion
        {
        public:
            static DiffusionSharedPtr create(std::string diffType)
            {
                return DiffusionSharedPtr(
                    new Diffusion3DHomogeneous1D(diffType));
            }
            static std::string type[];

        protected:
            Diffusion3DHomogeneous1D(std::string diffType);

            LibUtilities::TranspositionSharedPtr               m_trans;
            std::string                                        m_diffType;
            SolverUtils::DiffusionSharedPtr                    m_planeDiff;
            NekDouble                                          m_homoLen;
            std::size_t                                        m_numPoints;
            std::size_t                                        m_numPointsPlane;
            std::size_t                                        m_numPlanes;
            std::size_t                                        m_planeCounter;
            Array<OneD, unsigned int>                          m_planes;
            Array<OneD, unsigned int>                          m_planePos;
            Array<OneD, Array<OneD, NekDouble> >               m_homoDerivStore;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_homoDerivPlane;
            Array<OneD, Array<OneD, NekDouble> >               m_inarrayPlane;
            Array<OneD, Array<OneD, NekDouble> >               m_outarrayPlane;
            Array<OneD, MultiRegions::ExpListSharedPtr>        m_fieldsPlane;
            Array<OneD, Array<OneD, NekDouble> >               m_advVelPlane;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);
            virtual void v_Diffuse(
                const std::size_t                                  nConvective,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);
        };
    }
}

#endif
