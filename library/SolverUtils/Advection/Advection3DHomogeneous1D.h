///////////////////////////////////////////////////////////////////////////////
//
// File: Advection3DHomogeneous1D.h
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
// Description: FR 3DHomogeneous1D advection 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTION3DHOMOGENEOUS1D
#define NEKTAR_SOLVERUTILS_ADVECTION3DHOMOGENEOUS1D

#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/Advection/HomogeneousRSScalar.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        class Advection3DHomogeneous1D : public Advection
        {
        public:
            static AdvectionSharedPtr create(std::string advType)
            {
                return AdvectionSharedPtr(
                    new Advection3DHomogeneous1D(advType));
            }
            static std::string type[];

        protected:
            Advection3DHomogeneous1D(std::string advType);

            std::string                                        m_advType;
            SolverUtils::AdvectionSharedPtr                    m_planeAdv;
            int                                                m_numPoints;
            int                                                m_numPointsPlane;
            int                                                m_numPlanes;
            int                                                m_planeCounter;
            Array<OneD, unsigned int>                          m_planes;
            Array<OneD, unsigned int>                          m_planePos;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fluxVecStore;
            Array<OneD, Array<OneD, NekDouble> >               m_inarrayPlane;
            Array<OneD, Array<OneD, NekDouble> >               m_outarrayPlane;
            Array<OneD, MultiRegions::ExpListSharedPtr>        m_fieldsPlane;
            Array<OneD, Array<OneD, NekDouble> >               m_advVelPlane;
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                               m_fluxVecPlane;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr               pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>        pFields);

            virtual void v_Advect(
                const int                                          nConvField,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const NekDouble                                   &time,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

        private:
            void ModifiedFluxVector(
                const Array<OneD, Array<OneD, NekDouble> >         &physfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        };
    }
}

#endif
