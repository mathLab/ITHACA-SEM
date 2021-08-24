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
// Description: Wrapper for Riemann solver scalar.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace SolverUtils
    {
        /**
         * @brief Wrapper class for Riemann solver scalars.
         */
        class HomoRSScalar
        {
        public:
            HomoRSScalar(RSScalarFuncType func,
                         int              nPlanes)
                : m_func       (func),
                  m_planeNumber(0),
                  m_numPlanes  (nPlanes),
                  m_tmp        ()
            {
            }

            const Array<OneD, const NekDouble>& Exec()
            {
                if (m_planeNumber == 0)
                {
                    m_tmp = m_func();
                }

                const int nPts   = m_tmp.size() / m_numPlanes;
                const int offset = m_planeNumber * nPts;

                m_tmp2 = Array<OneD, NekDouble>(nPts, m_tmp + offset);
                m_planeNumber = (m_planeNumber + 1) % m_numPlanes;

                return m_tmp2;
            }

        private:
            RSScalarFuncType             m_func;
            int                          m_planeNumber;
            int                          m_numPlanes;
            Array<OneD, const NekDouble> m_tmp;
            Array<OneD, const NekDouble> m_tmp2;
        };

        /**
         * @brief Wrapper class for Riemann solver scalars.
         */
        class HomoRSVector
        {
        public:
            HomoRSVector(RSVecFuncType func,
                         int           nPlanes,
                         std::string   desc = "")
                : m_func       (func),
                  m_planeNumber(0),
                  m_numPlanes  (nPlanes),
                  m_tmp        (),
                  m_desc       (desc)
            {
            }

            const Array<OneD, const Array<OneD, NekDouble> >& Exec()
            {
                if (m_planeNumber == 0)
                {
                    m_tmp = m_func();
                }

                const int nDim   = m_tmp.size();
                const int nPts   = m_tmp[0].size() / m_numPlanes;
                const int offset = m_planeNumber * nPts;
                m_tmp2 = Array<OneD, Array<OneD, NekDouble> >(nDim);

                for (int i = 0; i < m_tmp.size(); ++i)
                {
                    m_tmp2[i] = Array<OneD, NekDouble>(nPts, m_tmp[i] + offset);
                }

                m_planeNumber = (m_planeNumber + 1) % m_numPlanes;

                return m_tmp2;
            }

        private:
            RSVecFuncType                        m_func;
            int                                  m_planeNumber;
            int                                  m_numPlanes;
            Array<OneD, Array<OneD, NekDouble> > m_tmp;
            Array<OneD, Array<OneD, NekDouble> > m_tmp2;
            std::string                          m_desc;
        };
    }
}
