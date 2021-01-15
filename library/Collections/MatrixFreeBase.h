///////////////////////////////////////////////////////////////////////////////
//
// File: MatrixFreeBase.h
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
// Description: base class definiitin for matrrix free type 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_COLLECTIONS_MATRIXFREEBASE_H
#define NEKTAR_LIBRARY_COLLECTIONS_MATRIXFREEBASE_H

#include <boost/core/ignore_unused.hpp>
#include <StdRegions/StdExpansion.h>

namespace Nektar
{
namespace Collections
{

    class MatrixFreeBase
    {
    public:
        // Default constructor
        MatrixFreeBase(){}
        
    protected:
        /// flag for padding
        bool m_isPadded{false};
        ///  size after padding 
        unsigned int  m_nElmtPad; 

    };
    
    class MatrixFreeOneInOneOut: protected MatrixFreeBase
    {
    public:
        /// Constructor
        MatrixFreeOneInOneOut(const unsigned int nIn, const unsigned int nOut,
                              const unsigned int nCollSize)
        {
            // Padding if needed
            using vec_t = tinysimd::simd<NekDouble>;
            const auto nElmtNoPad = nCollSize;
            m_nElmtPad = nElmtNoPad;

            if (nElmtNoPad % vec_t::width != 0)
            {
                m_isPadded = true;
                m_nElmtPad = nElmtNoPad + vec_t::width -
                    (nElmtNoPad % vec_t::width);
                m_input = Array<OneD, NekDouble>{nIn * m_nElmtPad, 0.0};
                m_output = Array<OneD, NekDouble>{nOut * m_nElmtPad, 0.0};
            }
        }
    protected:
        /// padded input/output vectors
        Array<OneD, NekDouble> m_input, m_output;
    };


    class MatrixFreeMultiInOneOut: protected MatrixFreeBase
    {
    public:
        /// Constructor
        MatrixFreeMultiInOneOut(const unsigned int coordim,
                                const unsigned int nIn,
                                const unsigned int nOut,
                                const unsigned int nCollSize)
        {
            m_coordim = coordim;

            // Padding if needed
            using vec_t = tinysimd::simd<NekDouble>;
            const auto nElmtNoPad = nCollSize;
            m_nElmtPad = nElmtNoPad;

            if (nElmtNoPad % vec_t::width != 0)
            {
                m_isPadded = true;
                m_nElmtPad = nElmtNoPad + vec_t::width -
                    (nElmtNoPad % vec_t::width);

                m_input = Array<OneD, Array<OneD, NekDouble>> (m_coordim);
                m_input[0] = Array<OneD, NekDouble>{nIn * m_nElmtPad, 0.0};
                m_input[1] = Array<OneD, NekDouble>{nIn * m_nElmtPad, 0.0};
                if (m_coordim == 3)
                {
                    m_input[2] = Array<OneD, NekDouble>{nIn * m_nElmtPad, 0.0};
                }
                m_output = Array<OneD, NekDouble>{nOut * m_nElmtPad, 0.0};
            }
        }
    protected:
        /// coordinates dimension
        unsigned short m_coordim;
        /// padded input/output vectors
        Array<OneD, Array<OneD, NekDouble> > m_input;
        Array<OneD, NekDouble> m_output;
    };

    class MatrixFreeOneInMultiOut: protected MatrixFreeBase
    {
    public:
        /// Constructor
        MatrixFreeOneInMultiOut(const unsigned int coordim,
                                const unsigned int nIn,
                                const unsigned int nOut,
                                const unsigned int nCollSize)
        {
            m_coordim = coordim;

            // Padding if needed
            using vec_t = tinysimd::simd<NekDouble>;
            const auto nElmtNoPad = nCollSize;
            m_nElmtPad = nElmtNoPad;

            if (nElmtNoPad % vec_t::width != 0)
            {
                m_isPadded = true;
                m_nElmtPad = nElmtNoPad + vec_t::width -
                    (nElmtNoPad % vec_t::width);

                m_input = Array<OneD, NekDouble>{nIn * m_nElmtPad, 0.0};

                m_output = Array<OneD, Array<OneD, NekDouble>> (m_coordim);
                m_output[0] = Array<OneD, NekDouble>{nOut * m_nElmtPad, 0.0};
                m_output[1] = Array<OneD, NekDouble>{nOut * m_nElmtPad, 0.0};
                if (m_coordim == 3)
                {
                    m_output[2] = Array<OneD, NekDouble>{nOut * m_nElmtPad, 0.0};
                }
            }
        }
    protected:
        /// coordinates dimension
        unsigned short m_coordim;
        /// padded input/output vectors
        Array<OneD, NekDouble> m_input;
        Array<OneD, Array<OneD, NekDouble> > m_output;
    };
}
}
#endif // NEKTAR_LIBRARY_COLLECTIONS_MATRIXFREEBASE_H
