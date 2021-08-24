///////////////////////////////////////////////////////////////////////////////
//
// File StdLinSysKey.hpp
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
// Description: Headers/definitions of StdLinSysKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDLINSYSKEY_HPP
#define STDLINSYSKEY_HPP


#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdMatrixKey.h>

namespace Nektar
{
    namespace StdRegions
    {
    
        class StdLinSysKey: public StdMatrixKey
        {  
        public:
            StdLinSysKey( MatrixType matrixType, ExpansionType expansionType, 
              StdExpansion &stdExpansion):
        StdMatrixKey(matrixType,expansionType,stdExpansion)
         {
            }

//            StdLinSysKey( MatrixType matrixType, ExpansionType expansionType,
//              StdExpansion &stdExpansion, LibUtilities::PointsType Ntype):
//        StdMatrixKey(matrixType,expansionType,stdExpansion,StdRegions::NNtype)
//         {
//            }
        
            virtual ~StdLinSysKey()
            {
            }
        
        private:

        };

    } // end of namespace
} // end of namespace

#endif //STDLINSYSKEY_HPP
