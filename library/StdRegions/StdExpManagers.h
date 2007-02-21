///////////////////////////////////////////////////////////////////////////////
//
// File StdExpManagers.h
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
// Description: Headers for StdExpansion Managers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDEXPMANAGERS_H
#define STDEXPMANAGERS_H


#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace StdRegions
    {
	
        class StdExpansion;

        class StdMatrixKey
        {
        public:
            StdMatrixKey( MatrixType matrixType, StdExpansion &stdExpansion);
	    
            virtual ~StdMatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                bool operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            };

            /// Used for finding value given the key in NekManager.
            friend bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            friend bool opLess::operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs);

	    MatrixType GetMatrixType() const
	    {
		return m_matrixType;
	    }

        private:

            StdMatrixKey();
	    
            ShapeType   m_shapeType;
            boost::shared_array<LibUtilities::BasisSharedPtr> m_base;

            unsigned int m_ncoeffs;
            MatrixType   m_matrixType;
        };


    } // end of namespace
} // end of namespace

#endif //STDEXPMANAGERSS_H

/**
* $Log: StdExpManagers.h,v $
***/
