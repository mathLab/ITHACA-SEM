///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysKeys.h
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
// Description: Headers for GlobalLinSysKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GLOBALLINSYSKEY_H
#define GLOBALLINSYSKEY_H

#include <MultiRegions/MultiRegions.hpp>

namespace Nektar
{
    namespace MultiRegions
    {

        class GlobalLinSysKey
        {
        public:
            GlobalLinSysKey(const StdRegions::MatrixType matrixType, 
                            const double scalefactor = 1.0,
#if 0
                            const GlobalSysSolnType solnType = eDirectFullMatrix);
#else
                            const GlobalSysSolnType solnType = eDirectStaticCond);
#endif

            GlobalLinSysKey(const GlobalLinSysKey &key);

            ~GlobalLinSysKey()
            {
            }

            friend bool operator<(const GlobalLinSysKey &lhs, 
                                  const GlobalLinSysKey &rhs);

            const StdRegions::MatrixType GetLinSysType() const
            {
                return m_linSysType; 
            }

            const GlobalSysSolnType  GetGlobalSysSolnType() const
            {
                return m_solnType; 
            }
            
            const NekDouble GetScaleFactor() const 
            {
                return m_scaleFactor;
            }

        protected:
            GlobalLinSysKey(); 
            GlobalSysSolnType      m_solnType;
            StdRegions::MatrixType m_linSysType;
            NekDouble              m_scaleFactor;
            
        private:
        };

        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs);

    } // end of namespace
} // end of namespace

#endif //GLOBALMATRIXKEY_H

/**
* $Log: GlobalLinSysKey.h,v $
* Revision 1.1  2007/07/19 20:02:26  sherwin
* Generalised global matrix solver
*
***/
