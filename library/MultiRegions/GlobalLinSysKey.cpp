///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysKey.cpp
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
// Description: Definition of GlobalLinSysKey 
//
///////////////////////////////////////////////////////////////////////////////


#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/GlobalLinSysKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Register Mass Matrix creator. 
        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const NekDouble    factor1, 
                                         const NekDouble    factor2, 
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_linSysType(matrixType),
            m_factor1(factor1),
            m_factor2(factor2)
        {
        }

        GlobalLinSysKey::GlobalLinSysKey(const GlobalLinSysKey &key):
            m_solnType(key.m_solnType),
            m_linSysType(key.m_linSysType),
            m_factor1(key.m_factor1),
            m_factor2(key.m_factor2)
        {
        }

        bool operator<(const GlobalLinSysKey &lhs, const GlobalLinSysKey &rhs)
        {
            if(lhs.m_linSysType < rhs.m_linSysType)
            {
                return true;
            }

            if(lhs.m_linSysType > rhs.m_linSysType)
            {
                return false;
            }

            if(lhs.m_solnType < rhs.m_solnType)
            {
                return true;
            }

            if(lhs.m_solnType > rhs.m_solnType)
            {
                return false;
            }

            if(lhs.m_factor1 > rhs.m_factor1)
            {
                return false;
            }

            if(lhs.m_factor1 < rhs.m_factor1)
            {
                return true;
            }


            if(lhs.m_factor2 > rhs.m_factor2)
            {
                return false;
            }

            if(lhs.m_factor2 < rhs.m_factor2)
            {
                return true;
            }

            return false;
        }

        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs)
        {
            os << "MatrixType: " << rhs.GetLinSysType() << ", factor1 (ScaleFactor): "
               <<rhs.GetFactor1() << ", Solution Type: "
               << GlobalSysSolnTypeMap[rhs.GetGlobalSysSolnType()]  << 
                ", factor2 (tau value): " << rhs.GetFactor2()  << std::endl;
            
            return os;
        }
    }
}

/**
* $Log: GlobalLinSysKey.cpp,v $
* Revision 1.2  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.1  2007/07/19 20:02:26  sherwin
* Generalised global matrix solver
*
***/

