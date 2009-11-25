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


#include <MultiRegions/GlobalLinSysKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Register Mass Matrix creator. 
        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,locToGloMap))
        {
        }
        
        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap, 
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,varcoeffs,locToGloMap))
        {
        } 

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const NekDouble factor,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,factor,locToGloMap))
        {
        }
        
        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const NekDouble factor1,
                                         const NekDouble factor2,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,factor1,factor2,locToGloMap))
        {
        }

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap, 
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,varcoeffs,locToGloMap))
        {
        }          
            
        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const NekDouble factor,
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,factor,varcoeffs,locToGloMap))
        {
        }     

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const NekDouble factor1,
                                         const NekDouble factor2, 
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,factor1,factor2,varcoeffs,locToGloMap))
        {
        }   

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                                         const Array<OneD, NekDouble> &factor1,
                                         const NekDouble factor2, 
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const GlobalSysSolnType solnType):
            m_solnType(solnType),
            m_globMatKey(MemoryManager<GlobalMatrixKey>::AllocateSharedPtr(matrixType,factor1,factor2,varcoeffs,locToGloMap))
        {
        }   

        GlobalLinSysKey::GlobalLinSysKey(const GlobalLinSysKey &key):
            m_solnType(key.m_solnType),
            m_globMatKey(key.m_globMatKey)            
        {
        }

        bool operator<(const GlobalLinSysKey &lhs, const GlobalLinSysKey &rhs)
        {
            if(lhs.m_solnType < rhs.m_solnType)
            {
                return true;
            }

            if(lhs.m_solnType > rhs.m_solnType)
            {
                return false;
            }

            if( *(lhs.m_globMatKey) < *(rhs.m_globMatKey))
            {
                return true;
            }

            if( *(rhs.m_globMatKey) < *(lhs.m_globMatKey))
            {
                return false;
            }

            return false;
        }

        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs)
        {
            int i;
            os << "MatrixType: " << rhs.GetMatrixType() << endl;
            os << "Solution Type: " << GlobalSysSolnTypeMap[rhs.GetGlobalSysSolnType()] << endl;
            os << "Number of constants: " << rhs.GetNconstants() << endl;
            for(i = 0; i < rhs.GetNconstants();i++) 
            {
                os << "  Constant " << i << ": " << rhs.GetConstant(i) << endl;
            }
            os << "Number of variable coefficients: " << rhs.GetNvariableCoefficients() << endl;
            
            return os;
        }
    }
}

/**
* $Log: GlobalLinSysKey.cpp,v $
* Revision 1.9  2009/11/07 21:11:30  sehunchun
* Variable coefficients parameters are added
*
* Revision 1.8  2009/07/09 21:39:18  sehunchun
* Add another constructor which deals with varcoeffs
*
* Revision 1.7  2009/03/23 10:51:52  pvos
* Added BlockMatrix support
*
* Revision 1.6  2009/02/08 09:10:15  sherwin
* Added member of LocalToGlobalBaseMap so that we can discern matrices of different boundary condition type
*
* Revision 1.5  2008/11/21 10:36:17  pvos
* Added (limited) support for matrix types: ePOSITIVE_DEFINITE_SYMMETRIC and ePOSITIVE_DEFINITE_SYMMETRIC_BANDED
*
* Revision 1.4  2008/11/19 16:02:33  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.3  2007/11/20 16:27:16  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
* Revision 1.2  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.1  2007/07/19 20:02:26  sherwin
* Generalised global matrix solver
*
***/

