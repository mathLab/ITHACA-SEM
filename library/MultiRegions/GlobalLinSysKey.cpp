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
        /**
         * @class GlobalLinSysKey
         * 
         * This class represents a global linear system and is in essence a
         * wrapper around a global matrix key, augmenting it with a specific
         * solution type from GlobalSysSolnType. Each constructor accepts a
         * MatrixType, describing the matrix to be constructed, a
         * LocalToGlobalBaseMap, defining the mapping from the local elemental
         * expansions to a global system, and a GlobalSysSolnType, defining the
         * type of solution (e.g. full matrix, static condenstation). Some
         * constructors include additional parameters for customising the 
         * global operator matrix.
         */

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varCoeffs) :
            GlobalMatrixKey(matrixType, locToGloMap, factors, varCoeffs),
            m_solnType(locToGloMap->GetGlobalSysSolnType())
        {

        }

//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                                    ::AllocateSharedPtr(matrixType,locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   varcoeffs   Matrix of coefficients.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                         ::AllocateSharedPtr(matrixType,varcoeffs,locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   factor      Scalar coefficient.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                        const NekDouble factor):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                            ::AllocateSharedPtr(matrixType,factor,locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   factor1     Scalar coefficient.
//         * @param   factor2     Scalar coefficient.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                        const NekDouble factor1,
//                        const NekDouble factor2):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                            ::AllocateSharedPtr(matrixType,factor1,factor2,
//                                                locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   factor      Scalar coefficient.
//         * @param   varcoeffs   Matrix of coefficients.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                        const NekDouble factor,
//                        const Array<OneD, Array<OneD, NekDouble> >& varcoeffs):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                            ::AllocateSharedPtr(matrixType,factor,varcoeffs,
//                                                locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   factor1     Scalar factor.
//         * @param   factor2     Scalar factor.
//         * @param   varcoeffs   Matrix of coefficients.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                        const NekDouble factor1,
//                        const NekDouble factor2,
//                        const Array<OneD, Array<OneD,NekDouble> >& varcoeffs):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                            ::AllocateSharedPtr(matrixType,factor1,factor2,
//                                                varcoeffs,locToGloMap))
//        {
//        }
//
//
//        /**
//         * @param   matrixType  Specify the type of global matrix to construct.
//         * @param   locToGloMap Mapping from local elements to the global
//         *                      system.
//         * @param   solnType    Type of solution to construct from
//         *                      GlobalSysSolnType.
//         */
//        GlobalLinSysKey::GlobalLinSysKey(
//                        const StdRegions::MatrixType matrixType,
//                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                        const Array<OneD, NekDouble> &factor1,
//                        const NekDouble factor2,
//                        const Array<OneD, Array<OneD,NekDouble> >& varcoeffs):
//            m_solnType(locToGloMap->GetGlobalSysSolnType()),
//            m_globMatKey(MemoryManager<GlobalMatrixKey>
//                            ::AllocateSharedPtr(matrixType,factor1,factor2,
//                                                varcoeffs,locToGloMap))
//        {
//        }


        /**
         * @param   key         Existing key to duplicate.
         */
        GlobalLinSysKey::GlobalLinSysKey(const GlobalLinSysKey &key):
            m_solnType(key.m_solnType),
            GlobalMatrixKey(key)
        {
        }


        /**
         *
         */
        GlobalLinSysKey::~GlobalLinSysKey()
        {
        }
        
        
        /**
         * Compares two GlobalLinSysKeys by comparing their solution types and
         * matrix keys.
         * @param   lhs         First operand.
         * @param   rhs         Second operand.
         * @returns true if the first operand is considered less than the
         *          second operand.
         */
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

            return (*dynamic_cast<const GlobalMatrixKey*>(&lhs)
                    < *dynamic_cast<const GlobalMatrixKey*>(&rhs));
        }


        /**
         * Writes the vital statistics of a global linear system to a stream.
         * @param   os          Output stream.
         * @param   rhs         GlobalLinSys object to use.
         * @returns Reference to the output stream \a os.
         */
        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs)
        {
            int i;
            os << "MatrixType: " << StdRegions::MatrixTypeMap[rhs.GetMatrixType()] << ", ShapeType: "
                            << StdRegions::ExpansionTypeMap[rhs.GetExpansionType()]
                            << std::endl;
            os << "Solution Type: " 
               << GlobalSysSolnTypeMap[rhs.GetGlobalSysSolnType()] << endl;
            os << "Number of constants: " << rhs.GetNConstFactors() << endl;
            StdRegions::ConstFactorMap::const_iterator x;
            for (x = rhs.GetConstFactors().begin(); x != rhs.GetConstFactors().end(); ++x)
            {
                os << "  Constant " << StdRegions::ConstFactorTypeMap[x->first]
                   << ": " << x->second << endl;
            }
            os << "Number of variable coefficients: " 
               << rhs.GetNVarCoeffs() << endl;
            
            return os;
        }
    }
}

/**
* $Log: GlobalLinSysKey.cpp,v $
* Revision 1.10  2009/11/25 17:15:45  sehunchun
* Add a function when factor1 is a vector
*
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

