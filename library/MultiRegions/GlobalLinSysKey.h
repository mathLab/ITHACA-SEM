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

#ifndef NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H
#define NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Describe a linear system.
        class GlobalLinSysKey : public GlobalMatrixKey
        {
        public:
            /// Default constructor.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey()
//            {
//            }

            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varCoeffs = StdRegions::NullVarCoeffMap);

//            /// Linear system with no parameters.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap);
//            /// Linear system with variable coefficients
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap);
//            /// Linear system with a single real factor.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                            const NekDouble factor);
//
//            /// Linear system with two real factors.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                            const NekDouble factor1,
//                            const NekDouble factor2);
//            /// Linear system with a real factor and variable coefficients.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                            const NekDouble factor,
//                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs);
//            /// Linear system with two real factors and variable coefficients.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                            const NekDouble factor1,
//                            const NekDouble factor2,
//                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs);
//             /// Linear system with an array of real factors, a second real
//            /// factor and variable coefficients.
//            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
//                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
//                            const Array<OneD,NekDouble> &factor1,
//                            const NekDouble factor2,
//                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs);
            /// Copy constructor.
            MULTI_REGIONS_EXPORT GlobalLinSysKey(const GlobalLinSysKey &key);
            /// Destructor.
            MULTI_REGIONS_EXPORT ~GlobalLinSysKey();

            /// Less-than operator for GlobalLinSysKey comparison.
            MULTI_REGIONS_EXPORT friend bool operator<(const GlobalLinSysKey &lhs, 
                                  const GlobalLinSysKey &rhs);
//            /// Return the associated global matrix key.
//            inline const GlobalMatrixKeySharedPtr& GetGlobalMatrixKey() const;
//            /// Return the type of the associated matrix.
//            inline const StdRegions::MatrixType GetMatrixType() const;
            /// Return the associated solution type.
            inline const GlobalSysSolnType  GetGlobalSysSolnType() const;
//            /// Determine if a local to global mapping exists.
//            inline const bool LocToGloMapIsDefined(void) const;
//            /// Return the number of constants.
//            inline int GetNconstants() const;
//            /// Retrieve a specific constant.
//            inline NekDouble GetConstant(int i) const;
//            /// Retrieve the vector of constants.
//            inline const Array<OneD,NekDouble>& GetConstants() const;
//            /// Returns the number of variable coefficients.
//            inline int GetNvariableCoefficients() const;
//            /// Retrieve a particular variable coefficient.
//            inline const Array<OneD,NekDouble>&
//                                            GetVariableCoefficient(int i) const;
//            /// Retrieve the array of variable coefficients.
//            inline const Array<OneD, Array<OneD,NekDouble> >&
//                                            GetVariableCoefficients() const;

        protected:
            /// Store the solution type associated with the linear system. This
            /// may be none, full matrix, static condensation or multi-level
            /// static condensation.
            GlobalSysSolnType        m_solnType;
            
//            /// Stores a shared pointer to the associated global matrix key.
//            GlobalMatrixKeySharedPtr m_globMatKey;
            
        private:
        
        };

        /// Writes information about the object to a given stream.
        MULTI_REGIONS_EXPORT std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs);

//        inline const GlobalMatrixKeySharedPtr&
//                                    GlobalLinSysKey::GetGlobalMatrixKey() const
//        {
//            return m_globMatKey;
//        }

//        inline const StdRegions::MatrixType GlobalLinSysKey::GetMatrixType()
//                                                                        const
//        {
//            return m_globMatKey->GetMatrixType();
//        }

        inline const GlobalSysSolnType GlobalLinSysKey::GetGlobalSysSolnType() 
                                                                        const
        {
            return m_solnType; 
        }

//        inline const bool GlobalLinSysKey::LocToGloMapIsDefined(void) const
//        {
//            return m_globMatKey->LocToGloMapIsDefined();
//        }
//
//        inline int GlobalLinSysKey::GetNconstants() const
//        {
//            return m_globMatKey->GetNconstants();
//        }
//
//        inline NekDouble GlobalLinSysKey::GetConstant(int i) const
//        {
//            return m_globMatKey->GetConstant(i);
//        }
//
//        inline const Array<OneD,NekDouble>& GlobalLinSysKey::GetConstants()
//                                                                        const
//        {
//            return m_globMatKey->GetConstants();
//        }
//
//        inline int GlobalLinSysKey::GetNvariableCoefficients() const
//        {
//            return m_globMatKey->GetNvariableCoefficients();
//        }
//
//        inline const Array<OneD,NekDouble>&
//                        GlobalLinSysKey::GetVariableCoefficient(int i) const
//        {
//            return m_globMatKey->GetVariableCoefficient(i);
//        }

//        inline const Array<OneD, Array<OneD,NekDouble> >&
//                        GlobalLinSysKey::GetVariableCoefficients() const
//        {
//            return m_globMatKey->GetVariableCoefficients();
//        }

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H

/**
* $Log: GlobalLinSysKey.h,v $
* Revision 1.13  2009/11/25 17:15:45  sehunchun
* Add a function when factor1 is a vector
*
* Revision 1.12  2009/11/07 21:11:30  sehunchun
* Variable coefficients parameters are added
*
* Revision 1.11  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
* Revision 1.10  2009/07/09 21:39:18  sehunchun
* Add another constructor which deals with varcoeffs
*
* Revision 1.9  2009/03/23 10:51:52  pvos
* Added BlockMatrix support
*
* Revision 1.8  2009/02/09 16:12:08  sherwin
* .
*
* Revision 1.7  2009/02/09 16:11:26  sherwin
* made LocToGloMapIsDefined return a const value
*
* Revision 1.6  2009/02/08 09:10:15  sherwin
* Added member of LocalToGlobalBaseMap so that we can discern matrices of different boundary condition type
*
* Revision 1.5  2008/11/19 16:02:33  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.4  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
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
