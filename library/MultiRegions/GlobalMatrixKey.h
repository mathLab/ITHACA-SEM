///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrixKey.h
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
// Description: Headers for GlobalMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H
#define NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Describes a matrix with ordering defined by a local to global map.
        class GlobalMatrixKey
        {
        public:
            /// Matrix without any parameters.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with a single real parameter.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with two real parameters.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor1,
                            const NekDouble factor2,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with a single variable coefficient parameter.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const Array<OneD,NekDouble>& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with multiple variable coefficient parameters.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with a single real parameter and multiple variable
            /// coefficient parameters.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor,
                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Matrix with two real parameters and multiple variable
            /// coefficient parameters.
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor1,
                            const NekDouble factor2,
                            const Array<OneD,Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap
                                        = NullLocalToGlobalBaseMapSharedPtr);

            /// Copy constructor.
            GlobalMatrixKey(const GlobalMatrixKey &key);

            /// Destructor
            ~GlobalMatrixKey();

            /// Provides ordering of GlobalMatrixKey objects.
            friend bool operator<(const GlobalMatrixKey &lhs,
                                  const GlobalMatrixKey &rhs);

            /// Return the matrix type.
            const StdRegions::MatrixType GetMatrixType() const;
            /// Returns true if a local to global map is defined.
            const bool LocToGloMapIsDefined() const;
            /// Returns the number of constants defined for this matrix.
            int GetNconstants() const;
            /// Returns the requested constant.
            NekDouble GetConstant(int i) const;
            /// Returns all the constants.
            const Array<OneD,NekDouble>& GetConstants() const;
            /// Returns the number of variable coefficient fields.
            int GetNvariableCoefficients() const;
            /// Returns the requested variable coefficient parameter.
            const Array<OneD,NekDouble>& GetVariableCoefficient(int i) const;
            /// Returns all the variable coefficient parameters.
            const Array<OneD, Array<OneD,NekDouble> >&
                                            GetVariableCoefficients() const;

        protected:
            /// Default constructor.
            GlobalMatrixKey();

            /// Stores the matrix type based on the enum StdRegions::MatrixType.
            StdRegions::MatrixType m_matrixType;

            /// The number of real parameters for the matrix.
            int                   m_nconstants;
            /// The real parameters for the matrix.
            Array<OneD,NekDouble> m_constant;

            /// The number of variable coefficients for the matrix.
            int m_nvariablecoefficients;
            /// The variable coefficients for the matrix.
            Array<OneD, Array<OneD,NekDouble> >  m_variablecoefficient;

            /// Pointer to the local to global mapping.
            LocalToGlobalBaseMapSharedPtr m_locToGloMap;

        private:

        };

        /// Writes statistics about the matrix key to an output stream.
        std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs);

        /// A pointer to a GlobalMatrixKey.
        typedef  boost::shared_ptr<GlobalMatrixKey> GlobalMatrixKeySharedPtr;

        inline const StdRegions::MatrixType
                        GlobalMatrixKey::GetMatrixType() const
        {
            return m_matrixType;
        }

        inline const bool GlobalMatrixKey::LocToGloMapIsDefined(void) const
        {
            if( m_locToGloMap.get() == 0) //NullLocalToGlobalBaseMapSharedPtr)
            {
                return false;
            }

            return true;
        }

        inline int GlobalMatrixKey::GetNconstants() const
        {
            return m_nconstants;
        }

        inline NekDouble GlobalMatrixKey::GetConstant(int i) const
        {
            ASSERTL1(i < m_nconstants,
                     "requesting constant which has not been definied");
            return m_constant[i];
        }

        inline const Array<OneD,NekDouble>&
                        GlobalMatrixKey::GetConstants() const
        {
            return m_constant;
        }

        inline int GlobalMatrixKey::GetNvariableCoefficients() const
        {
            return m_nvariablecoefficients;
        }

        inline const Array<OneD,NekDouble>&
                        GlobalMatrixKey::GetVariableCoefficient(int i) const
        {
            ASSERTL1(i < m_nvariablecoefficients,
                     "requesting a coefficient which has not been defined");
            return m_variablecoefficient[i];
        }

        inline const Array<OneD, Array<OneD,NekDouble> >&
                        GlobalMatrixKey::GetVariableCoefficients() const
        {
            return m_variablecoefficient;
        }


    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H

/**
* $Log: GlobalMatrixKey.h,v $
* Revision 1.3  2009/11/19 23:30:36  cantwell
* Documentation for ExpList2D and GlobalMatrixKey
* Updated doxygen pages.
*
* Revision 1.2  2009/11/07 21:11:30  sehunchun
* Variable coefficients parameters are added
*
* Revision 1.1  2009/03/23 10:46:41  pvos
* Added GlobalMatrixKey
*
**/
