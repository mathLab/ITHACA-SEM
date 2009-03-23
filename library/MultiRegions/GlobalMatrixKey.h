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
        class GlobalMatrixKey
        {
        public:
            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap = NullLocalToGlobalBaseMapSharedPtr);

            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap = NullLocalToGlobalBaseMapSharedPtr);

            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor1,
                            const NekDouble factor2,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap = NullLocalToGlobalBaseMapSharedPtr);

            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap = NullLocalToGlobalBaseMapSharedPtr);

            GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const NekDouble factor,
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap = NullLocalToGlobalBaseMapSharedPtr);

            GlobalMatrixKey(const GlobalMatrixKey &key);

            ~GlobalMatrixKey()
            {
            }

            friend bool operator<(const GlobalMatrixKey &lhs, 
                                  const GlobalMatrixKey &rhs);

            inline const StdRegions::MatrixType GetMatrixType() const
            {
                return m_matrixType; 
            }

            inline const bool LocToGloMapIsDefined(void) const
            {
                if( m_locToGloMap.get() == 0) //NullLocalToGlobalBaseMapSharedPtr)
                {
                    return false;
                }

                return true;
            }
            
            inline int GetNconstants() const
            {
                return m_nconstants;
            }

            inline NekDouble GetConstant(int i) const
            {
                ASSERTL1(i < m_nconstants,"requesting constant which has not been definied");
                return m_constant[i];
            }

            inline const Array<OneD,NekDouble>& GetConstants() const
            {         
                return m_constant;
            }

            inline int GetNvariableCoefficients() const
            {
                return m_nvariablecoefficients;
            }

            inline const Array<OneD,NekDouble>& GetVariableCoefficient(int i) const
            {
                ASSERTL1(i < m_nvariablecoefficients,"requesting a coefficient which has not been defined");                
                return m_variablecoefficient[i];
            }

            inline const Array<OneD, Array<OneD,NekDouble> >& GetVariableCoefficients() const
            {       
                return m_variablecoefficient;
            }

        protected:
            GlobalMatrixKey(); 

            StdRegions::MatrixType m_matrixType;
            
            int                   m_nconstants;
            Array<OneD,NekDouble> m_constant;

            int m_nvariablecoefficients;
            Array<OneD, Array<OneD,NekDouble> >  m_variablecoefficient;
            
            LocalToGlobalBaseMapSharedPtr m_locToGloMap; 

        private:
        };

        std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs);

        typedef  boost::shared_ptr<GlobalMatrixKey> GlobalMatrixKeySharedPtr;

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H

/**
* $Log: GlobalMatrixKey.h,v $
**/
