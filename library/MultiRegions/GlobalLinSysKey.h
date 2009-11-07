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

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class GlobalLinSysKey
        {
        public:
            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const NekDouble factor,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const NekDouble factor1,
                            const NekDouble factor2,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const NekDouble factor,
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                            const NekDouble factor1,
                            const NekDouble factor2, 
                            const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                            const GlobalSysSolnType solnType);

            GlobalLinSysKey(const GlobalLinSysKey &key);

            ~GlobalLinSysKey()
            {
            }

            friend bool operator<(const GlobalLinSysKey &lhs, 
                                  const GlobalLinSysKey &rhs);

            inline const GlobalMatrixKeySharedPtr& GetGlobalMatrixKey() const
            {
                return m_globMatKey;
            }

            inline const StdRegions::MatrixType GetMatrixType() const
            {
                return m_globMatKey->GetMatrixType(); 
            }

            inline const GlobalSysSolnType  GetGlobalSysSolnType() const
            {
                return m_solnType; 
            }

            inline const bool LocToGloMapIsDefined(void) const
            {
                return m_globMatKey->LocToGloMapIsDefined();
            }
            
            inline int GetNconstants() const
            {
                return m_globMatKey->GetNconstants();
            }

            inline NekDouble GetConstant(int i) const
            {
                return m_globMatKey->GetConstant(i);
            }

            inline const Array<OneD,NekDouble>& GetConstants() const
            {         
                return m_globMatKey->GetConstants();
            }

            inline int GetNvariableCoefficients() const
            {
                return m_globMatKey->GetNvariableCoefficients();
            }

            inline const Array<OneD,NekDouble>& GetVariableCoefficient(int i) const
            {
                return m_globMatKey->GetVariableCoefficient(i);
            }

            inline const Array<OneD, Array<OneD,NekDouble> >& GetVariableCoefficients() const
            {       
                return m_globMatKey->GetVariableCoefficients();
            }

        protected:
            GlobalLinSysKey(); 

            GlobalSysSolnType        m_solnType;
            GlobalMatrixKeySharedPtr m_globMatKey;
            
        private:
        };

        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs);

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H

/**
* $Log: GlobalLinSysKey.h,v $
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
