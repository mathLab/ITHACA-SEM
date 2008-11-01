///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalBaseMap.h
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
// Description: Local to Global base mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALBASEMAP_H
#define MULTIREGIONS_LOCALTOGLOBALBASEMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        class LocalToGlobalBaseMap
        {
        public:
            LocalToGlobalBaseMap(); 
            
            ~LocalToGlobalBaseMap(); 
                
            int GetLocalToGlobalBndMap(const int i) const
            {
                return m_localToGlobalBndMap[i];
            }
            
            inline const Array<OneD,const int>&  GetLocalToGlobalBndMap(void)
            {
                return m_localToGlobalBndMap;
            }

            NekDouble GetLocalToGlobalBndSign(const int i) const
            {
                if(m_signChange)
                {
                    return m_localToGlobalBndSign[i];
                }
                else
                {
                    return 1.0;
                }
            }
            
            inline int GetBndCondCoeffsToGlobalCoeffsMap(const int i)
            {
                return m_bndCondCoeffsToGlobalCoeffsMap[i];
            }
            
            inline NekDouble GetBndCondCoeffsToGlobalCoeffsSign(const int i)
            {
                if(m_signChange)
                {
                return m_bndCondCoeffsToGlobalCoeffsSign[i];
                }
                else
                {
                    return 1.0;
                }
            }
                        
            inline const Array<OneD,const int>& GetBndCondCoeffsToGlobalCoeffsMap(void)
            {
                return m_bndCondCoeffsToGlobalCoeffsMap;
            }

            inline int GetNumDirichletBndCoeffs() const
            {
                return m_numDirichletBndCoeffs;
            }

            inline int GetNumLocalBndCoeffs() const
            {
                return m_numLocalBndCoeffs;
            }

            inline int GetNumGlobalBndCoeffs() const
            {
                return m_numGlobalBndCoeffs;
            }



            inline void GlobalToLocalBnd(const NekVector<const NekDouble>& global, NekVector<NekDouble>& loc, int offset)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");
                
                // offset input data by length "offset" for Dirichlet boundary conditions.
                Array<OneD,NekDouble> tmp(global.GetDimension()+offset,0.0);
                Vmath::Vcopy(global.GetDimension(), global.GetRawPtr(), 1, tmp.get() + offset, 1);
                
                Vmath::Gathr(m_numLocalBndCoeffs, tmp.get(), m_localToGlobalBndMap.get(), loc.GetRawPtr());

                if(m_signChange)
                {
                    Vmath::Vmul(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), 1, loc.GetRawPtr(), 1, loc.GetRawPtr(), 1);
                }
            }     

            inline void GlobalToLocalBnd(const NekVector<const NekDouble>& global, NekVector<NekDouble>& loc)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");
                
                Vmath::Gathr(m_numLocalBndCoeffs, global.GetRawPtr(), m_localToGlobalBndMap.get(), loc.GetRawPtr());

                if(m_signChange)
                {
                    Vmath::Vmul(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), 1, loc.GetRawPtr(), 1, loc.GetRawPtr(), 1);
                }
            } 
            
            inline void AssembleBnd(const NekVector<const NekDouble>& loc, NekVector<NekDouble>& global, int offset)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");
                Array<OneD,NekDouble> tmp(global.GetDimension()+offset,0.0);

                if(m_signChange)
                {
                    Array<OneD, NekDouble> tmp2(m_numLocalBndCoeffs);
                    Vmath::Vmul(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), 1, loc.GetRawPtr(), 1, tmp2.get(), 1);
                    Vmath::Assmb(m_numLocalBndCoeffs, tmp2.get(), m_localToGlobalBndMap.get(), tmp.get());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.GetRawPtr(), m_localToGlobalBndMap.get(), tmp.get());
                }
                Vmath::Vcopy(global.GetDimension(), tmp.get() + offset, 1, global.GetRawPtr(), 1);
            }   
            
            inline void AssembleBnd(const NekVector<const NekDouble>& loc, NekVector<NekDouble>& global)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

                Vmath::Zero(m_numGlobalBndCoeffs, global.GetRawPtr(), 1);

                if(m_signChange)
                {
                    Array<OneD, NekDouble> tmp2(m_numLocalBndCoeffs);
                    Vmath::Vmul(m_numLocalBndCoeffs,  m_localToGlobalBndSign.get(), 1, loc.GetRawPtr(), 1, tmp2.get(), 1);
                    Vmath::Assmb(m_numLocalBndCoeffs, tmp2.get(), m_localToGlobalBndMap.get(), global.GetRawPtr());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.GetRawPtr(), m_localToGlobalBndMap.get(), global.GetRawPtr());
                }
            }    

        protected:
            int m_numLocalBndCoeffs;      //< number of local Bnd coefficients
            int m_numGlobalBndCoeffs;     // Total number of global boundary coefficients
            int m_numDirichletBndCoeffs;  // Number of Dirichlet Boundary Coefficient
            bool m_signChange;
            Array<OneD,int> m_localToGlobalBndMap;  //< integer map of local boundary coeffs to global space 
            Array<OneD,NekDouble> m_localToGlobalBndSign; //< integer sign of local boundary coeffs to global space 
 
            Array<OneD,int> m_bndCondCoeffsToGlobalCoeffsMap;  //< integer map of bnd cond coeffs to global coefficients
            Array<OneD,NekDouble> m_bndCondCoeffsToGlobalCoeffsSign;  //< integer map of bnd cond coeffs to global coefficients
        private:
        };
        typedef boost::shared_ptr<LocalToGlobalBaseMap>  LocalToGlobalBaseMapSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** $Log: LocalToGlobalBaseMap.h,v $
 Revision 1.4  2008/09/23 18:21:00  pvos
 Updates for working ProjectContField3D demo

 Revision 1.3  2008/09/17 13:46:40  pvos
 Added LocalToGlobalC0ContMap for 3D expansions

 Revision 1.2  2008/09/16 13:36:06  pvos
 Restructured the LocalToGlobalMap classes

 Revision 1.1  2008/08/18 08:16:23  sherwin
 First version of this new class container for mappings

 */

