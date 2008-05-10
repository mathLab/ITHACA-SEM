///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalBndryMap2D.h
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
// Description: Local to Global mapping routines in 1D, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOBNDMAP2D_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOBNDMAP2D_H

#include <MultiRegions/LocalToGlobalBndryMap.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {
        
    class LocalToGlobalBndryMap2D: 
        public LocalToGlobalBndryMap
        {
        public:
            LocalToGlobalBndryMap2D(){};
            
            LocalToGlobalBndryMap2D(const StdRegions::StdExpansionVector &locexp, 
                                    const SpatialDomains::MeshGraph1D &graph2D,
                                    const ConstArray<OneD,LocalRegions::ExpList1DSharedPtr> &bndCondExp,
                                    const ConstArray<OneD,SpatialDomains::BoundaryConditionType> &bndCondTypes);
            
            virtual ~LocalToGlobalBndryMap2D();
            
            inline void ContToLocalBnd(const DNekVec &cont, DNekVec &loc, int offset = 0)
            {
                ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
                ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");
                
                // offset input data by length "offset" for Dirichlet boundary conditions.
                Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);
                Vmath::Vcopy(cont.GetDimension(),cont.GetPtr(),1,&tmp[offset],1);
                
                Vmath::Gathr(m_totLocBndDofs,&tmp[0],&m_locToContBndMap[0], loc.GetPtr());
            }           
            
            inline void AssembleBnd(const DNekVec &loc, DNekVec &cont, int offset = 0)
            {
                ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
                ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");
                Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);
                
                Vmath::Assmb(m_totLocBndDofs,loc.GetPtr(), &m_locToContBndMap[0], &tmp[0]);
                Vmath::Vcopy(cont.GetDimension(),&tmp[offset],1,cont.GetPtr(),1);
            }
            
        protected:
        
        private:

            Array<OneD,NekDouble> m_bndSign;
            bool  m_sign_change;
            int m_totEdgeDofs;
            Array<OneD,NekDouble> m_edgeToElmtMap;
            
           virtual void v_ContToLocalBnd(const DNekVec &cont, DNekVec &loc, int offset = 0)
            {
                ContToLocalBnd(cont,loc,offset);
            }
                        
            virtual void v_AssembleBnd(const DNekVec &loc, DNekVec &cont, int offset = 0)
            {
                AssembleBnd(loc,cont,offset);
            }

            virtual NekDouble v_GetBndSign(int i) 
            {
                return 1.0;
            }
            
        };
        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP2D_H


/** $Log: LocalToGlobalBndryMap2D.h,v $
/** */

