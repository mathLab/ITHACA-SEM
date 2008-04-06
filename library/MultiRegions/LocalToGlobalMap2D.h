///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap2D.h
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
// Description: Local to Global mapping routines in 2D, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP2D_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP2D_H

#include <MultiRegions/LocalToGlobalMap.h>
#include <MultiRegions/ExpList1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace MultiRegions
    {
        
    class LocalToGlobalMap2D: 
        public LocalToGlobalMap
        {
        public:
            LocalToGlobalMap2D();

            LocalToGlobalMap2D(const int loclen, 
                               const StdRegions::StdExpansionVector &locexp, 
                               const SpatialDomains::MeshGraph2D &graph2D);

            LocalToGlobalMap2D(const int loclen, 
                                                   const StdRegions::StdExpansionVector &locexp, 
                                                   const SpatialDomains::MeshGraph2D &graph2D,
                                                   const Array<OneD, const MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                                   const Array<OneD, const SpatialDomains::BoundaryConditionType> &bndCondTypes);
            
            virtual ~LocalToGlobalMap2D();
            
            inline NekDouble GetSign(int i) 
            {
                if(m_sign_change)
                {
                    return m_sign[i];
                }
                else
                {
                    return 1.0;
                }
            }

            inline NekDouble GetBndSign(int i) 
            {
                if(m_sign_change)
                {
                    return m_bndSign[i];
                }
                else
                {
                    return 1.0;
                }
            }
            
            inline void LocalToCont(const Array<OneD, const NekDouble> &loc, 
                                    Array<OneD, NekDouble> &cont)
            {
                Array<OneD, NekDouble> tmp(m_totLocDofs);
                if(m_sign_change)
                {
                    Vmath::Vmul(m_totLocDofs,&m_sign[0],1,&loc[0],1,&tmp[0],1);
                    Vmath::Scatr(m_totLocDofs, &tmp[0],&m_locToContMap[0],&cont[0]);
                }
                else
                {                
                    Vmath::Scatr(m_totLocDofs, &loc[0],&m_locToContMap[0],&cont[0]);
                }
            }
            
            inline void ContToLocal(const Array<OneD, const NekDouble> &cont, 
                                    Array<OneD, NekDouble> &loc)
            {
                Vmath::Gathr(m_totLocDofs,&cont[0],&m_locToContMap[0], &loc[0]);
                
                if(m_sign_change)
                {
                    Vmath::Vmul(m_totLocDofs,&m_sign[0],1,&loc[0],1,&loc[0],1);
                }
            }
            
            inline void Assemble(const Array<OneD, const NekDouble> &loc, 
                                 Array<OneD, NekDouble> &cont)
            {
                Vmath::Zero(m_totGloDofs,&cont[0],1);
                
                if(m_sign_change)
                {
                    Array<OneD, NekDouble> tmp(m_totLocDofs);
                    Vmath::Vmul(m_totLocDofs,&m_sign[0],1,&loc[0],1,&tmp[0],1);
                    Vmath::Assmb(m_totLocDofs,&tmp[0],&m_locToContMap[0],&cont[0]);
                }
                else
                {
                    Vmath::Assmb(m_totLocDofs,&loc[0],&m_locToContMap[0],&cont[0]);
                }
            }    

            inline void ContToLocalBnd(const DNekVec &cont, DNekVec &loc, int offset = 0)
            {
                ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
                ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");
                
                // offset input data by length "offset" for Dirichlet boundary conditions.
                Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);
                Vmath::Vcopy(cont.GetDimension(),cont.GetRawPtr(),1,&tmp[offset],1);
                
                Vmath::Gathr(m_totLocBndDofs,&tmp[0],&m_locToContBndMap[0], loc.GetRawPtr());

                if(m_sign_change)
                {
                    Vmath::Vmul(m_totLocBndDofs,&m_bndSign[0],1,loc.GetRawPtr(),1,loc.GetRawPtr(),1);
                }
            }           
            
            inline void AssembleBnd(const DNekVec &loc, DNekVec &cont, int offset = 0)
            {
                ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
                ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");
                Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);

                if(m_sign_change)
                {
                    Array<OneD, NekDouble> tmp2(m_totLocBndDofs);
                    Vmath::Vmul(m_totLocBndDofs,&m_bndSign[0],1,loc.GetRawPtr(),1,&tmp2[0],1);
                    Vmath::Assmb(m_totLocBndDofs,&tmp2[0], &m_locToContBndMap[0], &tmp[0]);
                }
                else
                {                
                    Vmath::Assmb(m_totLocBndDofs,loc.GetRawPtr(), &m_locToContBndMap[0], &tmp[0]);
                }
                Vmath::Vcopy(cont.GetDimension(),&tmp[offset],1,cont.GetRawPtr(),1);
            }        
            
        protected:
            
        private:
            Array<OneD,NekDouble> m_sign;
            Array<OneD,NekDouble> m_bndSign;
            bool  m_sign_change;

            virtual void v_LocalToCont(const Array<OneD, const NekDouble> &loc, 
                                     Array<OneD, NekDouble> &cont)
            {
                LocalToCont(loc,cont);
            }                
            
            
            virtual void v_ContToLocal(const Array<OneD, const NekDouble> &cont, 
                                     Array<OneD, NekDouble> &loc)
            {
                ContToLocal(cont,loc);
            }
            
            virtual void v_Assemble(const Array<OneD, const NekDouble> &loc, 
                                        Array<OneD, NekDouble> &cont)
            {
                Assemble(loc,cont);
            }

            virtual void v_ContToLocalBnd(const DNekVec &cont, DNekVec &loc, int offset = 0)
            {
                ContToLocalBnd(cont,loc,offset);
            }
            
            
            virtual void v_AssembleBnd(const DNekVec &loc, DNekVec &cont, int offset = 0)
            {
                AssembleBnd(loc,cont,offset);
            }

            virtual NekDouble v_GetSign(int i) 
            {
                return GetSign(i);
            }

            virtual NekDouble v_GetBndSign(int i) 
            {
                return GetBndSign(i);
            }
        };
        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP2D_H


/** $Log: LocalToGlobalMap2D.h,v $
 Revision 1.11  2008/04/06 06:00:07  bnelson
 Changed ConstArray to Array<const>

 Revision 1.10  2008/03/18 14:14:13  pvos
 Update for nodal triangular helmholtz solver

 Revision 1.9  2008/01/25 05:50:46  bnelson
 Changed NekVector::GetPtr to NekVector::GetRawPtr and added a new NekVector::GetPtr that returns an Array.  This makes the calls consistent with NekMatrix.

 Revision 1.8  2008/01/20 16:31:11  bnelson
 Fixed linux compile errors.

 Revision 1.7  2007/12/06 22:52:30  pvos
 2D Helmholtz solver updates

 Revision 1.6  2007/10/03 11:37:51  sherwin
 Updates relating to static condensation implementation

 Revision 1.5  2007/07/22 23:04:22  bnelson
 Backed out Nektar::ptr.

 Revision 1.4  2007/07/20 02:04:13  bnelson
 Replaced boost::shared_ptr with Nektar::ptr

 Revision 1.3  2007/06/07 15:54:19  pvos
 Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 Also made corrections to various ASSERTL2 calls

 Revision 1.2  2007/06/05 16:36:55  pvos
 Updated Explist2D ContExpList2D and corresponding demo-codes

 Revision 1.1  2006/07/02 17:16:17  sherwin

 Modifications to make MultiRegions work for a connected domain in 2D (Tris)

 Revision 1.3  2006/06/05 00:14:33  bnelson
 Fixed a compiler error (couldn't find boost::shared_ptr<) and a couple of formatting updates for the standard.
 */

