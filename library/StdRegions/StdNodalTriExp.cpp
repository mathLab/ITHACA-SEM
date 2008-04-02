///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.cpp
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
// Description: Nodal triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTriExp.h>

namespace Nektar
{
    namespace StdRegions
    {
        StdNodalTriExp::StdNodalTriExp(void):
            StdTriExp(),
            m_nodalPointsKey()
        {
        }

        StdNodalTriExp::StdNodalTriExp(const LibUtilities::BasisKey &Ba, 
                                       const LibUtilities::BasisKey &Bb, 
                                       LibUtilities::PointsType Ntype):
            StdTriExp(Ba,Bb),
            m_nodalPointsKey()
        {
            ASSERTL0(m_base[0]->GetNumModes() == m_base[1]->GetNumModes(),
                     "Nodal basis initiated with different orders in the a "
                     "and b directions");   
            int nummodes =  Ba.GetNumModes();
            m_nodalPointsKey = MemoryManager<LibUtilities::PointsKey>::AllocateSharedPtr(nummodes,Ntype);
        }

        StdNodalTriExp::StdNodalTriExp(const StdNodalTriExp &T):
            StdTriExp(T),
            m_nodalPointsKey(T.m_nodalPointsKey)
        {
        }

        // Destructor
        StdNodalTriExp::~StdNodalTriExp()
        { 
        }
        
        DNekMatSharedPtr StdNodalTriExp::GenMatrix(const StdMatrixKey &mkey)
        {
            DNekMatSharedPtr Mat;

            switch(mkey.GetMatrixType())
            {
            case eNBasisTrans:
                Mat =  GenNBasisTransMatrix();
                break;
            default:
                Mat = StdExpansion::CreateGeneralMatrix(mkey);
                break;
            }

            return Mat;
        }

        DNekMatSharedPtr StdNodalTriExp::GenNBasisTransMatrix()
        {
            int             i,j;
            ConstArray<OneD, NekDouble>  r, s; 
            Array<OneD, NekDouble> c(2);
            DNekMatSharedPtr Mat;

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);
            GetNodalPoints(r,s);

            for(i = 0; i < m_ncoeffs; ++i)
            {
                // fill physical space with mode i
                StdTriExp::FillMode(i,m_phys);

                // interpolate mode i to the Nodal points 'j' and
                // store in outarray
                for(j = 0; j < m_ncoeffs; ++j)
                {
                    c[0] = r[j];
                    c[1] = s[j];
                    (*Mat)(j,i) = StdTriExp::PhysEvaluate(c);
                }
            }
            return Mat;
        }

        void StdNodalTriExp::NodalToModal()
        {
            NodalToModal(m_coeffs,m_coeffs); 
        }

        void StdNodalTriExp::NodalToModal(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {
            StdMatrixKey   Nkey(eInvNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<const NekDouble> nodal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = (*inv_vdm) * nodal;
        }


        void StdNodalTriExp::NodalToModalTranspose()
        {
            NodalToModalTranspose(m_coeffs,m_coeffs); 
        }

        // Operate with transpose of NodalToModal transformation
        void StdNodalTriExp::NodalToModalTranspose(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {            
            StdMatrixKey   Nkey(eInvNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<const NekDouble> nodal(m_ncoeffs,inarray,eCopy);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = Transpose(*inv_vdm) * nodal;
        }


        void StdNodalTriExp::ModalToNodal()
        {
            ModalToNodal(m_coeffs,m_coeffs);
        }

        void StdNodalTriExp::ModalToNodal(const ConstArray<OneD, NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {
            StdMatrixKey      Nkey(eNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  vdm = GetStdMatrix(Nkey);

            // Multiply out matrix
            NekVector<const NekDouble> modal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> nodal(m_ncoeffs,outarray,eWrapper);
            nodal = (*vdm)*modal;
        }


        //////////////////////////////
        /// Integration Methods
        //////////////////////////////


        void StdNodalTriExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD,  NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),
                m_base[1]->GetBdata(), inarray, 
                outarray);
        }

        void StdNodalTriExp:: IProductWRTBase(const ConstArray<OneD, NekDouble>& base0,
            const ConstArray<OneD, NekDouble>& base1,
            const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            // Take inner product with respect to Orthgonal basis using
            // StdTri routine
            StdTriExp::IProductWRTBase(base0,base1,inarray,outarray);
            NodalToModalTranspose(outarray,outarray);     
        }

        void StdNodalTriExp::IProductWRTDerivBase(const int dir, 
                                                  const ConstArray<OneD, NekDouble>& inarray, 
                                                  Array<OneD, NekDouble> & outarray)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            
            Array<OneD, NekDouble> gfac0(nqtot);
            Array<OneD, NekDouble> tmp0(nqtot);
            
            ConstArray<OneD, NekDouble> z1 = ExpPointsProperties(1)->GetZ();
            
            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }
            
            for(i = 0; i < nquad1; ++i)  
            {
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
            }
            
            switch(dir)
            {
            case 0:
                {                    
                    IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),
                                    tmp0,outarray);
                }
                break;
            case 1:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);    
                    ConstArray<OneD, NekDouble> z0 = ExpPointsProperties(0)->GetZ();
                    
                    for(i = 0; i < nquad0; ++i)
                    {
                        gfac0[i] = 0.5*(1+z0[i]);
                    }        
                    
                    for(i = 0; i < nquad1; ++i) 
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }       
                    
                    IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp0,tmp3); 
                    IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),inarray,outarray);  
                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,&outarray[0],1);      
                }
                break;
            default:
                {
                    ASSERTL1(dir >= 0 &&dir < 2,"input dir is out of range");
                }
                break;
            }             
        }

        void StdNodalTriExp::FillMode(const int mode, 
            Array<OneD, NekDouble> &outarray)
        {

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs, outarray, 1);
            outarray[mode] = 1.0;
            BwdTrans(outarray,outarray);
        }

        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////

        // Currently convert nodal values into tranformed values and
        // backward transform

        void StdNodalTriExp::BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs);
            NodalToModal(inarray,tmp);
            StdTriExp::BwdTrans(tmp,outarray);
        }

        void StdNodalTriExp::FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);
            
            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this, 
                                    m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            NekVector<const NekDouble> in(m_ncoeffs,outarray,eWrapper);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
            
            out = (*matsys)*in;
        }

        void StdNodalTriExp::GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if(outarray.num_elements()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }
            
            for(i = 0; i < NumBndryCoeffs(); i++)
            {
                outarray[i] = i;
            }
        }

        void StdNodalTriExp::GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if(outarray.num_elements()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }

            for(i = NumBndryCoeffs(); i < GetNcoeffs(); i++)
            {
                outarray[i-NumBndryCoeffs()] = i;
            }
        }
 
        void StdNodalTriExp::GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                           Array<OneD, unsigned int> &maparray)
        {
            ASSERTL0((eid>=0)&&(eid<=2),"Local Edge ID must be between 0 and 2"); 
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;
            
            if(maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }
            
            for(int i = 0; i < nEdgeIntCoeffs; i++)
            {
                maparray[i] = eid*nEdgeIntCoeffs+3+i; 
            }  

            if(edgeOrient == eBackwards)
            {
                reverse( maparray.get() , maparray.get()+nEdgeIntCoeffs );
            }
                          
        }

        void StdNodalTriExp::GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                             Array<OneD, unsigned int> &maparray)
        {
            ASSERTL0((eid>=0)&&(eid<=2),"Local Edge ID must be between 0 and 2"); 
            const int nEdgeCoeffs = GetEdgeNcoeffs(eid);
            
            if(maparray.num_elements() != nEdgeCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
            }
            
            maparray[0] = eid;
            maparray[1] = (eid==2) ? 0 : eid+1;
            for(int i = 2; i < nEdgeCoeffs; i++)
            {
                maparray[i] = eid*(nEdgeCoeffs-2)+1+i; 
            }  

            if(edgeOrient == eBackwards)
            {
                reverse( maparray.get() , maparray.get()+nEdgeCoeffs );
            }
        }

        void  StdNodalTriExp::MapTo(const int edge_ncoeffs, 
            const LibUtilities::BasisType Btype,
            const int eid, 
            const EdgeOrientation eorient, 
            StdExpMap &Map)
        {

            int i;
            int *dir, order0,order1;
            Array<OneD, int> wsp; 

            ASSERTL2(eid>=0&&eid<=2,"eid must be between 0 and 2");

            ASSERTL2(Btype == m_base[0]->GetBasisType(),
                "Expansion type of edge and StdQuadExp are different");

            // make sure have correct memory storage
            if(edge_ncoeffs != Map.GetLen())
            {
                Map.SetMapMemory(edge_ncoeffs);
            }

            order0 = m_base[0]->GetNumModes();
            order1 = m_base[1]->GetNumModes();

            wsp = Array<OneD, int>(edge_ncoeffs);
            dir = wsp.get(); 

            if(eorient == eForwards)
            {
                for(i = 0; i < edge_ncoeffs; ++i)
                {
                    dir[i] = i;
                }
            }
            else
            {
                dir[1] = 0; 
                dir[0] = 1;

                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    dir[i] = edge_ncoeffs-i+1;
                }
            }

            // Set up Mapping details
            switch (eid)
            {
            case 0:
                Map[dir[0]] = 0;
                Map[dir[1]] = 1;


                for(i = 2; i < edge_ncoeffs;  ++i)
                {
                    Map[dir[i]] = i+1; 
                }
                break;
            case 1:
                Map[dir[0]] = 1;
                Map[dir[1]] = 2;

                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = i+edge_ncoeffs-1; 
                }
                break;
            case 2:
                Map[dir[0]] = 0;
                Map[dir[1]] = 2;

                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = i+2*edge_ncoeffs-3; 
                }
                break;
            }
        }
        
        void StdNodalTriExp::MapTo_ModalFormat(const int edge_ncoeffs, 
            const LibUtilities::BasisType Btype, 
            const int eid, 
            const EdgeOrientation eorient,
            StdExpMap &Map)
        {        
            MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
        }


    } // end of namespace
} // end of namespace

/** 
* $Log: StdNodalTriExp.cpp,v $
* Revision 1.20  2008/03/18 14:15:45  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.19  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.18  2007/08/11 23:42:26  sherwin
* A few changes
*
* Revision 1.17  2007/07/27 00:22:54  bnelson
* Memory manager now accepts non-const parameters to the allocate methods.
*
* Revision 1.16  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.15  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.14  2007/07/11 13:35:17  kirby
* *** empty log message ***
*
* Revision 1.13  2007/07/11 06:35:24  sherwin
* Update after manager reshuffle
*
* Revision 1.12  2007/07/10 21:05:17  kirby
* even more fixes
*
* Revision 1.11  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.10  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.9  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.8  2007/04/18 09:44:01  sherwin
* Working version for StdNodalTri. Removed lapack.cpp from compile.
*
* Revision 1.7  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.6  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.5  2007/01/15 21:13:46  sherwin
* Nodal stuff correction and added Field Classes
*
* Revision 1.4  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.3  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.19  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.18  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.17  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.16  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
**/ 

