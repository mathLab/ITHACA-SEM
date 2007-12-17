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

        StdNodalTriExp::StdNodalTriExp(const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb, 
            LibUtilities::PointsType Ntype):
        StdTriExp(Ba,Bb)
        {

            ASSERTL0(m_base[0]->GetNumModes() == m_base[1]->GetNumModes(),
                "Nodal basis initiated with different orders in the a "
                "and b directions");

            int numModes = Ba.GetNumModes();
            m_nodalPointsKey = MemoryManager<LibUtilities::PointsKey>::AllocateSharedPtr (numModes,Ntype);

        }


        StdNodalTriExp::StdNodalTriExp(const StdNodalTriExp &T):
        StdTriExp(T)
        {
            m_nodalPointsKey = T.m_nodalPointsKey;
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
            Array<OneD, NekDouble> c = Array<OneD, NekDouble>(2);
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
                    // define matrix in row major format to have rows
                    // of all the different expansion bases defined at
                    // the nodal point
                    (*Mat)(j,i) = StdTriExp::PhysEvaluate(c);
                }
            }
            return Mat;
        }

        void StdNodalTriExp::NodalToModal()
        {
            NodalToModal(m_coeffs); 
        }

        void StdNodalTriExp::NodalToModal(Array<OneD, NekDouble> &in_out_array)
        {
            ASSERTL0(false,"This function needs fixing");
            //StdMatrixKey   Nkey(eInvNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            //DNekMatSharedPtr  matsys = m_stdMatManager[Nkey];

            // solve inverse of system
            //DNekVec   v(m_ncoeffs,in_out_array,eWrapper);
        }


        void StdNodalTriExp::NodalToModalTranspose()
        {
            NodalToModalTranspose(m_coeffs); 
        }

        // Operate with transpose of NodalToModal transformation
        void StdNodalTriExp::NodalToModalTranspose(Array<OneD, NekDouble> &in_out_array)
        {
            StdMatrixKey Nkey(eInvNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  mat = GetStdMatrix(Nkey);

            // solve inverse of system
            DNekVec   v(m_ncoeffs,in_out_array,eWrapper);
            v = (*mat)*v;
        }


        void StdNodalTriExp::ModalToNodal()
        {
            ModalToNodal(m_coeffs);
        }

        void StdNodalTriExp::ModalToNodal(Array<OneD, NekDouble> &in_out_array)
        {
            StdMatrixKey      Nkey(eNBasisTrans,DetShapeType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  mat = GetStdMatrix(Nkey);

            // Multiply out matrix
            DNekVec  v(m_ncoeffs,in_out_array,eWrapper);
            v = (*mat)*v;

            // this line should not be needed
            Vmath::Vcopy(m_ncoeffs,&v[0],1,&in_out_array[0],1);
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

            NodalToModalTranspose(outarray);
        }

        void StdNodalTriExp::FillMode(const int mode, 
            Array<OneD, NekDouble> &outarray)
        {

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs,&outarray[0],1);

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
            Array<OneD, NekDouble> tmp  = Array<OneD, NekDouble>(m_ncoeffs);

            // save nodal values
            Blas::Dcopy(m_ncoeffs,&inarray[0],1,&tmp[0],1);
            NodalToModal(tmp);
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
            DNekVec in (m_ncoeffs,outarray);
            DNekVec out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
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

