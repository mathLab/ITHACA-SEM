///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.cpp
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
// Description: Quadrilateral routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
    namespace StdRegions
    {


        StdQuadExp::StdQuadExp()
        {
        }

        StdQuadExp::StdQuadExp(const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb):
        StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb)
        { 
        }

        StdQuadExp::StdQuadExp(const StdQuadExp &T):
        StdExpansion2D(T)
        {
        }

        StdQuadExp::~StdQuadExp()
        {
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        NekDouble StdQuadExp::Integral(const ConstArray<OneD, NekDouble>& inarray)
        {
            ConstArray<OneD, NekDouble> w0, w1;

            w0 = ExpPointsProperties(0)->GetW();
            w1 = ExpPointsProperties(1)->GetW();

            return StdExpansion2D::Integral(inarray,w0,w1);
        }


        void StdQuadExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                inarray,outarray,1);
        }

        void StdQuadExp::IProductWRTDerivBase(const int dir, 
                                             const ConstArray<OneD, NekDouble>& inarray, 
                                             Array<OneD, NekDouble> & outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),
                                    inarray,outarray,1);
                }
                break;
            case 1:
                {
                    IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),
                                    inarray,outarray,1);  
                }
                break;
            default:
                {
                    ASSERTL1(dir >= 0 &&dir < 2,"input dir is out of range");
                }
                break;
            }             
        }
        

        void StdQuadExp:: IProductWRTBase(const ConstArray<OneD, NekDouble>& base0, 
            const ConstArray<OneD, NekDouble>& base1,
            const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray,
            int coll_check)
        {
            int i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();
            ConstArray<OneD, NekDouble> w0,w1;
            Array<OneD, NekDouble> tmp  = Array<OneD, NekDouble>(nquad0*nquad1);
            Array<OneD, NekDouble> tmp1 = Array<OneD, NekDouble>(nquad0*nquad1);

#if FULLDEBUG
            if((m_base[0]->GetAlpha() != 0.0)||(m_base[1]->GetAlpha() != 0.0))
            {
                ErrorUtil::Error(ErrorUtil::ewarning,"StdQuadExp::IProduct_WRT_B",
                    "Basis has non-zero alpha weight");
            }

            if((m_base[0]->GetBeta() != 0.0)||(m_base[1]->GetBeta() != 0.0))
            {
                ErrorUtil::Error(ErrorUtil::ewarning,"StdQuadExp::IProduct_WRT_B",
                    "Basis has non-zero beta weight");
            }
#endif

            w0 = ExpPointsProperties(0)->GetW();
            w1 = ExpPointsProperties(1)->GetW();
            // Note cannot use outarray as tmp space since dimensions are not always
            // guarenteed to be sufficient 

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,(NekDouble*)&inarray[0]+i*nquad0,1,
                    w0.get(),1,&tmp[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,&tmp[0]+i,nquad0,w1.get(),1,
                    &tmp[0]+i,nquad0);
            }

            if(coll_check&&m_base[0]->Collocation())
            {
                Vmath::Vcopy(order0*nquad1,&tmp[0],1,&tmp1[0],1);
            }
            else
            {
                 Blas::Dgemm('T','N',order0,nquad1,nquad0,1.0,base0.get(),
                             //                    Blas::Dgemm('N','N',order0,nquad1,nquad0,1.0,base0.get(), // make column major matrix
                    nquad0,&tmp[0],nquad0,0.0,&tmp1[0],order0);
            }

            if(coll_check&&m_base[1]->Collocation())
            {
                Vmath::Vcopy(order0*order1,&tmp1[0],1,&outarray[0],1);
            }
            else
            {
                 Blas::Dgemm('N', 'N',order0,order1, nquad1,1.0, &tmp1[0],
                             //                    Blas::Dgemm('N', 'T',order0,order1, nquad1,1.0, &tmp1[0],//make column major matrix
                    order0, base1.get(), nquad1, 0.0, 
                    &outarray[0],order0);
            }

        }

        void StdQuadExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD, NekDouble> base0  = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> base1  = m_base[1]->GetBdata();
            int   btmp0 = m_base[0]->GetNumModes();
            int   mode0 = mode%btmp0;
            int   mode1 = mode/btmp0;


            ASSERTL2(mode1 == (int)floor((1.0*mode)/btmp0),
                "Integer Truncation not Equiv to Floor");

            ASSERTL2(m_ncoeffs <= mode, 
                "calling argument mode is larger than total expansion order");

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vcopy(nquad0,(NekDouble *)(base0.get() + mode0*nquad0),
                    1, &outarray[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode1*nquad1),1,
                    &outarray[0]+i,nquad0,&outarray[0]+i,nquad0);
            }
        }

        DNekMatSharedPtr StdQuadExp::GenMatrix(const StdMatrixKey &mkey)
        {
            int      i;
            int      order0    = GetBasisNumModes(0);
            int      order1    = GetBasisNumModes(1);
            MatrixType mtype   = mkey.GetMatrixType();
            
            DNekMatSharedPtr Mat = StdExpansion::CreateGeneralMatrix(mkey);

            switch(mtype)
            {
            case eMass:
                // For Fourier basis set the imaginary component of mean mode
                // to have a unit diagonal component in mass matrix 
                if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
                {
                    for(i = 0; i < order1; ++i)
                    {
                        (*Mat)(order0*i+1,i*order0+1) = 1.0;
                    }
                }
                
                if(m_base[1]->GetBasisType() == LibUtilities::eFourier)
                {
                    for(i = 0; i < order0; ++i)
                    {
                        (*Mat)(order0+i ,order0+i) = 1.0;
                    }
                }
            }

            return Mat;
        }


        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        void StdQuadExp::PhysDeriv(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &out_d0,
            Array<OneD, NekDouble> &out_d1,
            Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(inarray, out_d0, out_d1);
        }

        //------------------------------
        // Evaluation Methods
        //-----------------------------

        void StdQuadExp::BwdTrans(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            int           nquad0 = m_base[0]->GetNumPoints();
            int           nquad1 = m_base[1]->GetNumPoints();
            int           order0 = m_base[0]->GetNumModes();
            int           order1 = m_base[1]->GetNumModes();
            ConstArray<OneD, NekDouble> base0 = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> base1 = m_base[1]->GetBdata();
            Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(nquad0*std::max(order1,nquad1));

            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(nquad0*order1,&inarray[0],1,&tmp[0],1);
            }
            else
            {
                 Blas::Dgemm('N','N', nquad0,order1,order0,1.0, base0.get(),
                             //                    Blas::Dgemm('N','T', nquad0,order1,order0,1.0, base0.get(), // make column major matrix 
                    nquad0, &inarray[0], order0,0.0,&tmp[0], nquad0);
            }

            if(m_base[1]->Collocation())
            {
                Vmath::Vcopy(nquad0*nquad1,&tmp[0],1,&outarray[0],1);
            }
            else
            {
                Blas::Dgemm('N','T', nquad0, nquad1,order1, 1.0, &tmp[0],
                            //                    Blas::Dgemm('N','N', nquad0, nquad1,order1, 1.0, &tmp[0],//make column major matrix
                    nquad0, base1.get(), nquad1, 0.0, &outarray[0], 
                    nquad0);
            }    
        }

        void StdQuadExp::FwdTrans(const ConstArray<OneD, NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&outarray[0],1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
                DNekMatSharedPtr matsys = GetStdMatrix(masskey);

                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }

        /// Single Point Evaluation
        NekDouble StdQuadExp::PhysEvaluate(ConstArray<OneD, NekDouble>& coords)
        {
            return  StdExpansion2D::PhysEvaluate2D(coords); 
        }

        const ConstArray<OneD, int> StdQuadExp::GetBoundaryMap(void)
        {
            int i;
            int cnt=0;
            int nummodes0, nummodes1;
            int value1, value2;
            Array<OneD, int> returnval(NumBndryCoeffs());

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
                value1 = nummodes0;
                break;
            case LibUtilities::eModified_A:
                value1 = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < value1; i++)
            {
                returnval[i]=i;
            }
            cnt=value1;

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
                value2 = value1+nummodes0-1;
                break;
            case LibUtilities::eModified_A:
                value2 = value1+1;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                returnval[cnt++]=value1+i*nummodes0;
                returnval[cnt++]=value2+i*nummodes0;
            }


            if(Btype1 == LibUtilities::eGLL_Lagrange)
            {
                for(i = nummodes0*(nummodes1-1);i < GetNcoeffs(); i++)
                { 
                    returnval[cnt++] = i;
                }
            }

            return returnval;

        }
        const ConstArray<OneD, int> StdQuadExp::GetInteriorMap(void)
        {
            int i,j;
            int cnt=0;
            int nummodes0, nummodes1;
            int startvalue;
            Array<OneD, int> returnval(GetNcoeffs()-NumBndryCoeffs());

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue = nummodes0;
                break;
            case LibUtilities::eModified_A:
                startvalue = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue++;
                break;
            case LibUtilities::eModified_A:
                startvalue+=2;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                for(j = 0; j < nummodes0-2; j++)
                {
                    returnval[cnt++]=startvalue+j;
                }
                startvalue+=nummodes0;
            }

            return returnval;
        }

        // For a specified edge 'eid' this function updates a class
        // StdExpMap which contains the mapping of the edge degrees of
        // freedom back into the elemental domain which is also
        // dependent upon the edge orientation. The vertex and edge
        // ordering of the mapping is dependent upon which basis is
        // being considered, i.e. modal expansions the vertices will
        // be first, nodal expansions the vertices will be the two
        // end points
        void StdQuadExp::MapTo(const int edge_ncoeffs, 
            const LibUtilities::BasisType Btype, 
            const int eid, 
            const EdgeOrientation eorient,
            StdExpMap &Map)
        {

            int i, start, skip;
            int *dir, order0,order1;
            Array<OneD, int> wsp; 

            ASSERTL2(eid>=0&&eid <=3,"eid must be between 0 and 3");
            // make sure have correct memory storage
            if(edge_ncoeffs != Map.GetLen())
            {
                Map.SetMapMemory(edge_ncoeffs);
            }

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
                if(Btype == LibUtilities::eGLL_Lagrange)
                {
                    for(i = 0; i < edge_ncoeffs; ++i)
                    {
                        dir[i] = edge_ncoeffs-i-1;
                    }
                }
                else{
                    dir[1] = 0; 
                    dir[0] = 1;
                    for(i = 2; i < edge_ncoeffs; ++i)
                    {
                        dir[i] = i;
                    }
                }
            }

            order0 = m_base[0]->GetNumModes();
            order1 = m_base[1]->GetNumModes();

            // Set up Mapping details
            if((eid == 0)||(eid == 2))
            { 
                ASSERTL2(Btype == m_base[0]->GetBasisType(),
                    "Expansion type of edge and StdQuadExp are different");

                switch(Btype)
                {
                case LibUtilities::eGLL_Lagrange:
                    ASSERTL2(edge_ncoeffs == order0,
                        "Expansion order of edge and StdQuadExp are different");

                    if(eid == 0)
                    {
                        start = 0;
                        skip  = 1;
                    }
                    else
                    {
                        start = order0*(order1-1);
                        skip = 1;
                    }
                    break;

                case LibUtilities::eModified_A:
                    if(eid == 0)
                    {
                        start = 0;
                        skip  = 1;
                    }
                    else
                    {
                        start = order0;
                        skip = 1;
                    }
                    break;
                default:
                    ASSERTL0(0,"Mapping array is not defined for this expansion");
                    break;
                }
            }
            else
            {
                ASSERTL2(Btype == m_base[1]->GetBasisType(),
                    "Expansion type of edge and StdQuadExp are different");      

                switch(Btype)
                {
                case LibUtilities::eGLL_Lagrange:
                    ASSERTL2(edge_ncoeffs == order1,
                        "Expansion order of edge and StdQuadExp are different");
                    if(eid == 1)
                    {
                        start = order0-1;
                        skip  = order0;
                    }
                    else
                    {
                        start = 0;
                        skip = order0;
                    }
                    break;

                case LibUtilities::eModified_A:    
                    if(eid == 1)
                    {
                        start = 1;
                        skip  = order0;
                    }
                    else
                    {
                        start = 0;
                        skip = order0;
                    }
                    break;
                default:
                    ASSERTL0(0,"Mapping array is not defined for this expansion");
                    break;
                }
            }

            for(i = 0; i < edge_ncoeffs; ++i)
            {
                Map[dir[i]] = start + i*skip; 
            }

        }

        // same as MapTo but assume that mapping is provided in modal
        // basis type format where vertices are listed first followed
        // by edges degrees of freedom

        void StdQuadExp::MapTo_ModalFormat(const int edge_ncoeffs, 
            const LibUtilities::BasisType Btype, 
            const int eid, 
            const EdgeOrientation eorient,
            StdExpMap &Map)
        {
            MapTo(edge_ncoeffs,Btype,eid,eorient,Map);

            if(Btype == LibUtilities::eGLL_Lagrange)
            {
                int i;
                int vert = Map[edge_ncoeffs-1];
                for(i = edge_ncoeffs-1; i > 1; --i)
                {
                    Map.SetMap(i,Map[i-1]);
                }
                Map.SetMap(1,vert);
            }
        }

        void StdQuadExp::GetCoords(Array<OneD, NekDouble> &coords_0, 
            Array<OneD, NekDouble> &coords_1)
        {
            ConstArray<OneD, NekDouble> z0 = ExpPointsProperties(0)->GetZ();
            ConstArray<OneD, NekDouble> z1 = ExpPointsProperties(1)->GetZ();
            int nq0 = GetNumPoints(0);
            int nq1 = GetNumPoints(1);
            int i;

            for(i = 0; i < nq1; ++i)
            {
                Blas::Dcopy(nq0,z0.get(), 1,&coords_0[0] + i*nq0,1);
                Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
            }
        }
    } //end namespace            
}//end namespace

/** 
* $Log: StdQuadExp.cpp,v $
* Revision 1.25  2007/12/06 22:44:47  pvos
* 2D Helmholtz solver updates
*
* Revision 1.24  2007/11/08 16:55:14  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.23  2007/10/15 20:38:54  ehan
* Make changes of column major matrix
*
* Revision 1.22  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.21  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.20  2007/07/11 13:35:18  kirby
* *** empty log message ***
*
* Revision 1.19  2007/07/10 21:05:17  kirby
* even more fixes
*
* Revision 1.17  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.16  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.15  2007/05/15 05:18:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.14  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.13  2007/04/06 08:44:43  sherwin
* Update to make 2D regions work at StdRegions level
*
* Revision 1.12  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.11  2007/04/05 11:39:47  pvincent
* quad_edited
*
* Revision 1.10  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.9  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.8  2007/01/18 18:44:45  bnelson
* Updates to compile on Visual Studio 2005.
*
* Revision 1.7  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.6  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.2  2006/06/01 14:13:36  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.38  2006/04/25 20:23:34  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.37  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.36  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.35  2006/03/05 22:11:03  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.34  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.33  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/ 





