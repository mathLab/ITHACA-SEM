///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion.cpp
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
// Description: Definition of methods in class StdExpansion which is
// the base class to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////


#include <StdRegions/StdExpansion.h>
#include <LibUtilities/BasicUtils/SharedArrayUtil.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        /** define list of number of vertices corresponding to each ShapeType */
        const int g_shapenverts[SIZE_ShapeType] = {0,2,3,4,4,5,6,8};

        /** define list of number of edges corresponding to each ShapeType */
        const int g_shapenedges[SIZE_ShapeType] = {0,1,3,4,6,8,9,12};

        /** define list of number of faces corresponding to each ShapeType */
        const int g_shapenfaces[SIZE_ShapeType] = {0,0,0,0,4,5,5,6};


        StdExpansion::StdExpansion(void): 
        m_ncoeffs(0),
            m_numbases(0)
        {
        }

        StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
            const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb, 
            const LibUtilities::BasisKey &Bc):
        m_ncoeffs(numcoeffs),
            m_numbases(numbases)
        {

            m_base = Array<OneD, LibUtilities::BasisSharedPtr>(m_numbases);

            switch(m_numbases)
            {
            case 3:
                ASSERTL2(Bc!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[2] = LibUtilities::BasisManager()[Bc];

            case 2:
                ASSERTL2(Bb!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");

                m_base[1] = LibUtilities::BasisManager()[Bb];
            case 1:
                ASSERTL2(Ba!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[0] = LibUtilities::BasisManager()[Ba];
                break;
            default:
                ASSERTL0(false, "numbases incorrectly specified");
            };

            //allocate memory for coeffs
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());

            // Register Creators for  Managers
            for(int i = 0; i < SIZE_MatrixType; ++i)
            {
                m_stdMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoShapeType,*this),
                                  boost::bind(&StdExpansion::CreateStdMatrix, this, _1));
            }

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T)
        {
            CopyArray(T.m_base, m_base);

            // NOTE: Copy Constructor produces a deep copy
            // allocate memory for coeffs
            // need to check allocation for variable order. 
            m_ncoeffs = T.m_ncoeffs;
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            Vmath::Vcopy(m_ncoeffs,&T.m_coeffs[0],1,&m_coeffs[0],1);

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());
            Vmath::Vcopy(GetTotPoints(),&T.m_phys[0],1,&m_phys[0],1);
        }

        StdExpansion::~StdExpansion()
        {
        }


        DNekMatSharedPtr StdExpansion::CreateStdMatrix(const StdMatrixKey &mkey) 
        {
            DNekMatSharedPtr returnval;
            MatrixType mtype = mkey.GetMatrixType();

            switch(mtype)
            {
            case eMass: 
            case eLaplacian:
            case eLaplacian00:
            case eLaplacian01:
            case eLaplacian02:
            case eLaplacian11:
            case eLaplacian12:
            case eLaplacian22:
            case eBwdTrans:
            case eNBasisTrans:
                returnval = StdExpansion::GenMatrix(mtype);
                break;

            case eInvMass:
                {
                    StdMatrixKey masskey(eMass,mkey.GetShapeType(),mkey.GetBase(),
                        mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr mmat = m_stdMatrixManager[masskey];
                    DNekLinSys invmass(mmat);

                    int dim = mmat->GetRows(); //assume square matrix 
                    Array<OneD,NekDouble> invdata = Array<OneD,NekDouble>(dim*dim);
                    Array<OneD,NekDouble> data_offset;

                    Vmath::Zero(dim*dim,&invdata[0],1);

                    for(int i = 0; i < dim; ++i)
                    {
                        // set array to be identity matrix
                        invdata[i*dim + i] = 1.0;

                        //call inverse on symmetric matrix
                        DNekVec   v(dim,data_offset = invdata+i*dim,eWrapper);
                        invmass.Solve(v,v);
                    }

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(dim,dim,&invdata[0]);
                }
                break;
                            
            case eInvNBasisTrans:
                {
                    StdMatrixKey tmpkey(eNBasisTrans,mkey.GetShapeType(),mkey.GetBase(),
                        mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr tmpmat = m_stdMatrixManager[tmpkey];
                    DNekLinSys invtmpmat(tmpmat);

                    int dim = tmpmat->GetRows(); //assume square matrix 
                    Array<OneD,NekDouble> invdata = Array<OneD,NekDouble>(dim*dim);
                    Array<OneD,NekDouble> data_offset;

                    Vmath::Zero(dim*dim,&invdata[0],1);

                    for(int i = 0; i < dim; ++i)
                    {
                        // set array to be identity matrix
                        invdata[i*dim + i] = 1.0;

                        //call inverse on symmetric matrix
                        DNekVec   v(dim,data_offset = invdata+i*dim,eWrapper);
                        invtmpmat.Solve(v,v);
                    }

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(dim,dim,&invdata[0]);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }


        NekDouble StdExpansion::Linf(const ConstArray<OneD, NekDouble>& sol)
        {
            int     ntot;
            NekDouble  val;
            Array<OneD, NekDouble>  wsp;

            ntot = GetTotPoints();
            wsp  = Array<OneD, NekDouble>(ntot);

            Vmath::Vsub(ntot,sol.get(),1,&m_phys[0],1,&wsp[0],1);
            Vmath::Vabs(ntot,&wsp[0],1,&wsp[0],1);
            val = Vmath::Vamax(ntot,&wsp[0],1);    

            return  val;
        }

        NekDouble StdExpansion::Linf()
        {
            return Vmath::Vamax(GetTotPoints(),&m_phys[0],1);    
        }

        NekDouble StdExpansion::L2(const ConstArray<OneD, NekDouble>& sol)
        {
            int     ntot = GetTotPoints();
            NekDouble  val;
            Array<OneD, NekDouble> wsp;

            wsp = Array<OneD, NekDouble>(ntot);

            Vmath::Vsub(ntot, sol.get(), 1, m_phys.get(), 1, wsp.get(), 1);
            Vmath::Vmul(ntot, wsp.get(), 1, wsp.get(),  1, wsp.get(), 1);

            val = sqrt(v_Integral(wsp));

            return val;
        }

        NekDouble StdExpansion::L2()
        {
            int     ntot = GetTotPoints();
            NekDouble  val;
            Array<OneD, NekDouble> wsp;

            wsp = Array<OneD, NekDouble>(ntot);

            Vmath::Vmul(ntot, &m_phys[0], 1, &m_phys[0], 1, &wsp[0], 1);
            val   = sqrt(v_Integral(wsp));

            return val;
        }

        

        DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(MatrixType mtype)
        {
            int     i;
            DNekMatSharedPtr  Mat;

            Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(m_ncoeffs);

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);            
            for(i=0; i<m_ncoeffs; ++i)
            {
                Vmath::Zero(m_ncoeffs,&tmp[0],1);
                tmp[i] = 1.0;
                
                GeneralMatrixOp(mtype,tmp,tmp);

                Vmath::Vcopy(m_ncoeffs,&tmp[0],1,&((*Mat).GetPtr())[0]+i*m_ncoeffs,1);
            }

            return Mat;
        }

       
        void StdExpansion::GeneralMatrixOp(MatrixType mtype, 
                                           const ConstArray<OneD,NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           double lambda)
        {
            switch(mtype)
            {
            case eMass:
                MassMatrixOp(inarray,outarray);
                break;
            case eLaplacian:
                LaplacianMatrixOp(inarray,outarray);
                break;
            case eLaplacian00:
                LaplacianMatrixOp(0,0,inarray,outarray);
                break;
            case eLaplacian11:
                LaplacianMatrixOp(1,1,inarray,outarray);
                break;
            case eLaplacian22:
                LaplacianMatrixOp(2,2,inarray,outarray);
                break;
            case eBwdTrans:
                BwdTransMatrixOp(inarray,outarray);
                break;
            case eHelmholtz:
                HelmholtzMatrixOp(inarray,outarray,lambda);
                break;
            case eNBasisTrans:
                NEKERROR(ErrorUtil::efatal,"This is a specialised matrix for nodal "
                         "expansions only and does not have an operator");
                break;

            default:
                NEKERROR(ErrorUtil::efatal, "Need to populate switch");
                break;
            }
        }
            
        void StdExpansion::MassMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp   = Array<OneD, NekDouble>(GetTotPoints());

            v_BwdTrans(inarray,tmp);
            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::LaplacianMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {

            switch(ShapeTypeDimMap[v_DetShapeType()])
            {
            case 1:
                LaplacianMatrixOp(0,0,inarray,outarray);
                break;

            case 2:
                {
                    Array<OneD, NekDouble> store = Array<OneD, NekDouble>(m_ncoeffs);
                    
                    LaplacianMatrixOp(0,0,inarray,store);
                    LaplacianMatrixOp(1,1,inarray,outarray);
                    
                    Vmath::Vadd(m_ncoeffs,&store[0],1,&outarray[0],1,&outarray[0],1);
                }
                break;
            case 3:
                {
                    Array<OneD, NekDouble> store0 = Array<OneD, NekDouble>(m_ncoeffs);
                    Array<OneD, NekDouble> store1 = Array<OneD, NekDouble>(m_ncoeffs);
                    
                    LaplacianMatrixOp(0,0,inarray,store0);
                    LaplacianMatrixOp(1,1,inarray,store1);
                    LaplacianMatrixOp(2,2,inarray,outarray);
                    
                    Vmath::Vadd(m_ncoeffs,&store0[0],1,&outarray[0],1,&outarray[0],1);
                    Vmath::Vadd(m_ncoeffs,&store1[0],1,&outarray[0],1,&outarray[0],1);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Dimenion not recognised.");
                break;
            }
        }


        void StdExpansion::LaplacianMatrixOp(const int k1, const int k2, 
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
                
            ASSERTL1(k1 >= 0 && k1 < ShapeTypeDimMap[v_DetShapeType()],"invalid first  argument");
            ASSERTL1(k2 >= 0 && k2 < ShapeTypeDimMap[v_DetShapeType()],"invalid second argument");
                                  
            Array<OneD, NekDouble> tmp   = Array<OneD, NekDouble>(GetTotPoints());
            Array<OneD, NekDouble> dtmp  = Array<OneD, NekDouble>(GetTotPoints());
            
            v_BwdTrans(inarray,tmp);
            v_StdPhysDeriv(k2,tmp,dtmp);
            v_IProductWRTDerivBase(k1,dtmp, outarray);
        }


        void StdExpansion::WeakDerivMatrixOp(const int i,
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
            NEKERROR(ErrorUtil::efatal, "Method not populated yet.");
        }


        void StdExpansion::BwdTransMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
            v_BwdTrans(inarray,outarray);
        }


        void StdExpansion::HelmholtzMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const double lambda)
        {
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(m_ncoeffs);
            LaplacianMatrixOp(inarray,tmp);
            MassMatrixOp(inarray,outarray);

            Blas::Daxpy(m_ncoeffs,lambda,&tmp[0],1,&outarray[0],1);
        }



        void StdExpansion::TensProdBwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      bwdtransmatkey(eBwdTrans,DetShapeType(),*this);
            DNekMatSharedPtr  bwdtransmat = m_stdMatrixManager[bwdtransmatkey];

            DNekVec v_in(m_ncoeffs,inarray);
            DNekVec v_out(nq,outarray,eWrapper);

            v_out = (*bwdtransmat) * v_in;
            // This line below should be removed once the eWrapper method of NekVEctor works properly
            Vmath::Vcopy(nq,&((v_out).GetPtr())[0],1,&outarray[0],1);
        }

        // 2D Interpolation
        void StdExpansion::Interp2D(const LibUtilities::BasisKey &fbasis0, 
            const LibUtilities::BasisKey &fbasis1, 
            const ConstArray<OneD, NekDouble>& from,  
            const LibUtilities::BasisKey &tbasis0,
            const LibUtilities::BasisKey &tbasis1,
            Array<OneD, NekDouble> &to)
        {
            Interp2D(fbasis0,fbasis1,from.data(),tbasis0,tbasis1,to.data());
        }

        void StdExpansion::Interp2D(const LibUtilities::BasisKey &fbasis0, 
            const LibUtilities::BasisKey &fbasis1, 
            const NekDouble *from,  
            const LibUtilities::BasisKey &tbasis0,
            const LibUtilities::BasisKey &tbasis1,
            NekDouble *to)
        {
            DNekMatSharedPtr I0,I1;
            Array<OneD, NekDouble> wsp = Array<OneD, NekDouble>(tbasis1.GetNumPoints()*
                fbasis0.GetNumPoints());

            I0 = LibUtilities::PointsManager()[fbasis0.GetPointsKey()]
            ->GetI(tbasis0.GetPointsKey());
            I1 = LibUtilities::PointsManager()[fbasis1.GetPointsKey()]
            ->GetI(tbasis1.GetPointsKey());

            Blas::Dgemm('T', 'T', tbasis1.GetNumPoints(), fbasis0.GetNumPoints(),
                fbasis1.GetNumPoints(), 1.0, I1->GetPtr().get(),  
                fbasis1.GetNumPoints(),
                from,fbasis0.GetNumPoints(), 0.0, wsp.get(),
                tbasis1.GetNumPoints());

            Blas::Dgemm('T', 'T',tbasis0.GetNumPoints(),tbasis1.GetNumPoints(),
                fbasis0.GetNumPoints(),1.0,I0->GetPtr().get(),
                fbasis0.GetNumPoints(),wsp.get(), tbasis1.GetNumPoints(),
                0.0,to, tbasis0.GetNumPoints());
        }

        // 1D Interpolation
        void StdExpansion::Interp1D(const LibUtilities::BasisKey &fbasis0, 
            const ConstArray<OneD, NekDouble>& from,  
            const LibUtilities::BasisKey &tbasis0, 
            Array<OneD, NekDouble> &to)
        {
            Interp1D(fbasis0, from.get(), tbasis0, to.get());
        }

        void StdExpansion::Interp1D(const LibUtilities::BasisKey &fbasis0, 
            const NekDouble *from,  
            const LibUtilities::BasisKey &tbasis0, 
            NekDouble *to)
        {
            DNekMatSharedPtr I0;

            I0 = LibUtilities::PointsManager()[fbasis0.GetPointsKey()]
            ->GetI(tbasis0.GetPointsKey());

            //DNekVec in(fbasis0.GetNumPoints(),from);
            //DNekVec out(tbasis0.GetNumPoints(),to,eWrapper);

            //out  = (*I0)*in;
            // this line should not be needed
            //Vmath::Vcopy(tbasis0.GetNumPoints(),&out[0],1,&to[0],1);

            Blas::Dgemv('T', fbasis0.GetNumPoints(), tbasis0.GetNumPoints(), 
                1.0, I0->GetPtr().get(), fbasis0.GetNumPoints(), 
                from, 1, 0.0, to, 1);
        }

        //   I/O routine
        void StdExpansion::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int i;
            for(i=0; i<m_ncoeffs; ++i)
            {
                outfile << m_coeffs[i] << std::endl;
            }
        }

    }//end namespace
}//end namespace

/**
* $Log: StdExpansion.cpp,v $
* Revision 1.42  2007/07/12 12:55:14  sherwin
* Simplified Matrix Generation
*
* Revision 1.41  2007/07/10 20:41:52  kirby
* more fixes
*
* Revision 1.40  2007/07/10 19:27:58  kirby
* Update for new matrix structures
*
* Revision 1.39  2007/07/09 15:19:14  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.38  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.37  2007/05/30 23:56:54  sherwin
* Silly errors
*
* Revision 1.36  2007/05/30 20:49:12  sherwin
* Updates to do with LocalRegions and SpatialDomains
*
* Revision 1.35  2007/05/23 15:12:45  pvos
* removed some obsolete lines
*
* Revision 1.34  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.33  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.32  2007/04/18 16:09:12  pvos
* Added some new Tensor Operations routines
*
* Revision 1.31  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.30  2007/04/08 03:36:57  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.29  2007/03/31 00:40:02  bnelson
* *** empty log message ***
*
* Revision 1.28  2007/03/29 19:35:08  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.27  2007/03/25 15:48:22  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.26  2007/03/21 20:56:42  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.25  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.24  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.23  2007/03/02 12:01:51  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.22  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.21  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.20  2007/02/23 19:26:07  jfrazier
* General bug fix and formatting.
*
* Revision 1.19  2007/02/22 22:02:27  sherwin
* Update with executing StdMatManager
*
* Revision 1.18  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.17  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.16  2007/02/17 04:03:22  jfrazier
* Added NekManager for holding matrices.  Need to finish the create function.
*
* Revision 1.15  2007/02/14 16:35:50  pvos
* Corrected an error in the code
*
* Revision 1.14  2007/02/13 09:52:27  sherwin
* Updates to fix mass matrix inverse issues
*
* Revision 1.13  2007/02/07 12:51:52  sherwin
* Compiling version of Project1D
*
* Revision 1.12  2007/02/06 02:23:28  jfrazier
* Minor cleanup.
*
* Revision 1.11  2007/01/30 20:01:35  sherwin
* Update for first compiling Project1D routine
*
* Revision 1.10  2007/01/29 15:04:53  sherwin
* StdBasis.h moved to LibUtilities. Other minor mods
*
* Revision 1.9  2007/01/28 18:34:18  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.8  2007/01/23 23:20:20  sherwin
* New version after Jan 07 update
*
* Revision 1.7  2007/01/20 22:35:20  sherwin
* Version with StdExpansion compiling
*
* Revision 1.6  2007/01/15 11:08:37  pvos
* Updating doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.2  2006/05/29 19:03:08  sherwin
* Modifications to wrap geometric information in shared_ptr
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.54  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.53  2006/04/01 21:59:26  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.52  2006/03/21 09:21:31  sherwin
* Introduced NekMemoryManager
*
* Revision 1.51  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.50  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.49  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.48  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.47  2006/03/02 16:20:20  sherwin
*
* Introduced method GetPointsTot
*
* Revision 1.46  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.45  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.44  2006/02/26 21:23:20  bnelson
* Fixed a variety of compiler errors caused by updates to the coding standard.
*
* Revision 1.43  2006/02/15 08:06:36  sherwin
*
* Put files into coding standard (although they do not compile)
*
* Revision 1.42  2006/02/12 21:51:42  sherwin
*
* Added licence
*
* Revision 1.41  2006/02/10 16:44:10  sherwin
*
* Updated to comply with coding standard
*
**/
