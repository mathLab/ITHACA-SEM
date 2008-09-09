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

namespace Nektar
{
    namespace StdRegions
    {

        /** define list of number of vertices corresponding to each ShapeType */
        const int g_shapenverts[SIZE_ExpansionType] = {0,2,3,4,4,5,6,8};

        /** define list of number of edges corresponding to each ShapeType */
        const int g_shapenedges[SIZE_ExpansionType] = {0,1,3,4,6,8,9,12};

        /** define list of number of faces corresponding to each ShapeType */
        const int g_shapenfaces[SIZE_ExpansionType] = {0,0,0,0,4,5,5,6};

        StdExpansion::StdExpansion(void): 
            m_elmt_id(0),
            m_numbases(0),
            m_ncoeffs(0)           
        {
        }

        StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
            const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb, 
            const LibUtilities::BasisKey &Bc):
            m_elmt_id(0),
            m_numbases(numbases),
            m_base(m_numbases),
            m_ncoeffs(numcoeffs),
            m_coeffs(m_ncoeffs,0.0),
            m_stdMatrixManager(std::string("StdExpansionStdMatrix")),
            m_stdStaticCondMatrixManager(std::string("StdExpansionStdStaticCondMatrix"))
        {
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

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());

            // Register Creators for  Managers
            for(int i = 0; i < SIZE_MatrixType; ++i)
            {
                m_stdMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoExpansionType,*this),
                                   boost::bind(&StdExpansion::CreateStdMatrix, this, _1));
                m_stdStaticCondMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoExpansionType,*this),  
                                   boost::bind(&StdExpansion::CreateStdStaticCondMatrix, this, _1));
            }

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T):
            m_elmt_id(T.m_elmt_id),
            m_numbases(T.m_numbases),
            m_base(T.m_base),
            m_ncoeffs(T.m_ncoeffs),
            m_coeffs(m_ncoeffs),
            m_phys((T.m_phys).num_elements()),
            m_stdMatrixManager(std::string("StdExpansion")),
            m_stdStaticCondMatrixManager(std::string("StdExpansionStaticCondMat"))
        {
            //CopyArray(T.m_base, m_base); 
            CopyArray(T.m_coeffs, m_coeffs);
            CopyArray(T.m_phys, m_phys);

            // Register Creators for  Managers
            for(int i = 0; i < SIZE_MatrixType; ++i)
            {
                m_stdMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoExpansionType,*this), boost::bind(&StdExpansion::CreateStdMatrix, this, _1));
                m_stdStaticCondMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoExpansionType,*this), boost::bind(&StdExpansion::CreateStdStaticCondMatrix, this, _1));
            }
        }
        
        StdExpansion::~StdExpansion()
        {
        }

        NekDouble StdExpansion::Linf(const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble>  wsp(ntot);

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, wsp, 1);
            Vmath::Vabs(ntot, wsp, 1, wsp, 1);
            val = Vmath::Vamax(ntot, wsp, 1);

            return  val;
        }

        DNekBlkMatSharedPtr StdExpansion::CreateStdStaticCondMatrix(const StdMatrixKey &mkey) 
        {
            DNekBlkMatSharedPtr returnval;
            MatrixType mtype = mkey.GetMatrixType();
            
            DNekMatSharedPtr&  mat = GetStdMatrix(mkey);
            int nbdry = NumBndryCoeffs(); // also checks to see if this is a boundary interior decomposed expansion
            int nint = m_ncoeffs - nbdry;
            DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
            DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
            DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
            DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);
            
            int i,j;

            Array<OneD,unsigned int> bmap(nbdry);
            Array<OneD,unsigned int> imap(nint);
            GetBoundaryMap(bmap);
            GetInteriorMap(imap);

            for(i = 0; i < nbdry; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*A)(i,j) = (*mat)(bmap[i],bmap[j]);
                }
                
                for(j = 0; j < nint; ++j)
                {
                    (*B)(i,j) = (*mat)(bmap[i],imap[j]);
                }
            }
            
            for(i = 0; i < nint; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*C)(i,j) = (*mat)(imap[i],bmap[j]);
                }
                
                for(j = 0; j < nint; ++j)
                {
                    (*D)(i,j) = (*mat)(imap[i],imap[j]);
                }
            }
            
            // Calculate static condensed system 
            if(nint)
            {
                D->Invert();
                (*B) = (*B)*(*D);
                (*A) = (*A) - (*B)*(*C);
            }

            // set up block matrix system
            Array<OneD, unsigned int> exp_size(2);
            exp_size[0] = nbdry;
            exp_size[1] = nint;
            returnval = MemoryManager<DNekBlkMat>::AllocateSharedPtr(exp_size,exp_size); 
            
            returnval->SetBlock(0,0,A);
            returnval->SetBlock(0,1,B);
            returnval->SetBlock(1,0,C);
            returnval->SetBlock(1,1,D);
            
            return returnval;
        }

        NekDouble StdExpansion::Linf()
        {
            return Vmath::Vamax(GetTotPoints(), m_phys, 1);    
        }

        NekDouble StdExpansion::L2(const Array<OneD, const NekDouble>& sol)
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble> wsp(ntot);

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, wsp, 1);
            Vmath::Vmul(ntot, wsp, 1, wsp, 1, wsp, 1);

            val = sqrt(v_Integral(wsp));

            return val;
        }

        NekDouble StdExpansion::L2()
        {
            NekDouble  val;
            int     ntot = GetTotPoints();
            Array<OneD, NekDouble> wsp(ntot);

            Vmath::Vmul(ntot, m_phys, 1, m_phys, 1, wsp, 1);
            val   = sqrt(v_Integral(wsp));

            return val;
        }

        DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(const StdMatrixKey &mkey)
        {
            int     i;
            DNekMatSharedPtr  returnval;

            switch(mkey.GetMatrixType())
            {
            case eInvMass:
                {
                    StdMatrixKey masskey(eMass,mkey.GetExpansionType(),mkey.GetBase(),  mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr& mmat = GetStdMatrix(masskey);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*mmat); //Populate standard mass matrix.
                    returnval->Invert();
                }
                break;
            case eInvNBasisTrans:
                {
                    StdMatrixKey tmpkey(eNBasisTrans,mkey.GetExpansionType(),mkey.GetBase(),
                        mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr& tmpmat = GetStdMatrix(tmpkey);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*tmpmat); //Populate standard mass matrix.
                    returnval->Invert();
                }
                break;
            case eBwdTrans:
                {
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmp(nq);
                    
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nq,m_ncoeffs);            
                    DNekMat &Mat = *returnval;                     
                    
                    for(int i=0; i<m_ncoeffs; ++i)
                    {                        
                        v_FillMode(i,tmp);                        
                        Vmath::Vcopy(nq,&tmp[0],1,&(Mat.GetPtr())[0]+i*nq,1);
                    }
                }
                break;                            
            case eMass: 
            case eHelmholtz:
            case eLaplacian:
            case eLaplacian00:
            case eLaplacian01:
            case eLaplacian02:
            case eLaplacian11:
            case eLaplacian12:
            case eLaplacian22:
            case eWeakDeriv0:
            case eWeakDeriv1:
            case eWeakDeriv2:
                {
                    Array<OneD, NekDouble> tmp(m_ncoeffs);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);            
                    DNekMat &Mat = *returnval; 
                    
                    for(i=0; i < m_ncoeffs; ++i)
                    {
                        Vmath::Zero(m_ncoeffs, tmp, 1);
                        tmp[i] = 1.0;
                        
                        GeneralMatrixOp(mkey,tmp,tmp);
                        
                        Vmath::Vcopy(m_ncoeffs,&tmp[0],1,
                                     &(Mat.GetPtr())[0]+i*m_ncoeffs,1);
                    }
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "This type of matrix can not be created using a general approach");
                }
                break;
            }
            
            return returnval;
        }

        void StdExpansion::GeneralMatrixOp(const StdMatrixKey &mkey, 
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
        {
            switch(mkey.GetMatrixType())
            {
            case eMass:
                MassMatrixOp(inarray,outarray);
                break;
            case eWeakDeriv0:
                WeakDerivMatrixOp(0,inarray,outarray);
                break;
            case eWeakDeriv1:
                WeakDerivMatrixOp(1,inarray,outarray);
                break;
            case eWeakDeriv2:
                WeakDerivMatrixOp(2,inarray,outarray);
                break;
            case eLaplacian:
                LaplacianMatrixOp(inarray,outarray);
                break;
            case eLaplacian00:
                LaplacianMatrixOp(0,0,inarray,outarray);
                break;
            case eLaplacian01:
                LaplacianMatrixOp(0,1,inarray,outarray);
                break;
            case eLaplacian11:
                LaplacianMatrixOp(1,1,inarray,outarray);
                break;
            case eLaplacian22:
                LaplacianMatrixOp(2,2,inarray,outarray);
                break;
            case eBwdTrans:
                BwdTrans(inarray,outarray);
                break;
            case eHelmholtz:
                HelmholtzMatrixOp(inarray,outarray,mkey.GetConstant(0));
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "This matrix does not have an operator");
                break;
            }
        }
            
        void StdExpansion::MassMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(GetTotPoints());

            v_BwdTrans(inarray,tmp);            
            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::LaplacianMatrixOp(const int k1, const int k2, 
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {                
            ASSERTL1(k1 >= 0 && k1 < ExpansionTypeDimMap[v_DetExpansionType()],"invalid first  argument");
            ASSERTL1(k2 >= 0 && k2 < ExpansionTypeDimMap[v_DetExpansionType()],"invalid second argument");
                                  
            Array<OneD, NekDouble> tmp(GetTotPoints());
            Array<OneD, NekDouble> dtmp(GetTotPoints());
            
            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k2,tmp,dtmp);
            v_IProductWRTDerivBase(k1, dtmp, outarray);
        }


        void StdExpansion::WeakDerivMatrixOp(const int k1,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
            ASSERTL1(k1 >= 0 && k1 < ExpansionTypeDimMap[v_DetExpansionType()],"invalid first  argument");
                                  
            Array<OneD, NekDouble> tmp(GetTotPoints());
            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k1,tmp,tmp);
            v_IProductWRTBase(tmp, outarray);
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
* Revision 1.75  2008/08/14 22:09:50  sherwin
* Modifications to remove HDG routines from StdRegions and removed StdExpMap
*
* Revision 1.74  2008/07/29 22:21:15  sherwin
* A bunch of mods for DG advection and separaring the GetGeom calls into GetGeom1D ...
*
* Revision 1.73  2008/07/19 21:12:54  sherwin
* Removed MapTo function and made orientation convention anticlockwise in UDG routines
*
* Revision 1.72  2008/07/12 16:30:07  sherwin
* Added an new member m_elmt_id so that there is an element number for use later in lists
*
* Revision 1.71  2008/07/02 14:08:56  pvos
* Implementation of HelmholtzMatOp and LapMatOp on shape level
*
* Revision 1.70  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.69  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.68  2008/04/06 06:04:14  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.67  2008/04/02 22:18:10  pvos
* Update for 2D local to global mapping
*
* Revision 1.66  2008/03/18 14:15:45  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.65  2008/03/12 15:25:09  pvos
* Clean up of the code
*
* Revision 1.63  2008/02/29 19:15:19  sherwin
* Update for UDG stuff
*
* Revision 1.62  2008/02/16 05:59:14  ehan
* Added interpolation 3D.
*
* Revision 1.61  2008/01/23 09:09:46  sherwin
* Updates for Hybrized DG
*
* Revision 1.60  2007/12/17 13:03:45  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.59  2007/12/06 22:44:46  pvos
* 2D Helmholtz solver updates
*
* Revision 1.58  2007/11/29 21:40:20  sherwin
* updates for MultiRegions and DG solver
*
* Revision 1.57  2007/11/08 16:55:12  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.56  2007/10/15 20:38:32  ehan
* Tested standard mass matrix
*
* Revision 1.55  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.54  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.53  2007/09/27 12:55:57  pvos
* Column major Blas calls corrections
*
* Revision 1.52  2007/09/25 14:25:56  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.51  2007/08/29 23:26:48  jfrazier
* Created non-static manager that shares data across instances.
*
* Revision 1.50  2007/07/27 16:56:50  jfrazier
* Changed manager to static.
*
* Revision 1.49  2007/07/27 00:22:53  bnelson
* Memory manager now accepts non-const parameters to the allocate methods.
*
* Revision 1.48  2007/07/22 23:04:25  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.47  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.46  2007/07/15 19:28:28  bnelson
* *** empty log message ***
*
* Revision 1.45  2007/07/13 15:20:19  kirby
* *** empty log message ***
*
* Revision 1.43  2007/07/13 09:02:25  sherwin
* Mods for Helmholtz solver
*
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
