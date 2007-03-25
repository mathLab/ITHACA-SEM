///////////////////////////////////////////////////////////////////////////////
//
// File StdSegExp.cpp
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
// Description: Routines within Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdSegExp.h>


namespace Nektar
{
    namespace StdRegions
    {

	StdSegExp::StdSegExp(const LibUtilities::BasisKey &Ba):
	    StdExpansion1D(Ba.GetNumModes(),Ba)
	{    
	}
	
	StdSegExp::StdSegExp(const StdSegExp &T):
	    StdExpansion1D(T)
	{
	}
	
	StdSegExp::~StdSegExp()
	{    
	}
	

	//----------------------------
	// Integration Methods
	//----------------------------
    
	double StdSegExp::Integral(const NekDoubleSharedArray &inarray)
	{
	    double Int = 0.0;
	    const NekDouble *z,*w0;
	    int    nquad0 = m_base[0]->GetNumPoints();
	    NekDoubleSharedArray tmp = GetDoubleTmpSpace(nquad0);
    
	    ExpPointsProperties(0)->GetZW(z,w0);

	    // multiply by integration constants 
	    Vmath::Vmul(nquad0,&inarray[0],1,(NekDouble*)w0,1,&tmp[0],1);
    
	    Int = Vmath::Vsum(nquad0,&tmp[0],1);
      
	    return Int;
	}
	
	void StdSegExp::IProductWRTBase(const NekDouble *base, 
					const NekDoubleSharedArray &inarray, 
					NekDoubleSharedArray &outarray, int coll_check)
	{
	    int    nquad = m_base[0]->GetNumPoints();
	    const NekDouble *z,*w;
	    NekDoubleSharedArray tmp  = GetDoubleTmpSpace(nquad);

	    ExpPointsProperties(0)->GetZW(z,w);

	    Vmath::Vmul(nquad,&inarray[0],1,(NekDouble*)w,1,&tmp[0],1);

	    if(coll_check&&m_base[0]->Collocation())
	    {
		Vmath::Vcopy(nquad,&tmp[0],1,&outarray[0],1);
	    }
	    else
	    {
		Blas::Dgemv('T',nquad,m_base[0]->GetNumModes(),1.0,base,nquad,
			    &tmp[0],1,0.0,&outarray[0],1);
	    }    
	}
  
	void StdSegExp::IProductWRTBase(const NekDoubleSharedArray &inarray, 
					NekDoubleSharedArray &outarray)
	{
	    IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
	}
	
	void StdSegExp::FillMode(const int mode, NekDoubleSharedArray &outarray)
	{
	    int   nquad = m_base[0]->GetNumPoints();
	    const NekDouble * base  = m_base[0]->GetBdata();
	    
	    ASSERTL2(modes <= m_ncoeffs , 
	     "calling argument mode is larger than total expansion order");

	    Vmath::Vcopy(nquad,(NekDouble *)base+mode*nquad,1, &outarray[0],1);
	}
	
	DNekMatSharedPtr StdSegExp::GenMassMatrix() 
	{
	    DNekMatSharedPtr Mat = StdExpansion::GenerateMassMatrix();

	    // For Fourier basis set the imaginary component of mean mode
	    // to have a unit diagonal component in mass matrix 
	    if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
	    {
		(*Mat)(1,1) = 1.0;
	    }

	    return Mat;
	}
	
	DNekMatSharedPtr StdSegExp::GenLapMatrix() 
	{
	    int    i;
	    int   nquad = m_base[0]->GetNumPoints();
	    const NekDouble * dbase  = m_base[0]->GetDbdata();
	    const NekDouble *z,*w;
	    NekDoubleSharedArray tmp  = GetDoubleTmpSpace(nquad);
	    int nummodes = m_base[0]->GetNumModes();
	    DNekMatSharedPtr Mat;

	    //	    Mat.reset(new DNekMat(nummodes,nummodes,
	    // new double [nummodes*nummodes]));
	    Mat = MemoryManager::AllocateSharedPtr<DNekMat>(nummodes,
		  nummodes,MemoryManager::AllocateArray<NekDouble>(nummodes*nummodes));
		      
	    ExpPointsProperties(0)->GetZW(z,w);
    
	    for(i = 0; i < nummodes; ++i)
	    {
		Vmath::Vcopy(nquad,(NekDouble *)dbase+i*nquad,1, &tmp[0],1);
		IProductWRTBase(m_base[0]->GetDbdata(), tmp, tmp, 0);
		Vmath::Vcopy(nquad,&tmp[0],1,&((*Mat).GetPtr())[0]+ i*m_base[0]->GetNumModes(),1);
	    }
	    
	    return Mat;
	}
	
	//----------------------------
	// Differentiation Methods
	//-----------------------------
	
	void StdSegExp::PhysDeriv(const NekDoubleSharedArray &inarray, 
				  NekDoubleSharedArray &out_d0,
				  NekDoubleSharedArray &out_d2,
				  NekDoubleSharedArray &out_d3)
	{
	    PhysTensorDeriv(inarray,out_d0);
	}
	
	//----------------------------
	// Evaluation Methods
	//----------------------------

	void StdSegExp::BwdTrans(const NekDoubleSharedArray &inarray, 
				 NekDoubleSharedArray &outarray)
	{
	    int  nquad = m_base[0]->GetNumPoints();
	    
	    if(m_base[0]->Collocation())
	    {
		Vmath::Vcopy(nquad,&inarray[0],1,&outarray[0],1);
	    }
	    else
	    {
		const NekDouble *base  = m_base[0]->GetBdata();

		Blas::Dgemv('N',nquad,m_base[0]->GetNumModes(),1.0,base,nquad,
			    &inarray[0],1,0.0,&outarray[0],1);
	    }
	}
	
	void StdSegExp::FwdTrans(const NekDoubleSharedArray &inarray, NekDoubleSharedArray &outarray)
	{
	    if(m_base[0]->Collocation())
	    {
		Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&outarray[0],1);
	    }
	    else{
		IProductWRTBase(inarray,outarray);

		// get Mass matrix
		StdLinSysKey         masskey(eMassMatrix,DetShapeType(),*this);
		DNekLinSysSharedPtr  matsys = m_stdLinSysManager[masskey];

		// solve inverse of system
		DNekVec   v(m_ncoeffs,outarray,eWrapper);
		matsys->Solve(v,v);
	    }
	}
 
	NekDouble StdSegExp::PhysEvaluate(const NekDoubleSharedArray &Lcoord)
	{
	    return StdExpansion1D::PhysEvaluate1D(Lcoord);
	}
    
    
	void StdSegExp::MapTo(EdgeOrientation dir, StdExpMap& Map)
	{
	    
	    if(Map.GetLen() < 2)
	    {
		Map.SetMapMemory(2);
	    }
    
	    switch(m_base[0]->GetBasisType())
	    {
	    case LibUtilities::eGLL_Lagrange:
		{
		    int order = m_base[0]->GetNumModes();
		    if(dir == eForwards)
		    {
			Map[0] = 0;
			Map[1] = order-1;
		    }
		    else
		    {
			Map[0] = order-1;
			Map[1] = 0;
		    }
		}
		break;
	    case LibUtilities::eModified_A:
		
		if(dir == eForwards)
		{
		    Map[0] = 0;
		    Map[1] = 1;
		}
		else
		{
		    Map[0] = 1;
		    Map[1] = 0;
		}
		break;
	    default:
		ASSERTL0(0,"Mapping not defined for this expansion");
	    }
	}    
	
	void StdSegExp::GetCoords(NekDoubleSharedArray &coords)
	{
	    Blas::Dcopy(GetNumPoints(0),ExpPointsProperties(0)->GetZ(),
			1,&coords[0],1);
	}
    
    }//end namespace
}//end namespace

/** 
 * $Log: StdSegExp.cpp,v $
 * Revision 1.22  2007/03/21 20:56:43  sherwin
 * Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
 *
 * Revision 1.21  2007/03/20 16:58:43  sherwin
 * Update to use NekDoubleSharedArray storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.20  2007/03/14 21:24:09  sherwin
 * Update for working version of MultiRegions up to ExpList1D
 *
 * Revision 1.19  2007/02/24 09:07:25  sherwin
 * Working version of stdMatrixManager and stdLinSysMatrix
 *
 * Revision 1.18  2007/02/23 19:26:08  jfrazier
 * General bug fix and formatting.
 *
 * Revision 1.17  2007/02/22 22:02:28  sherwin
 * Update with executing StdMatManager
 *
 * Revision 1.16  2007/02/22 18:11:32  sherwin
 * Version with some create functions introduced for StdMatManagers
 *
 * Revision 1.15  2007/02/21 22:55:16  sherwin
 * First integration of StdMatrixManagers
 *
 * Revision 1.14  2007/02/17 03:40:21  jfrazier
 * Couple changes to reflect additions and corrections to reflect linear algebra calls.
 *
 * Revision 1.13  2007/02/13 09:52:28  sherwin
 * Updates to fix mass matrix inverse issues
 *
 * Revision 1.12  2007/02/12 17:00:20  sherwin
 * Modifcations to make a working version of Project1D
 *
 * Revision 1.11  2007/02/07 12:51:53  sherwin
 * Compiling version of Project1D
 *
 * Revision 1.10  2007/01/28 18:34:24  sherwin
 * More modifications to make Demo Project1D compile
 *
 * Revision 1.9  2007/01/23 23:20:21  sherwin
 * New version after Jan 07 update
 *
 * Revision 1.8  2007/01/20 22:35:21  sherwin
 * Version with StdExpansion compiling
 *
 * Revision 1.7  2007/01/15 15:07:20  pvos
 * updating doxygen documentation
 *
 * Revision 1.6  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.5  2006/08/16 23:34:42  jfrazier
 * *** empty log message ***
 *
 * Revision 1.4  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.3  2006/06/06 15:25:21  jfrazier
 * Removed unreferenced variables and replaced ASSERTL0(false, ....) with
 * NEKERROR.
 *
 * Revision 1.2  2006/06/01 13:43:20  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.52  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.51  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.50  2006/03/05 22:11:03  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.49  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.48  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/ 




