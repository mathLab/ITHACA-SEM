///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansioneD.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion3D.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {
	
	
	StdExpansion3D::StdExpansion3D()
	{
	}
	
	StdExpansion3D::StdExpansion3D(int numcoeffs, const LibUtilities::BasisKey &Ba, 
				       const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc):
	    StdExpansion(numcoeffs,3,Ba,Bb,Bc)
	{
	}
	
	StdExpansion3D::StdExpansion3D(const StdExpansion3D &T):
	    StdExpansion(T)
	{
	}
	
	StdExpansion3D::~StdExpansion3D()
	{ 
	}
	
	void StdExpansion3D::PhysTensorDeriv(const ConstArray<OneD, NekDouble> &inarray, 
					     Array<OneD, NekDouble> &outarray_d0, 
					     Array<OneD, NekDouble> &outarray_d1, 
					     Array<OneD, NekDouble> &outarray_d2)
	{
	    int    i;
	    int    nquad0 = m_base[0]->GetNumPoints();
	    int    nquad1 = m_base[1]->GetNumPoints();
	    int    nquad2 = m_base[2]->GetNumPoints();
	    DNekMatSharedPtr D0,D1,D2;
	    Array<OneD, NekDouble> wsp;
            wsp = Array<OneD, NekDouble>(nquad0*nquad1*nquad2);
	    
	    // copy inarray to wsp in case inarray is used as outarray
            Vmath::Vcopy(nquad0*nquad1*nquad2,&inarray[0],1,&wsp[0],1);
	    
	    D0 = ExpPointsProperties(0)->GetD();
	    D1 = ExpPointsProperties(1)->GetD();
	    D2 = ExpPointsProperties(2)->GetD();
	    
	    // calculate du/dx_0
	    if(outarray_d0.num_elements() > 0)
	    {
		for(i=0; i < nquad2; ++i)
		{
		    Blas::Dgemm('T','N',nquad0,nquad1,nquad0,1.0,&(D0->GetPtr())[0],
				nquad0,&wsp[0]+i*nquad0*nquad1, nquad0,0.0,
				&outarray_d0[0]+i*nquad0*nquad1,nquad0);
		}
	    }
	    
	    // calculate du/dx_1
	    if(outarray_d1.num_elements() > 0 )
	    {
		for(i=0; i < nquad2; ++i)
		{
		    Blas:: Dgemm('N','N',nquad0,nquad1,nquad1,1.0,
				 &wsp[0]+i*nquad0*nquad1,nquad0,&(D1->GetPtr())[0],
				 nquad1,0.0, &outarray_d1[0]+i*nquad0*nquad1,nquad0);
		}
	    }
	    
	    // calculate du/dx_2
	    if(outarray_d2.num_elements() > 0)
	    {
		for(i=0; i < nquad0*nquad1; ++i)
		{
		    Blas:: Dgemv('T',nquad2,nquad2,1.0,&(D2->GetPtr())[0],
				 nquad2, &wsp[0]+i,nquad0*nquad1,0.0,
				 &outarray_d2[0]+i,nquad0*nquad1);
		}
	    }   
	}
	
	NekDouble StdExpansion3D::PhysEvaluate(ConstArray<OneD, NekDouble> &coords)
	{
	    NekDouble  val;
	    int i;
	    int nq1 = m_base[0]->GetNumPoints();
	    int nq2 = m_base[1]->GetNumPoints();
	    int nq3 = m_base[2]->GetNumPoints();
	    	    
            DNekMatSharedPtr I;

	    Array<OneD, NekDouble> wsp  = Array<OneD, NekDouble>(nq2*nq3);
	    Array<OneD, NekDouble> wsp1 = Array<OneD, NekDouble>(nq3);
	    
	    ASSERTL2(coords[0] < -1,"coord[0] < -1");
	    ASSERTL2(coords[0] >  1,"coord[0] >  1");
	    ASSERTL2(coords[1] < -1,"coord[1] < -1");
	    ASSERTL2(coords[1] >  1,"coord[1] >  1");
	    ASSERTL2(coords[2] < -1,"coord[2] < -1");
	    ASSERTL2(coords[2] >  1,"coord[2] >  1");
	    
	    // interpolate first coordinate direction
            I = ExpPointsProperties(0)->GetI(coords);

	    for(i = 0; i < nq2*nq3;++i)
	    {
		wsp[i] =  Blas::Ddot(nq1,&(I->GetPtr())[0],1,&m_phys[0]+i*nq1,1);
	    }
	    
	    // interpolate in second coordinate direction 
            I = ExpPointsProperties(1)->GetI(coords+1);
	    for(i =0; i < nq3; ++i)
	    {
		wsp1[i] = Blas::Ddot(nq2,&wsp[0]+i*nq2,1,&(I->GetPtr())[0],1);
	    }
	    
	    // interpolate in third coordinate direction 
            I = ExpPointsProperties(2)->GetI(coords+2);
	    val = Blas::Ddot(nq3,&wsp1[0],1,&(I->GetPtr())[0],1);
	    
	    return val;    
	}
	
    }//end namespace
}//end namespace

/** 
 * $Log: StdExpansion3D.cpp,v $
 * Revision 1.9  2007/05/15 05:18:23  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.8  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.7  2007/04/04 20:48:17  sherwin
 * Update to handle SharedArrays
 *
 * Revision 1.6  2007/03/29 19:35:09  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.5  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.4  2007/01/18 18:44:45  bnelson
 * Updates to compile on Visual Studio 2005.
 *
 * Revision 1.3  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.2  2006/06/01 14:46:16  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.9  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.8  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.7  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 **/ 



