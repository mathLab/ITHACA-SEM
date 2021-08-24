///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion0D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// which are common to 0d expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion0D.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdExpansion0D::StdExpansion0D()
        {
        }
	
        StdExpansion0D::StdExpansion0D(int numcoeffs, const LibUtilities::BasisKey &Ba):
            StdExpansion(numcoeffs,1,Ba)
        {
        }
		
        StdExpansion0D::StdExpansion0D(const StdExpansion0D &T):StdExpansion(T)
        {
        }
	
        StdExpansion0D::~StdExpansion0D()
        {
        }
        
	
        //----------------------------
        // Differentiation Methods
        //-----------------------------
	
        void StdExpansion0D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble>& outarray)
        {
            int nquad = GetTotPoints();
            DNekMatSharedPtr D = m_base[0]->GetD();

            if( inarray.data() == outarray.data())
            {
                Array<OneD, NekDouble> wsp(nquad);
                CopyArray(inarray, wsp);
                Blas::Dgemv('N',nquad,nquad,1.0,&(D->GetPtr())[0],nquad,
                            &wsp[0],1,0.0,&outarray[0],1);
			}
            else
            {
                Blas::Dgemv('N',nquad,nquad,1.0,&(D->GetPtr())[0],nquad,
                            &inarray[0],1,0.0,&outarray[0],1);
            }
        }
	
        NekDouble StdExpansion0D::v_PhysEvaluate(const Array<OneD, const NekDouble>& Lcoord, const Array<OneD, const NekDouble>& physvals)
        {
            int    nquad = GetTotPoints();
            NekDouble  val;
            DNekMatSharedPtr I = m_base[0]->GetI(Lcoord);
            
            ASSERTL2(Lcoord[0] >= -1,"Lcoord[0] < -1");
            ASSERTL2(Lcoord[0] <=  1,"Lcoord[0] >  1");
            
            val = Blas::Ddot(nquad, I->GetPtr(), 1, physvals, 1);
            
            return val;
        }
	
    }//end namespace
}//end namespace
