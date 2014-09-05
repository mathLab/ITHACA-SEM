///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion1D.cpp
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
// which are common to 1d expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion1D.h>

namespace Nektar
{
    namespace StdRegions
    {

    StdExpansion1D::StdExpansion1D()
    {
    }

    StdExpansion1D::StdExpansion1D(int numcoeffs, const LibUtilities::BasisKey &Ba):
        StdExpansion(numcoeffs,1,Ba)
    {
    }

    StdExpansion1D::StdExpansion1D(const StdExpansion1D &T):StdExpansion(T)
    {
    }

    StdExpansion1D::~StdExpansion1D()
    {
    }


    //----------------------------
    // Differentiation Methods
    //-----------------------------

    void StdExpansion1D::PhysTensorDeriv(const Array<OneD, const NekDouble>& inarray,
                         Array<OneD, NekDouble>& outarray)
    {
        int nquad = GetTotPoints();
        DNekMatSharedPtr D = m_base[0]->GetD();

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

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

#else //NEKTAR_USING_DIRECT_BLAS_CALLS

        NekVector<NekDouble> out(nquad,outarray,eWrapper);

        if(inarray.data() == outarray.data()) // copy intput array
        {
            NekVector<NekDouble> in(nquad,inarray,eCopy);
            out = (*D)*in;
        }
        else
        {
            NekVector<NekDouble> in (nquad,inarray,eWrapper);
            out = (*D)*in;
        }

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS
    }

        NekDouble StdExpansion1D::v_PhysEvaluate(
                const Array<OneD, const NekDouble>& Lcoord,
                const Array<OneD, const NekDouble>& physvals)
        {
        int    nquad = GetTotPoints();
        NekDouble  val;
        DNekMatSharedPtr I = m_base[0]->GetI(Lcoord);

        ASSERTL2(Lcoord[0] >= -1,"Lcoord[0] < -1");
        ASSERTL2(Lcoord[0] <=  1,"Lcoord[0] >  1");

        val = Blas::Ddot(nquad, I->GetPtr(), 1, physvals, 1);

        return val;
    }
	
	void StdExpansion1D::v_SetUpPhysNormals(const int vertex)
    {
		ComputeVertexNormal(vertex);
    }

        const NormalVector & StdExpansion1D::v_GetSurfaceNormal(const int id) const
        {
            return v_GetVertexNormal(id);
        }

		
	const NormalVector & StdExpansion1D::v_GetVertexNormal(const int vertex) const
    {
         std::map<int, NormalVector>::const_iterator x;
         x = m_vertexNormals.find(vertex);
         ASSERTL0 (x != m_vertexNormals.end(),
				  "vertex normal not computed.");
         return x->second;
	}

    }//end namespace
}//end namespace

/**
 * $Log: StdExpansion1D.cpp,v $
 * Revision 1.28  2008/08/20 09:14:57  sherwin
 * In TensorDeriv replaced comparison of arrays with comparison of stored pointer
 *
 * Revision 1.27  2008/08/18 08:31:26  sherwin
 * Update to fix PhysTensorDeriv  when inarray == outarray
 *
 * Revision 1.26  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.25  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.24  2008/04/06 06:04:14  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.23  2008/04/03 16:12:11  pvos
 * updates for NEKTAR_USING_DIRECT_BLAS_CALLS
 *
 * Revision 1.22  2008/03/12 15:25:09  pvos
 * Clean up of the code
 *
 * Revision 1.21  2007/11/29 21:40:22  sherwin
 * updates for MultiRegions and DG solver
 *
 * Revision 1.20  2007/11/08 14:23:33  ehan
 * Removed white space
 *
 * Revision 1.19  2007/07/20 02:16:53  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.18  2007/05/28 08:35:26  sherwin
 * Updated for localregions up to Project1D
 *
 * Revision 1.17  2007/05/17 17:59:28  sherwin
 * Modification to make Demos work after introducion of Array<>
 *
 * Revision 1.16  2007/05/15 05:18:22  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.15  2007/04/26 15:00:17  sherwin
 * SJS compiling working version using SHaredArrays
 *
 * Revision 1.14  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.13  2007/04/08 03:36:58  jfrazier
 * Updated to use SharedArray consistently and minor reformatting.
 *
 * Revision 1.12  2007/04/04 20:48:17  sherwin
 * Update to handle SharedArrays
 *
 * Revision 1.11  2007/03/29 19:35:09  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.10  2007/03/20 16:58:42  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.9  2007/03/14 21:24:09  sherwin
 * Update for working version of MultiRegions up to ExpList1D
 *
 * Revision 1.8  2007/02/07 12:51:53  sherwin
 * Compiling version of Project1D
 *
 * Revision 1.7  2007/01/30 20:01:35  sherwin
 * Update for first compiling Project1D routine
 *
 * Revision 1.6  2007/01/28 18:34:22  sherwin
 * More modifications to make Demo Project1D compile
 *
 * Revision 1.5  2007/01/21 02:28:08  sherwin
 * Compiling under new revision
 *
 * Revision 1.4  2007/01/20 22:35:21  sherwin
 * Version with StdExpansion compiling
 *
 * Revision 1.3  2007/01/15 11:30:20  pvos
 * Updating doxygen documentation
 *
 * Revision 1.2  2006/06/01 14:46:16  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:30  kirby
 * *** empty log message ***
 *
 * Revision 1.15  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.14  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.13  2006/03/13 18:29:35  sherwin
 *
 * Corrected error with definition of GetCoords
 *
 * Revision 1.12  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 * Revision 1.11  2006/02/26 23:37:29  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/

