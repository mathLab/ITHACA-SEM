///////////////////////////////////////////////////////////////////////////////
//
// File Stdexpansion.cpp
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
        const int g_shapenverts[SIZE_ShapeType] = {0,2,3,4,4,5,6,8};

        /** define list of number of edges corresponding to each ShapeType */
        const int g_shapenedges[SIZE_ShapeType] = {0,1,3,4,6,8,9,12};

        /** define list of number of faces corresponding to each ShapeType */
        const int g_shapenfaces[SIZE_ShapeType] = {0,0,0,0,4,5,5,6};

        StdExpansion::StdExpansion(void): 
	    m_numbases(0), 
            m_ncoeffs(0)
        {
        }


        StdExpansion::StdExpansion(int numbases, const LibUtilities::BasisKey &Ba, 
	   const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc,
	   int numcoeffs, double *coeffs, double *phys):
	    m_numbases(numbases)
	{
	    m_base = MemoryManager::AllocateArray<LibUtilities::BasisSharedPtr>(m_numbases);
	    
            switch(m_numbases)
            {
            case 3:
                ASSERTL2(Bc==NULL,"NULL Basis attempting to be used.");

                m_base[2] = LibUtilities::BasisManager()[Bc];

            case 2:
                ASSERTL2(Bb==NULL,"NULL Basis attempting to be used.");

                m_base[1] = LibUtilities::BasisManager()[Bb];
            case 1:
                ASSERTL2(Ba==NULL,"NULL Basis attempting to be used.");

                m_base[0] = LibUtilities::BasisManager()[Ba];
                break;
            default:
                ASSERTL0(false, "numbases incorrectly specified");
            };

            //allocate memory for coeffs
            m_ncoeffs = numcoeffs;
	    m_coeffs = MemoryManager::AllocateSharedArray<double>(m_ncoeffs);
	    Vmath::Zero(m_ncoeffs,&m_coeffs[0],1);

	    //allocate memory for phys
	    m_phys = MemoryManager::AllocateSharedArray<double>(GetTotPoints());
			 

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T):
        m_numbases(T.m_numbases)
        {
            int i,j;

	    m_base = MemoryManager::AllocateArray<LibUtilities::BasisSharedPtr>(m_numbases);

            for(j=0; j<m_numbases; j++)
            {
                m_base[j] = T.m_base[j];
            }


            // NOTE: Copy Constructor produces a deep copy
            // allocate memory for coeffs
            // need to check allocation for variable order. 
            m_ncoeffs = T.m_ncoeffs;
	    m_coeffs = MemoryManager::AllocateSharedArray<double>(m_ncoeffs);
            for(i=0; i<m_ncoeffs; i++)
            {
                m_coeffs[i] = T.m_coeffs[i];
            }

            //allocate memory for phys
            int numphys = GetTotPoints();
	    m_phys = MemoryManager::AllocateSharedArray<double>(GetTotPoints());
            for(j=0; j < numphys; j++)
            {
                m_phys[j] = T.m_phys[j];
            }
        }

        StdExpansion::~StdExpansion()
        {
        }

        double StdExpansion::Linf(const double *sol)
        {
            int     ntot;
            double  val;
            double *tmp;
            BstShrDArray  wsp;

	    ntot =  GetTotPoints();

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vsub(ntot,sol,1,&m_phys[0],1,tmp,1);
            Vmath::Vabs(ntot,tmp,1,tmp,1);
            val = Vmath::Vamax(ntot,tmp,1);    

            return  val;
        }

        double StdExpansion::Linf()
        {
            return Vmath::Vamax(GetTotPoints(),&m_phys[0],1);    
        }

        double StdExpansion::L2(const double *sol)
        {
            int     ntot = GetTotPoints();
            double  val;
            double *tmp;
            BstShrDArray wsp;

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vsub(ntot, sol, 1, &m_phys[0], 1, tmp, 1);
            Vmath::Vmul(ntot, tmp, 1, tmp, 1, tmp, 1);
	    
            val = sqrt(v_Integral(tmp));

            return val;
        }

        double StdExpansion::L2()
        {
            int     ntot = GetTotPoints();
            double  val;
            double *tmp;
            BstShrDArray wsp;

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vmul(ntot, &m_phys[0], 1, &m_phys[0], 1, tmp, 1);
            val   = sqrt(v_Integral(tmp));

            return val;
        }

        DNekMatSharedPtr StdExpansion::GenerateMassMatrix()
        {
            int     i;
            BstShrDArray store = GetDoubleTmpSpace(m_ncoeffs);
            BstShrDArray tmp   = GetDoubleTmpSpace(GetTotPoints());

	    DNekMatSharedPtr  Mat;
	    
	    Mat = MemoryManager::AllocateSharedPtr<DNekMat>(m_ncoeffs,m_ncoeffs,MemoryManager::AllocateArray<double>(m_ncoeffs*m_ncoeffs));
	    
	    Blas::Dcopy(m_ncoeffs,&m_coeffs[0],1,&store[0],1);
            for(i=0; i<m_ncoeffs; ++i)
            {
                v_FillMode(i, &tmp[0]);
                v_IProductWRTBase(&tmp[0], &((*Mat).GetPtr())[0]+i*m_ncoeffs);
            }
	    Blas::Dcopy(m_ncoeffs,&store[0],1,&m_coeffs[0],1);
        }



        // 2D Interpolation
        void StdExpansion::Interp2D(const  LibUtilities::BasisKey *fbasis0, 
            const LibUtilities::BasisKey *fbasis1, const double *from,  
	    const LibUtilities::BasisKey *tbasis0,
            const LibUtilities::BasisKey* tbasis1, double *to)
        {
	    DNekMatSharedPtr I0,I1;
            double *tmp;
            BstShrDArray wsp = GetDoubleTmpSpace(tbasis1->GetNumPoints()*
                fbasis0->GetNumPoints());
            tmp = wsp.get();

	    I0 = LibUtilities::PointsManager()[fbasis0->GetPointsKey()]
		->GetI(tbasis0->GetPointsKey());
	    I1 = LibUtilities::PointsManager()[fbasis1->GetPointsKey()]
		->GetI(tbasis1->GetPointsKey());

            Blas::Dgemm('T', 'T', tbasis1->GetNumPoints(), fbasis0->GetNumPoints(),
                fbasis1->GetNumPoints(), 1.0, &((*I1).GetPtr())[0],  
			fbasis1->GetNumPoints(),
			(double *) from,fbasis0->GetNumPoints(), 0.0, tmp,
			tbasis1->GetNumPoints());

            Blas::Dgemm('T', 'T',tbasis0->GetNumPoints(),tbasis1->GetNumPoints(),
			fbasis0->GetNumPoints(),1.0,&((*I0).GetPtr())[0],
			fbasis0->GetNumPoints(),tmp, tbasis1->GetNumPoints(),
			0.0,to, tbasis0->GetNumPoints());
        }

        // 1D Interpolation
        void StdExpansion::Interp1D(const  LibUtilities::BasisKey *fbasis0, const double *from,  
            const LibUtilities::BasisKey *tbasis0, double *to)
        {
	    DNekMatSharedPtr I0;

	    I0 = LibUtilities::PointsManager()[fbasis0->GetPointsKey()]
		->GetI(tbasis0->GetPointsKey());

            Blas::Dgemv('T', fbasis0->GetNumPoints(), tbasis0->GetNumPoints(), 
			1.0, &((*I0).GetPtr())[0], fbasis0->GetNumPoints(), 
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
