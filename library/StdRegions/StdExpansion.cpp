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

        /// define list of number of vertices corresponding to each ShapeType
        const int g_shapenverts[SIZE_ShapeType] = {0,2,3,4,4,5,6,8};

        /// define list of number of edges corresponding to each ShapeType
        const int g_shapenedges[SIZE_ShapeType] = {0,1,3,4,6,8,9,12};

        /// define list of number of faces corresponding to each ShapeType
        const int g_shapenfaces[SIZE_ShapeType] = {0,0,0,0,4,5,5,6};

        StdExpansion::StdExpansion(void): 
        m_numbases(0), 
            m_ncoeffs(0),
            m_coeffs(NULL),
            m_owncoeffs(false),
            m_phys(NULL),
            m_ownphys(false)
        {
        }


        StdExpansion::StdExpansion(int numbases, const BasisKey &Ba, 
            const BasisKey &Bb, const BasisKey &Bc,int numcoeffs, double *coeffs, 
            double *phys, bool spaceowner):
        m_numbases(numbases)
        {
            m_base     = new const Basis*[m_numbases];

            switch(m_numbases)
            {
            case 3:
                ASSERTL2(Bc==NULL,"NULL Basis attempting to be used.");

                m_base[2] = BasisManagerSingleton::Instance().GetBasis(Bc.GetBasisType(),
                    Bc.GetBasisOrder(),Bc.GetPointsType(), Bc.GetPointsOrder(),
                    Bc.GetAlpha(), Bc.GetBeta());

            case 2:
                ASSERTL2(Bb==NULL,"NULL Basis attempting to be used.");

                m_base[1] = BasisManagerSingleton::Instance().GetBasis(Bb.GetBasisType(), 
                    Bb.GetBasisOrder(),Bb.GetPointsType(), Bb.GetPointsOrder(),
                    Bb.GetAlpha(), Bb.GetBeta());
            case 1:
                ASSERTL2(Ba==NULL,"NULL Basis attempting to be used.");

                m_base[0] = BasisManagerSingleton::Instance().GetBasis(Ba.GetBasisType(),
                    Ba.GetBasisOrder(),Ba.GetPointsType(), Ba.GetPointsOrder(),
                    Ba.GetAlpha(), Ba.GetBeta());
                break;
            default:
                ASSERTL0(false, "numbases incorrectly specified");
            };

            //allocate memory for coeffs
            m_ncoeffs = numcoeffs;

            if(spaceowner)
            {
                int i,cnt;

                m_coeffs = new double[m_ncoeffs];

                Vmath::Zero(m_ncoeffs,m_coeffs,1);
                m_owncoeffs = true;

                //allocate memory for phys
                cnt = GetPointsOrder(0);

                for(i=1; i<numbases; ++i)
                {
                    cnt *= GetPointsOrder(i);
                }

                m_phys = new double[cnt];
                m_ownphys = true;
            }
            else
            {
                m_coeffs = coeffs;
                m_owncoeffs = false;
                m_phys = phys;
                m_ownphys = false;
            }

        }//end constructor


        StdExpansion::StdExpansion(const StdExpansion &T):
        m_numbases(T.m_numbases)
        {
            int i,j;
            m_base = new const Basis*[m_numbases];

            for(j=0; j<m_numbases; j++)
            {
                m_base[j] = T.m_base[j];
            }


            // NOTE: Copy Constructor produces a deep copy
            // allocate memory for coeffs
            // need to check allocation for variable order. 
            m_ncoeffs = T.m_ncoeffs;
            m_coeffs = new double[m_ncoeffs];
            for(i=0; i<m_ncoeffs; i++)
            {
                m_coeffs[i] = T.m_coeffs[i];
            }

            m_owncoeffs = true;  


            //allocate memory for phys
            int numphys = GetPointsOrder(0)*GetPointsOrder(1);
            m_phys = new double[numphys];
            m_ownphys = true;
            for(j=0; j < numphys; j++)
            {
                m_phys[j] = T.m_phys[j];
            }
        }

        StdExpansion::~StdExpansion()
        {
            if(m_base)
            {
                delete[] m_base;
                m_base = NULL;
            }

            if(m_owncoeffs && m_coeffs)
            {
                delete[] m_coeffs;
                m_coeffs = (double*)NULL;
            }

            if(m_ownphys && m_phys)
            {
                delete[] m_phys;
                m_phys = (double*)NULL;
            }
        }

        /// \brief Function to evaluate the discrete \f$ L_\infty\f$
        /// error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
        ///	u_{exact}\f$ is given by the array \a sol. 
        ///
        // 	The full function is defined in StdExpansion::Linf 
        ///
        ///	Input: 
        ///
        ///	- \a _phys: takes the physical value space array as
        ///      approximate solution
        ///
        ///  - \a sol: array of solution function  at physical quadrature points
        ///
        ///      output: 
        ///
        ///      - returns the \f$ L_\infty \f$ error as a double. 

        double StdExpansion::Linf(const double *sol)
        {
            int     i,ntot = 1;
            double  val;
            double *tmp;
            BstShrDArray  wsp;

            for(i=0; i<m_numbases; ++i)
            {
                ntot *= m_base[i]->GetPointsOrder();
            }

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vsub(ntot,sol,1,m_phys,1,tmp,1);
            Vmath::Vabs(ntot,tmp,1,tmp,1);
            val = Vmath::Vamax(ntot,tmp,1);    

            return  val;
        }

        /** \brief Function to evaluate the \f$ L_\infty \f$ norm of
        the function defined at the physical points \a (this)->_phys. 

        The full function is defined in StdExpansion::Linf 

        Input: 

        - \a _phys: uses the physical value space array as discrete
        function to be evaulated.

        output: 

        - returns the \f$ L_\infty \f$  as a double. 
        */

        double StdExpansion::Linf()
        {
            int  i,ntot = 1;

            for(i=0; i<m_numbases; ++i)
            {
                ntot *= m_base[i]->GetPointsOrder();
            }

            return Vmath::Vamax(ntot,m_phys,1);    
        }

        /** \brief Function to evaluate the \f$ L_2\f$, \f$ | \epsilon
        |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2 dx \right]^{1/2}
        d\xi_1 \f$ where \f$ u_{exact}\f$ is given by the array sol.

        The full function is defined in StdExpansion::L2 

        Input: 

        - \a _phys: takes the physical value space array as
        approximate solution
        - \a sol: array of solution function  at physical quadrature points

        output: 

        - returns the \f$ L_2 \f$ error as a double. 
        */

        double StdExpansion::L2(const double *sol)
        {
            int     i,ntot = 1;
            double  val;
            double *tmp;
            BstShrDArray wsp;

            for(i=0; i<m_numbases; ++i)
            {
                ntot *= m_base[i]->GetPointsOrder();
            }

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vsub(ntot, sol, 1, m_phys, 1, tmp, 1);
            Vmath::Vmul(ntot, tmp, 1, tmp, 1, tmp, 1);

            val = sqrt(v_Integral(tmp));

            return val;
        }


        /// \brief Function to evaluate the \f$ L_2\f$ norm of the
        ///  function defined at the physical points \a (this)->_phys.  
        ///    
        ///  The full function is defined in StdExpansion::L2 
        ///
        ///      
        /// \param _phys: uses the physical value space array as discrete
        ///    function to be evaulated.
        ///    
        ///  \return value of the \f$ L_2 \f$  as a double. 

        double StdExpansion::L2()
        {
            int     i,ntot = 1;
            double  val;
            double *tmp;
            BstShrDArray wsp;

            for(i=0; i<m_numbases; ++i)
            {
                ntot *= m_base[i]->GetPointsOrder();
            }

            wsp = GetDoubleTmpSpace(ntot);
            tmp = wsp.get();

            Vmath::Vmul(ntot, m_phys, 1, m_phys, 1, tmp, 1);
            val   = sqrt(v_Integral(tmp));

            return val;
        }

        void StdExpansion::GenerateMassMatrix(double *outarray)
        {
            int     i;
            BstShrDArray tmp = GetDoubleTmpSpace(GetPointsTot());

            for(i=0; i<m_ncoeffs; ++i)
            {
                v_FillMode(i, tmp.get());
                v_IProductWRTBase(tmp.get(), outarray+i*m_ncoeffs);
            }
        }



        // 2D Interpolation
        void StdExpansion::Interp2D(const  BasisKey *fbasis0, 
            const BasisKey *fbasis1, const double *from,  const BasisKey *tbasis0,
            const BasisKey* tbasis1, double *to)
        {
            const double *I0,*I1;
            double *tmp;
            BstShrDArray wsp = GetDoubleTmpSpace(tbasis1->GetPointsOrder()*
                fbasis0->GetPointsOrder());
            tmp = wsp.get();

            BasisManagerSingleton::Instance().GetI(fbasis0, tbasis0, I0);
            BasisManagerSingleton::Instance().GetI(fbasis1, tbasis1, I1);

            Blas::Dgemm('T', 'T', tbasis1->GetPointsOrder(), fbasis0->GetPointsOrder(),
                fbasis1->GetPointsOrder(), 1.0, I1,  fbasis1->GetPointsOrder(),
                (double *) from,fbasis0->GetPointsOrder(), 0.0, tmp,
                tbasis1->GetPointsOrder());

            Blas::Dgemm('T', 'T',tbasis0->GetPointsOrder(),tbasis1->GetPointsOrder(),
                fbasis0->GetPointsOrder(),1.0,I0, fbasis0->GetPointsOrder(),
                tmp, tbasis1->GetPointsOrder(),0.0,to,
                tbasis0->GetPointsOrder());
        }

        // 1D Interpolation
        void StdExpansion::Interp1D(const  BasisKey *fbasis0, const double *from,  
            const BasisKey *tbasis0, double *to)
        {
            const double *I0;

            BasisManagerSingleton::Instance().GetI(fbasis0, tbasis0, I0);

            Blas::Dgemv('T', fbasis0->GetPointsOrder(), tbasis0->GetPointsOrder(), 
                1.0, I0, fbasis0->GetPointsOrder(), from, 1, 0.0, to, 1);
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
