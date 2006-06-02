///////////////////////////////////////////////////////////////////////////////
//
// File StdTriExp.cpp
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
// Description: Triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdTriExp.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdMatrix StdTriExp::s_elmtmats;

        StdTriExp::StdTriExp() // default constructor of StdExpansion is directly called. 
        {
        } //default constructor


        StdTriExp::StdTriExp(const BasisKey &Ba, const BasisKey &Bb):
        StdExpansion2D(Ba,Bb,Ba.GetBasisOrder()*(Ba.GetBasisOrder()+1)/2+
            Ba.GetBasisOrder()*(Bb.GetBasisOrder()-Ba.GetBasisOrder()), 
            NULL,NULL,true)
        {    

            if(Ba.GetBasisOrder() >  Bb.GetBasisOrder())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'b' direction");
            }
        }

        StdTriExp::StdTriExp(const BasisKey &Ba,  const BasisKey &Bb, 
            double *coeffs,  double *phys):
        StdExpansion2D(Ba,Bb,Ba.GetBasisOrder()*(Ba.GetBasisOrder()+1)/2+
            Ba.GetBasisOrder()*(Bb.GetBasisOrder()-Ba.GetBasisOrder())
            ,coeffs,phys,false)
        {    
            if(Ba.GetBasisOrder() >  Bb.GetBasisOrder())
            {
                ASSERTL0(false, "order in 'a' direction is higher than order in 'b' direction");
            }
        }

        StdTriExp::StdTriExp(const StdTriExp &T):
        StdExpansion2D(T)
        {
        }

        StdTriExp::~StdTriExp()
        {
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////
        /** 
        This is just a wrapper around StdExpansoin2D which multiplies by
        0.5 which is due to the factor (1-\xi_2)/2 in the integral
        weight
        */
        double StdTriExp::Integral(const double *inarray)
        {
            const double *z0,*z1,*w0,*w1;
            int    i,nquad1 = m_base[1]->GetPointsOrder();
            BstShrDArray wsp = GetDoubleTmpSpace(nquad1);
            double  *w1_tmp  = wsp.get();

            BasisManagerSingleton::Instance().GetZW(m_base[0],z0,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z1,w1);

            switch((int)m_base[1]->GetAlpha())
            {
            case 0: // Legendre inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    w1_tmp[i] = 0.5*(1-z1[i])*w1[i];
                }
                break;
            case 1: // (1,0) Jacobi Inner product 
                Vmath::Smul(nquad1,0.5,(double *)w1,1,w1_tmp,1);      
                break;
            }

            return StdExpansion2D::Integral(inarray,w0,w1_tmp);
        }


        void StdTriExp::IProductWRTBase(const double * inarray, double * outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),inarray,
                outarray);
        }

        /** \brief Calculate the inner product of inarray with respect to
        the basis B=base0[p]*base1[pq] and put into outarray:

        \f$ 
        \begin{array}{rcl}
        I_{pq} = (\phi^A_q \phi^B_{pq}, u) &=& \sum_{i=0}^{nq_0}\sum_{j=0}^{nq_1}
        \phi^A_p(\eta_{0,i})\phi^B_{pq}(\eta_{1,j}) w^0_i w^1_j 
        u(\xi_{0,i} \xi_{1,j})\\
        & = & \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i})
        \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) \tilde{u}_{i,j} 
        \end{array}
        \f$ 

        where

        \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$

        which can be implemented as

        \f$  f_{pj} = \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i}) \tilde{u}_{i,j} 
        \rightarrow {\bf B_1 U}  \f$
        \f$  I_{pq} = \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) f_{pj} 
        \rightarrow {\bf B_2[p*skip] f[skip]}  \f$

        \b Recall: \f$ \eta_{1} = \frac{2(1+\xi_1)}{(1-\xi_2)}-1, \, 
        \eta_2 = \xi_2\f$

        \b Note: For the orthgonality of this expansion to be realised
        the 'q' ordering must run fastest in contrast to the Quad and
        Hex ordering where 'p' index runs fastest to be consistent with
        the quadrature ordering. 

        In the triangular space the i (i.e. \f$\eta_1\f$ direction)
        ordering still runs fastest by convention.
        **/

        void StdTriExp:: IProductWRTBase(const double *base0, const double *base1, 
            const double *inarray, double *outarray)
        {
            int    i,mode;
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            int    order0 = m_base[0]->GetBasisOrder();
            int    order1 = m_base[1]->GetBasisOrder();
            const  double *z0,*z1,*w0,*w1;
            BstShrDArray tmp  = GetDoubleTmpSpace(nquad0*nquad1);
            BstShrDArray tmp1 = GetDoubleTmpSpace(nquad0*nquad1);


            BasisManagerSingleton::Instance().GetZW(m_base[0],z0,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z1,w1);

            ASSERTL2((m_base[1]->GetBasisType() == eOrtho_B)||
                (m_base[1]->GetBasisType() == eModified_B), 
                "Basis[1] is not of general tensor type");

            ASSERTL2((m_base[0]->GetAlpha() == 0.0)&&(m_base[1]->GetAlpha() > 1.0),
                "Basis[0] has illegal alpha weight");

            ASSERTL2((m_base[1]->GetBeta() == 0.0)&&(m_base[1]->GetBeta() == 0.0),
                "Basis[1] has non-zero beta weight");

            // Note cannot use outarray as tmp space since dimensions are not always
            // guarenteed to be sufficient 

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,(double*)inarray+i*nquad0,1,(double*)w0,1,
			    tmp.get()+i*nquad0,1);
            }

            switch((int)m_base[1]->GetAlpha())
            {
            case 0: // Legendre inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i], tmp.get()+i*nquad0,1);
                }
                break;
            case 1: // (1,0) Jacobi Inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*w1[i], tmp.get()+i*nquad0,1);      
                }
                break;
            }

            // Inner product with respect to 'a' direction 
            Blas::Dgemm('T','N',nquad1,order0,nquad0,1.0,tmp.get(),nquad0,base0,nquad0,
                0.0,tmp1.get(),nquad1);

            // Inner product with respect to 'b' direction 
            for(mode=i=0; i < order0; ++i)
            {
                Blas::Dgemv('T',nquad1,order1-i,1.0, base1+mode*nquad1,nquad1,
                    tmp1.get()+i*nquad1,1, 0.0, outarray + mode,1);
                mode += order1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == eModified_A)
            {
                outarray[1] += Blas::Ddot(nquad1,base1+nquad1,1,tmp1.get()+nquad1,1);
            }

        }


        /** Fill outarray with mode 'mode' of expansion

        Note for quadrilateral expansions _base[0] (i.e. p)  modes run fastest

        */

        void StdTriExp::FillMode(const int mode, double *outarray)
        {
            int    i,m;
            int   nquad0 = m_base[0]->GetPointsOrder();
            int   nquad1 = m_base[1]->GetPointsOrder();
            int   order0 = m_base[0]->GetBasisOrder();
            int   order1 = m_base[1]->GetBasisOrder();
            const double * base0  = m_base[0]->GetBdata();
            const double * base1  = m_base[1]->GetBdata();
            int   mode0;

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            m= order1;
            for(i = 0; i < order0; ++i, m+=order1-i)
            {
                if(m > mode)
                {
                    mode0 = i;
                    break;
                }
            }

            // deal with top vertex mode in modified basis
            if((mode == 1)&&(m_base[0]->GetBasisType() == eModified_A))
            {
                Vmath::Fill(nquad0*nquad1,1.0,outarray,1);
            }
            else
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad0,(double *)(base0 + mode0*nquad0),1,
                        outarray+i*nquad0,1);
                }
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(double *)(base1 + mode*nquad1),1,outarray+i,
                    nquad0,outarray+i,nquad0);
            }
        }


        StdMatContainer * StdTriExp::GetMassMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalMass(this);
            return tmp;
        }

        StdMatContainer * StdTriExp::GetLapMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalLap(this);
            return tmp;
        }

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        /** 
        \brief Calculate the deritive of the physical points 

        \f$ \frac{\partial u}{\partial  x_1} =  \left . 
        \frac{2.0}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
        \right |_{\eta_2}\f$

        \f$ \frac{\partial u}{\partial  x_2} =  \left . 
        \frac{1+\eta_1}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
        \right |_{\eta_2}  + \left . \frac{\partial u}{\partial d\eta_2}
        \right |_{\eta_1}  \f$

        **/

        void StdTriExp::Deriv(const double *inarray, double *outarray_d0, 
            double *outarray_d1)
        {
            int    i;
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            const  double *D0,*D1,*z0,*z1,*w;
            double *d0, *d1;
            BstShrDArray wsp;
            BstShrDArray wsp1  = GetDoubleTmpSpace(nquad0*nquad1);
            double *gfac = wsp1.get();

            d0 = outarray_d0;
            d1 = outarray_d1;

            BasisManagerSingleton::Instance().GetZW(m_base[0],z0,w);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z1,w);

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac[i] = 2.0/(1-z1[i]);
            }

            if(!outarray_d1)// if no d1 required do not need to calculate both deriv
            {
                TensorDeriv(inarray, d0,d1);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],d0+i*nquad0,1);
                }
            }
            else
            {
                if(!outarray_d0)// need other local callopsed derivative for d1 
                {
                    wsp = GetDoubleTmpSpace(nquad0*nquad1);
                    d0 = wsp.get();
                }

                TensorDeriv(inarray, d0,d1);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],d0+i*nquad0,1);
                }

                // set up geometric factor: (1_z0)/(1-z1)
                for(i = 0; i < nquad0; ++i)
                {
                    gfac[i] = 0.5*(1+z0[i]);
                }

                for(i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,gfac,1,d0+i*nquad0,1,d1+i*nquad0,1,
                        d1+i*nquad0,1);
                }	
            }
        }

        /** 
        \brief Calculate the deritive of the physical points 

        **/

        inline void StdTriExp::Deriv(double *outarray_d0, double *outarray_d1)
        {
            Deriv(this->m_phys,outarray_d0,outarray_d1);
        }

        ///////////////////////////////
        // Evaluation Methods
        ///////////////////////////////

        /** 
        \brief Backward tranform for triangular elements

        Note: That 'q' (base[1]) runs fastest in this element 
        */

        void StdTriExp::BwdTrans(double * outarray)
        {
            int           i,mode;
            int           nquad0 = m_base[0]->GetPointsOrder();
            int           nquad1 = m_base[1]->GetPointsOrder();
            int           order0 = m_base[0]->GetBasisOrder();
            int           order1 = m_base[1]->GetBasisOrder();
            const double *base0  = m_base[0]->GetBdata();
            const double *base1  = m_base[1]->GetBdata();
            BstShrDArray tmp  = GetDoubleTmpSpace(order0*nquad1);


            ASSERTL2((m_base[1]->GetBasisType() != eOrtho_B)||
                (m_base[1]->GetBasisType() != eModified_B),
                "Basis[1] is not of general tensor type");

            for(i = mode = 0; i < order0; ++i)
            {
                Blas::Dgemv('N', nquad1,order1-i,1.0,base1+mode*nquad1,
                    nquad1,m_coeffs+mode,1,0.0,tmp.get()+i*nquad1,1);
                mode += order1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == eModified_A)
            {
                Blas::Daxpy(nquad1,m_coeffs[1],base1+nquad1,1,tmp.get()+nquad1,1);
            }

            Blas::Dgemm('N','T', nquad0,nquad1,order0,1.0, base0,nquad0, 
                tmp.get(), nquad1,0.0,outarray, nquad0);
        }

        void StdTriExp::FwdTrans(const double *inarray)
        {
            StdMatContainer *M;

            IProductWRTBase(inarray,m_coeffs);
            M = GetMassMatrix();
            M->Solve(m_coeffs,1);
        }

        /// Single Point Evaluation
        double StdTriExp::Evaluate(const double * coords)
        {
            double coll[2];

            // set up local coordinate system 
            if((fabs(coords[0]+1.0) < NekConstants::kEvaluateTol)
                &&(fabs(coords[1]-1.0) < NekConstants::kEvaluateTol))
            {
                coll[0] = 0.0;
                coll[1] = 1.0;
            }
            else
            {
                coll[0] = 2*(1+coords[0])/(1-coords[1])-1.0; 
                coll[1] = coords[1]; 
            }

            return  PhysEvaluate(coll); 
        }


        void  StdTriExp::MapTo(StdSegExp& edge, const int eid, 
            const EdgeOrientation eorient, StdExpMap &Map)
        {

            int i, start, skip;
            const int ncoeffs = edge.GetNcoeffs();
            int *dir, order0,order1;
            BasisType Btype = edge.GetBasisType(0);
            BstShrIArray wsp; 

            ASSERTL2(eid>=0&&eid<=2,"eid must be between 0 and 2");
            ASSERTL2(Btype == eModified_A,"Mapping only set up "
                "for Modified_A edges");
            ASSERTL2(Btype == m_base[0]->GetBasisType(),
                "Expansion type of edge and StdQuadExp are different");

            // make sure have correct memory storage
            if(ncoeffs != Map.GetLen())
            {
                Map.SetMap(ncoeffs);
            }

            order0 = m_base[0]->GetBasisOrder();
            order1 = m_base[0]->GetBasisOrder();

            wsp = GetIntTmpSpace(ncoeffs);
            dir = wsp.get(); 

            if(eorient == eForwards)
            {
                for(i = 0; i < ncoeffs; ++i)
                {
                    dir[i] = i;
                }
            }
            else
            {
                dir[1] = 0; 
                dir[0] = 1;

                for(i = 2; i < ncoeffs; ++i)
                {
                    dir[i] = i;
                }
            }

            // Set up Mapping details
            switch (eid)
            {
            case 0:
                {
                    int cnt = 0;

                    for(i = 0; i < ncoeffs; cnt+=order1-i, ++i)
                    {
                        Map[dir[i]] = cnt; 
                    }
                }
                break;
            case 1:
                Map[dir[0]] = order0;
                Map[dir[1]] = 1;

                for(i = 2; i < ncoeffs; ++i)
                {
                    Map[dir[i]] = order0+1+i; 
                }
                break;
            case 2:
                for(i = 0; i < ncoeffs; ++i)
                {
                    Map[dir[i]] = i; 
                }
                break;
            }
        }

        void StdTriExp::WriteToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  nquad0 = m_base[0]->GetPointsOrder();
            int  nquad1 = m_base[1]->GetPointsOrder();
            const double *z0,*z1,*w0,*w1;

            BasisManagerSingleton::Instance().GetZW(m_base[0],z0,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z1,w1);

            outfile << "Variables = z1,  z2, Coeffs \n" << std::endl;      
            outfile << "Zone, I=" << nquad0 <<", J=" << nquad1 <<", F=Point" << std::endl;

            for(j = 0; j < nquad1; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    outfile << 0.5*(1+z0[i])*(1.0-z1[j])-1 <<  " " << 
                        z1[j] << " " << m_phys[j*nquad0+i] << std::endl;
                }
            }

        }

        //   I/O routine
        void StdTriExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  order0 = m_base[0]->GetBasisOrder();
            int  order1 = m_base[1]->GetBasisOrder();
            int  cnt = 0;
            BstShrDArray wsp  = GetDoubleTmpSpace(order0*order1);

            double *mat = wsp.get(); 

            // put coeffs into matrix and reverse order so that p index is fastest
            // recall q is fastest for tri's

            Vmath::Zero(order0*order1,mat,1);

            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i; ++j,cnt++)
                {
                    mat[i+j*order1] = m_coeffs[cnt];
                }
            }

            outfile <<"Coeffs = [" << " "; 

            for(j = 0; j < order1; ++j)
            {
                for(i = 0; i < order0; ++i)
                {
                    outfile << mat[j*order0+i] <<" ";
                }
                outfile << std::endl; 
            }
            outfile << "]" ; 
        }

    }//end namespace
}//end namespace


/** 
* $Log: StdTriExp.cpp,v $
* Revision 1.2  2006/06/01 14:13:37  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:33  kirby
* *** empty log message ***
*
* Revision 1.50  2006/04/25 20:23:34  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.49  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.48  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.47  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.46  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.45  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.44  2006/03/01 08:25:05  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.43  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/ 


