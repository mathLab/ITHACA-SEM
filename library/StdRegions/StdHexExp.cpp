///////////////////////////////////////////////////////////////////////////////
//
// File StdexExp.cpp
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
// Description: Heaxhedral methods
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdHexExp.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdMatrix StdHexExp::s_elmtmats;


        StdHexExp::StdHexExp(const BasisKey &Ba, const BasisKey &Bb, 
            const BasisKey &Bc):
        StdExpansion3D(Ba,Bb,Bc,Ba.GetBasisOrder()*Bb.GetBasisOrder()*Bc.GetBasisOrder(),
            NULL,NULL,true)
        {    
        }

        StdHexExp::StdHexExp(const BasisKey &Ba, const BasisKey &Bb, 
            const BasisKey &Bc, double *coeffs, double *phys):
        StdExpansion3D(Ba,Bb,Bc,Ba.GetBasisOrder()*Bb.GetBasisOrder()*Bc.
            GetBasisOrder(),coeffs,phys,false)
        {
        }

        StdHexExp::StdHexExp(const StdHexExp &T):
        StdExpansion3D(T)
        {
        }

        // Destructor
        StdHexExp::~StdHexExp()
        {   
        }  

        //////////////////////////////
        /// Integration Methods
        //////////////////////////////

        double StdHexExp::Integral(const double *inarray)
        {
            int    i,j;
            double Int = 0.0;
            const double *z,*w0,*w1,*w2;
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            int    nquad2 = m_base[2]->GetPointsOrder();
            BstShrDArray tmp = GetDoubleTmpSpace(nquad0*nquad1*nquad2);

            BasisManagerSingleton::Instance().GetZW(m_base[0],z,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z,w1);
            BasisManagerSingleton::Instance().GetZW(m_base[2],z,w2);

            // multiply by integration constants

            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0,(double*)inarray+i*nquad0,1,(double*)w0,1,
                    tmp.get()+i*nquad0,1);
            }

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,tmp.get()+i+j*nquad0*nquad1,nquad0,(double*)w1,1,
                        tmp.get()+i+j*nquad0*nquad1, nquad0);
                }
            }

            for(i = 0; i < nquad2; ++i)
            {
                Blas::Dscal(nquad0*nquad1,(double) w2[i],tmp.get()+i*nquad0*nquad1,1);
            }

            Int = Vmath::Vsum(nquad0*nquad1*nquad2,tmp.get(),1);

            return Int;
        }


        /**
        f_pq = (phi_p phi_q, u)

        **/

        void StdHexExp::IProductWRTBase(const double * inarray, double * outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                m_base[2]->GetBdata(),inarray,outarray,1);
        }



        void StdHexExp:: IProductWRTBase(const double *base0, const double *base1, 
            const double *base2, const double *inarray,
            double *outarray, int coll_check)
        {
            int i,j;
            int    nquad0 = m_base[0]->GetPointsOrder();
            int    nquad1 = m_base[1]->GetPointsOrder();
            int    nquad2 = m_base[2]->GetPointsOrder();
            int    order0 = m_base[0]->GetBasisOrder();
            int    order1 = m_base[1]->GetBasisOrder();
            int    order2 = m_base[2]->GetBasisOrder();

            const double *z,*w0,*w1,*w2;

            int size = std::max(order0,nquad0)*std::max(order1,nquad1)*
                std::max(order2,nquad2);
            BstShrDArray tmp  = GetDoubleTmpSpace(size);
            BstShrDArray tmp1 = GetDoubleTmpSpace(size);

#if FULLDEBUG
            if((m_base[0]->GetAlpha() != 0.0)||(m_base[1]->GetAlpha() != 0.0)||
                m_base[2]->GetAlpha() != 0.0)
            {
                ErrorUtil::Error(ErrorUtil::ewarning,"StdHexExp::IProductWRTBase",
                    "Basis has non-zero alpha weight");
            }

            if((m_base[0]->GetBeta() != 0.0)||(m_base[1]->GetBeta() != 0.0)||
                m_base[2]->GetBeta() != 0.0)
            {
                ASSERTL0(false, "Basis has non-zero beta weight");
            }
#endif

            BasisManagerSingleton::Instance().GetZW(m_base[0],z,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z,w1);
            BasisManagerSingleton::Instance().GetZW(m_base[2],z,w2);

            // Note cannot use outarray as tmp space since dimensions are not always
            // guarenteed to be sufficient 

            // multiply by integration constants 

            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0,(double*)inarray+i*nquad0,1,(double*)w0,1,
                    tmp.get()+i*nquad0,1);
            }

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,tmp.get()+i+j*nquad0*nquad1,nquad0,(double*)w1,1,
                        tmp.get()+i+j*nquad0*nquad1, nquad0);
                }
            }

            for(i = 0; i < nquad2; ++i)
            {
                Blas::Dscal(nquad0*nquad1,(double) w2[i],tmp.get()+i*nquad0*nquad1,1);
            }


            if(coll_check&&m_base[0]->Collocation())
            {
                Vmath::Vcopy(order0*nquad1*nquad2,tmp.get(),1,tmp1.get(),1);
            }
            else
            {
                Blas::Dgemm('T','N',order0,nquad1*nquad2,nquad0,1.0,base0,nquad0,
                    tmp.get(),nquad0,0.0,tmp1.get(),order0);
            }


            if(coll_check&&m_base[1]->Collocation())
            {
                Vmath::Vcopy(order0*order1*nquad2,tmp.get(),1,tmp.get(),1);
            }
            else
            {
                for(i=0;i<nquad2;++i)
                {
                    Blas::Dgemm('N','N',order0,order1,nquad1,1.0,tmp1.get()+i*order0*nquad1,
                        order0,base1, nquad1, 0.0, tmp.get()+i*order0*order1,order0);
                }
            }

            if(coll_check&&m_base[2]->Collocation())
            {
                Vmath::Vcopy(order0*order1*order2,tmp.get(),1,outarray,1);
            }
            else
            {
                Blas::Dgemm('N','N',order0*order1,order2,nquad2,1.0,tmp.get(), 
                    order0*order1, base2, nquad2, 0.0, outarray, 
                    order0*order1);
            }
        }


        /** Fill outarray with mode 'mode' of expansion

        Note for hexahedral expansions _base[0] (i.e. p)  modes run fastest

        */

        void StdHexExp::FillMode(const int mode, double *outarray)
        {
            int    i,j;
            int   nquad0 = m_base[0]->GetPointsOrder();
            int   nquad1 = m_base[1]->GetPointsOrder();
            int   nquad2 = m_base[2]->GetPointsOrder();
            const double * base0  = m_base[0]->GetBdata();
            const double * base1  = m_base[1]->GetBdata();
            const double * base2  = m_base[2]->GetBdata();
            int   btmp0 = m_base[0]->GetBasisOrder();
            int   btmp1 = m_base[1]->GetBasisOrder();
            int   mode2 = mode/(btmp0*btmp1);
            int   mode1 = (mode-mode2*btmp0*btmp1)/btmp0;
            int   mode0 = (mode-mode2*btmp0*btmp1)%btmp0;

            ASSERTL2(mode2 == (int)floor((1.0*mode)/(btmp0*btmp1)),
                "Integer Truncation not Equiv to Floor");
            ASSERTL2(mode1 == (int)floor((1.0*mode-mode2*btmp0*btmp1)/(btmp0*btmp1)),
                "Integer Truncation not Equiv to Floor");
            ASSERTL2(m_ncoeffs <= modes,
                "calling argument mode is larger than total expansion order");


            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vcopy(nquad0,(double *)(base0 + mode0*nquad0),1,
                    outarray+i*nquad0,1);
            }

            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,(double *)(base1 + mode1*nquad1),1,
                        outarray+i+j*nquad0*nquad1,nquad0,
                        outarray+i+j*nquad0*nquad1,nquad0);
                }
            }

            for(i = 0; i < nquad2; i++)
            {
                Blas::Dscal(nquad0*nquad1,base2[mode2*nquad2+i],
                    outarray+i*nquad0*nquad1,1);
            }
        }

        void StdHexExp::GenMassMatrix(double * outarray)
        {
            int      i,j;
            int      order0    = GetBasisOrder(0);
            int      order1    = GetBasisOrder(1);
            int      order2    = GetBasisOrder(2);
            int      tot_order = GetNcoeffs();

            StdExpansion::GenerateMassMatrix(outarray);


            // For Fourier basis set the imaginary component of mean mode
            // to have a unit diagonal component in mass matrix 
            if(m_base[0]->GetBasisType() == eFourier)
            {
                for(i = 0; i < order1*order2; ++i)
                {
                    outarray[(order0*i+1)*tot_order+i*order0+1] = 1.0;
                }
            }

            if(m_base[1]->GetBasisType() == eFourier)
            {
                for(j = 0; j < order2; ++j)
                {
                    for(i = 0; i < order0; ++i)
                    {
                        outarray[(order0+i)*tot_order+order0+i + 
                            j*(order0*order1)*(tot_order+1)] = 1.0;
                    }
                }
            }

            if(m_base[2]->GetBasisType() == eFourier)
            {
                for(i = 0; i < order0*order1; ++i)
                {
                    outarray[(order0*order1)*(tot_order+1)+i*tot_order +i] = 1.0;
                }
            }
        }

        void StdHexExp::GenLapMatrix(double * outarray)
        {
            ASSERTL0(false, "Not implemented");
        }

        StdMatContainer * StdHexExp::GetMassMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalMass(this);
            return tmp;
        }

        StdMatContainer * StdHexExp::GetLapMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalLap(this);
            return tmp;
        }


        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        void StdHexExp::Deriv(double *outarray_d0, double *outarray_d1, 
            double *outarray_d2)
        {
            TensorDeriv(this->m_phys, outarray_d0, outarray_d1, outarray_d2);
        }

        void StdHexExp::Deriv(const double *inarray, double *outarray_d0, 
            double *outarray_d1, double *outarray_d2)
        {
            TensorDeriv(inarray, outarray_d0, outarray_d1, outarray_d2);
        }

        //------------------------------
        // Evaluation Methods
        //-----------------------------

        void StdHexExp::BwdTrans(double * outarray)
        {
            int           i;
            int           nquad0 = m_base[0]->GetPointsOrder();
            int           nquad1 = m_base[1]->GetPointsOrder();
            int           nquad2 = m_base[2]->GetPointsOrder();

            int           order0 = m_base[0]->GetBasisOrder();
            int           order1 = m_base[1]->GetBasisOrder();
            int           order2 = m_base[2]->GetBasisOrder();

            const double *base0  = m_base[0]->GetBdata();
            const double *base1  = m_base[1]->GetBdata();
            const double *base2  = m_base[2]->GetBdata();

            double *tmp  = new double [nquad0*std::max(nquad1,order1)*order2];
            double *tmp1 = new double [nquad0*nquad1*std::max(nquad2,order2)];

            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(nquad0*order1*order2,m_coeffs,1,tmp,1);
            }
            else
            {
                Blas::Dgemm('N','N', nquad0,order1*order2,order0,1.0, base0, nquad0,
                    m_coeffs, order0,0.0,tmp,nquad0);
            }


            if(m_base[1]->Collocation())
            {
                Vmath::Vcopy(nquad0*nquad1*order2,tmp,1,tmp1,1);
            }
            else
            {
                for(i = 0; i < order2; ++i)
                {
                    Blas::Dgemm('N','T',nquad0,nquad1,order1,1.0,tmp+i*nquad0*order1,
                        nquad0,base1,nquad1,0.0,tmp1+i*nquad0*nquad1,nquad0);
                }
            }

            if(m_base[2]->Collocation())
            {
                Vmath::Vcopy(nquad0*nquad1*nquad2,tmp1,1,outarray,1);
            }
            else
            {
                Blas::Dgemm('N','T',nquad0*nquad1,nquad2,order2,1.0,tmp1,
                    nquad0*nquad1,base2,nquad2,0.0,outarray,nquad0*nquad1);
            }
        }

        void StdHexExp::FwdTrans(const double *inarray)
        {
            StdMatContainer *M;

            if((m_base[0]->Collocation())&&(m_base[1]->Collocation())&&
                (m_base[2]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(),inarray,1,m_coeffs,1);
            }
            else
            {
                IProductWRTBase(inarray,m_coeffs);
                M = GetMassMatrix();
                M->Solve(m_coeffs,1);
            }
        }

        /// Single Point Evaluation
        double StdHexExp::Evaluate(const double * coords)
        {
            return  PhysEvaluate(coords); 
        }

	void StdHexExp::SetInvInfo(StdMatContainer *mat, MatrixType Mform)
	{
	    mat->SetLda(m_ncoeffs);
	    mat->SetMatForm(eSymmetric_Positive);
      
	    if(GeoFacType() == eRegular)
	    {
		switch(Mform)
		{
		case eMassMatrix:
		    switch(m_base[2]->GetBasisType())
		    {
		    case eOrtho_A: case eLegendre:
			if(m_base[2]->ExactIprodInt())
			{
			    goto eHexOrtho1;
			}
			break;
		    case eGLL_Lagrange:
			if(m_base[1]->Collocation())
			{
			    goto eHexOrtho1;
			}
			break;
		    case eFourier:
			goto eHexOrtho1;
			break;
		    eHexOrtho1:
			{
			    switch(m_base[1]->GetBasisType())
			    {
			    case eOrtho_A: case eLegendre:
				if(m_base[1]->ExactIprodInt())
				{
				    goto eHexOrtho2;
				}
				break;
			    case eGLL_Lagrange:
				if(m_base[1]->Collocation())
				{
				    goto eHexOrtho2;
                                        }
				break;
			    case eFourier:
				goto eHexOrtho2;
				break;
			    default:
				mat->SetMatForm(eSymmetric_Positive_Banded);
				mat->SetBwidth(m_base[0]->GetBasisOrder()*m_base[1]->GetBasisOrder());
				break;
			    eHexOrtho2:
				{
				    switch(m_base[0]->GetBasisType())
				    {
				    case eOrtho_A: case eLegendre:
					if(m_base[0]->ExactIprodInt())
					{
					    goto eHexOrtho3;
					}
					break;
				    case eGLL_Lagrange:
					if(m_base[0]->Collocation())
					{
					    goto eHexOrtho3;
					}
					break;
				    case eFourier:
					goto eHexOrtho3;
					break;
				    default:
					mat->SetMatForm(eSymmetric_Positive_Banded);
					mat->SetBwidth(m_base[0]->GetBasisOrder());
					break;
				    eHexOrtho3:
					{
					    mat->SetMatForm(eSymmetric_Positive_Banded);
					    mat->SetBwidth(1);
					    break;
					}
				    }
				}
			    }
			}
		    }
		    break;
		case eLapMatrix:
		    mat->SetMatForm(eSymmetric);	
		    break;
		default:
		    ASSERTL0(false, "MatrixType not known");
		    break;
	    
		}
	    }
	}
  
    }//end namespace
}//end namespace

/** 
* $Log: StdHexExp.cpp,v $
* Revision 1.3  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.2  2006/06/01 14:13:36  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.23  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.22  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.21  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.20  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.19  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.18  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
*
**/ 


