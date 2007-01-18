///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.cpp
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
// Description: Quadrilateral routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdQuadExp.h>

#ifdef max
#undef max
#endif

namespace Nektar
{
    namespace StdRegions
    {

        StdMatrix StdQuadExp::s_elmtmats;

        StdQuadExp::StdQuadExp()
        {
        }

        StdQuadExp::StdQuadExp(const BasisKey &Ba, const BasisKey &Bb):
        StdExpansion2D(Ba,Bb,Ba.GetBasisOrder()*Bb.GetBasisOrder(),NULL,NULL,true)
        { 
        }

        StdQuadExp::StdQuadExp(const BasisKey &Ba, const BasisKey &Bb,
            double *coeffs, double *phys):
        StdExpansion2D(Ba,Bb,Ba.GetBasisOrder()*Bb.GetBasisOrder(),coeffs,phys,false)
        { 
        }

        StdQuadExp::StdQuadExp(const StdQuadExp &T):
        StdExpansion2D(T)
        {
        }

        StdQuadExp::~StdQuadExp()
        {
        }



        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        double StdQuadExp::Integral(const double *inarray)
        {
            const double *z, *w0, *w1;

            BasisManagerSingleton::Instance().GetZW(m_base[0],z,w0);
            BasisManagerSingleton::Instance().GetZW(m_base[1],z,w1);

            return StdExpansion2D::Integral(inarray,w0,w1);
        }


        void StdQuadExp::IProductWRTBase(const double * inarray, double * outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                inarray,outarray,1);
        }

        void StdQuadExp:: IProductWRTBase(const double *base0, const double *base1,
            const double *inarray, double *outarray,
            int coll_check){
                int i;
                int    nquad0 = m_base[0]->GetPointsOrder();
                int    nquad1 = m_base[1]->GetPointsOrder();
                int    order0 = m_base[0]->GetBasisOrder();
                int    order1 = m_base[1]->GetBasisOrder();
                const double *z,*w0,*w1;
                BstShrDArray tmp  = GetDoubleTmpSpace(nquad0*nquad1);
                BstShrDArray tmp1 = GetDoubleTmpSpace(nquad0*nquad1);

#if FULLDEBUG
                if((m_base[0]->GetAlpha() != 0.0)||(m_base[1]->GetAlpha() != 0.0))
                {
                    ErrorUtil::Error(ErrorUtil::ewarning,"StdQuadExp::IProduct_WRT_B",
                        "Basis has non-zero alpha weight");
                }

                if((m_base[0]->GetBeta() != 0.0)||(m_base[1]->GetBeta() != 0.0))
                {
                    ErrorUtil::Error(ErrorUtil::ewarning,"StdQuadExp::IProduct_WRT_B",
                        "Basis has non-zero beta weight");
                }
#endif

                BasisManagerSingleton::Instance().GetZW(m_base[0],z,w0);
                BasisManagerSingleton::Instance().GetZW(m_base[1],z,w1);

                // Note cannot use outarray as tmp space since dimensions are not always
                // guarenteed to be sufficient 

                // multiply by integration constants 
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Vmul(nquad0,(double*)inarray+i*nquad0,1,(double*)w0,1,
                        tmp.get()+i*nquad0,1);
                }

                for(i = 0; i < nquad0; ++i)
                {
                    Vmath::Vmul(nquad1,tmp.get()+i,nquad0,(double*)w1,1,
                        tmp.get()+i,nquad0);
                }

                if(coll_check&&m_base[0]->Collocation())
                {
                    Vmath::Vcopy(order0*nquad1,tmp.get(),1,tmp1.get(),1);
                }
                else
                {
                    Blas::Dgemm('T','N',order0,nquad1,nquad0,1.0,base0,nquad0,tmp.get(),
                        nquad0,0.0,tmp1.get(),order0);
                }

                if(coll_check&&m_base[1]->Collocation())
                {
                    Vmath::Vcopy(order0*order1,tmp1.get(),1,outarray,1);
                }
                else
                {
                    Blas::Dgemm('N', 'N',order0,order1, nquad1,1.0, tmp1.get(),order0, 
                        base1, nquad1, 0.0, outarray,order0);
                }

            }
	
	void StdQuadExp::FillMode(const int mode, double *outarray)
	{
	    int    i;
	    int   nquad0 = m_base[0]->GetPointsOrder();
	    int   nquad1 = m_base[1]->GetPointsOrder();
	    const double * base0  = m_base[0]->GetBdata();
	    const double * base1  = m_base[1]->GetBdata();
	    int   btmp0 = m_base[0]->GetBasisOrder();
	    int   mode0 = mode%btmp0;
	    int   mode1 = mode/btmp0;
	    
	    
	    ASSERTL2(mode1 == (int)floor((1.0*mode)/btmp0),
		     "Integer Truncation not Equiv to Floor");
	    
	    ASSERTL2(m_ncoeffs <= modes, 
		     "calling argument mode is larger than total expansion order");
	    
	    for(i = 0; i < nquad1; ++i)
	    {
		Vmath::Vcopy(nquad0,(double *)(base0 + mode0*nquad0),1,
			     outarray+i*nquad0,1);
	    }
	    
	    for(i = 0; i < nquad0; ++i)
	    {
		Vmath::Vmul(nquad1,(double *)(base1 + mode1*nquad1),1,
			    outarray+i,nquad0,outarray+i,nquad0);
	    }
	}
	
	void StdQuadExp::GenMassMatrix(double * outarray)
	{
	    int      i;
	    int      order0    = GetBasisOrder(0);
	    int      order1    = GetBasisOrder(1);

	    StdExpansion::GenerateMassMatrix(outarray);
	    
	    // For Fourier basis set the imaginary component of mean mode
	    // to have a unit diagonal component in mass matrix 
	    if(m_base[0]->GetBasisType() == eFourier)
	    {
		for(i = 0; i < order1; ++i)
		{
		    outarray[(order0*i+1)*m_ncoeffs+i*order0+1] = 1.0;
		}
	    }
	    
	    if(m_base[1]->GetBasisType() == eFourier)
	    {
		for(i = 0; i < order0; ++i)
		{
		    outarray[(order0+i)*m_ncoeffs+order0+i] = 1.0;
		}
	    }
	}

	void StdQuadExp::GenLapMatrix(double * outarray)
	{
	    ASSERTL0(false, "Not implemented");
	}
	
	StdMatContainer * StdQuadExp::GetMassMatrix() 
	{
	    StdMatContainer * mat;
	    mat = s_elmtmats.GetLocalMass(this);
	    return mat;
	}

	StdMatContainer * StdQuadExp::GetLapMatrix() 
	{
	    StdMatContainer * mat;
	    mat = s_elmtmats.GetLocalLap(this);
	    return mat;
	}
	

	///////////////////////////////
	/// Differentiation Methods
	///////////////////////////////
	
	void StdQuadExp::Deriv(double *outarray_d0, double *outarray_d1)
	{
	    TensorDeriv(this->m_phys, outarray_d0, outarray_d1);
	}

	void StdQuadExp::Deriv(const double *inarray, double *outarray_d0, 
			       double *outarray_d1)
	{
	    TensorDeriv(inarray, outarray_d0, outarray_d1);
	}
	
	//------------------------------
	// Evaluation Methods
	//-----------------------------
	
	void StdQuadExp::BwdTrans(double * outarray)
	{
	    int           nquad0 = m_base[0]->GetPointsOrder();
	    int           nquad1 = m_base[1]->GetPointsOrder();
	    int           order0 = m_base[0]->GetBasisOrder();
	    int           order1 = m_base[1]->GetBasisOrder();
	    const double *base0  = m_base[0]->GetBdata();
	    const double *base1  = m_base[1]->GetBdata();
	    BstShrDArray tmp  = GetDoubleTmpSpace(nquad0*std::max(order1,nquad1));
	    
	    if(m_base[0]->Collocation())
	    {
		Vmath::Vcopy(nquad0*order1,m_coeffs,1,tmp.get(),1);
	    }
	    else
	    {
		Blas::Dgemm('N','N', nquad0,order1,order0,1.0, base0, nquad0, 
			    m_coeffs, order0,0.0,tmp.get(), nquad0);
	    }
	    
	    if(m_base[1]->Collocation())
	    {
		Vmath::Vcopy(nquad0*nquad1,tmp.get(),1,outarray,1);
	    }
	    else
	    {
		Blas::Dgemm('N','T', nquad0, nquad1,order1, 1.0, tmp.get(), 
			    nquad0, base1, nquad1, 0.0,outarray, nquad0);
	    }    
	}

	void StdQuadExp::FwdTrans(const double *inarray)
	{
	    StdMatContainer *M;
	    
	    if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
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
	double StdQuadExp::Evaluate(const double * coords)
	{
	    return  PhysEvaluate(coords); 
	}



	// For a specified edge 'eid' this function updates a class
	// StdExpMap which contains the mapping of the edge degrees of
	// freedom back into the elemental domain which is also
	// dependent upon the edge orientation. The vertex and edge
	// ordering of the mapping is dependent upon which basis is
	// being considered, i.e. modal expansions the vertices will
	// be first, nodal expansions the vertices will be the two
	// end points
	void StdQuadExp::MapTo(const int edge_ncoeffs, const BasisType Btype, 
			       const int eid, const EdgeOrientation eorient,
			       StdExpMap &Map)
	{
	    
	    int i, start, skip;
	    int *dir, order0,order1;
	    BstShrIArray wsp; 
	    
	    ASSERTL2(eid>=0&&eid <=3,"eid must be between 0 and 3");
	    // make sure have correct memory storage
	    if(edge_ncoeffs != Map.GetLen())
	    {
		Map.SetMapMemory(edge_ncoeffs);
	    }
	    
	    wsp = GetIntTmpSpace(edge_ncoeffs);
	    dir = wsp.get(); 
	    
	    if(eorient == eForwards)
	    {
		for(i = 0; i < edge_ncoeffs; ++i)
		{
		    dir[i] = i;
		}
	    }
	    else
	    {
		if(Btype == eGLL_Lagrange)
		{
		    for(i = 0; i < edge_ncoeffs; ++i)
		    {
			dir[i] = edge_ncoeffs-i-1;
		    }
		}
		else{
		    dir[1] = 0; 
		    dir[0] = 1;
		    for(i = 2; i < edge_ncoeffs; ++i)
		    {
			dir[i] = i;
		    }
		}
	    }
	    
	    order0 = m_base[0]->GetBasisOrder();
	    order1 = m_base[0]->GetBasisOrder();
	    
	    // Set up Mapping details
	    if((eid == 0)||(eid == 2))
	    { 
		ASSERTL2(Btype == m_base[0]->GetBasisType(),
			 "Expansion type of edge and StdQuadExp are different");
		
		switch(Btype)
		{
		case eGLL_Lagrange:
		    ASSERTL2(edge_ncoeffs == order0,
		      "Expansion order of edge and StdQuadExp are different");
		    
		    if(eid == 0)
		    {
			start = 0;
			skip  = 1;
		    }
		    else
		    {
			start = order0*(order1-1);
			skip = 1;
		    }
		    break;
		    
		case eModified_A:
		    if(eid == 0)
		    {
			start = 0;
			skip  = 1;
		    }
		    else
		    {
			start = order0;
			skip = 1;
		    }
		    break;
		default:
		    ASSERTL0(0,"Mapping array is not defined for this expansion");
		    break;
		}
	    }
	    else
	    {
		ASSERTL2(Btype == m_base[1]->GetBasisType(),
			 "Expansion type of edge and StdQuadExp are different");      
		
		switch(Btype)
		{
		case eGLL_Lagrange:
		    ASSERTL2(edge_ncoeffs == order1,
			     "Expansion order of edge and StdQuadExp are different");
		    if(eid == 1)
		    {
			start = order0-1;
			skip  = order0;
		    }
		    else
		    {
			start = 0;
			skip = order0;
		    }
		    break;
		    
		case eModified_A:	
		    if(eid == 1)
		    {
			start = 1;
			skip  = order0;
		    }
		    else
		    {
			start = 0;
			skip = order0;
		    }
		    break;
		default:
		    ASSERTL0(0,"Mapping array is not defined for this expansion");
		    break;
		}
	    }
	    
	    for(i = 0; i < edge_ncoeffs; ++i)
	    {
		Map[dir[i]] = start + i*skip; 
	    }
	    
	}

	// same as MapTo but assume that mapping is provided in modal
	// basis type format where vertices are listed first followed
	// by edges degrees of freedom

	void StdQuadExp::MapTo_ModalFormat(const int edge_ncoeffs, 
					   const BasisType Btype, 
					   const int eid, 
					   const EdgeOrientation eorient,
					   StdExpMap &Map)
	{
	    MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
	    
	    if(Btype == eGLL_Lagrange)
	    {
		int i;
		int vert = Map[edge_ncoeffs-1];
		for(i = edge_ncoeffs-1; i > 1; --i)
		{
		    Map.SetMap(i,Map[i-1]);
		}
		Map.SetMap(1,vert);
	    }
	}
	
	void StdQuadExp::SetInvInfo(StdMatContainer *mat, MatrixType Mform)
	{
	    mat->SetLda(m_ncoeffs);
	    mat->SetMatForm(eSymmetric_Positive);
      
	    if(GeoFacType() == eRegular)
	    {
		switch(Mform)
		{
		case eMassMatrix:
		    {
			switch(m_base[1]->GetBasisType())
			{
			case eOrtho_A: case eLegendre:
			    if(m_base[1]->ExactIprodInt())
			    {
				goto eQuadOrtho1;
			    }
			    break;
			case eGLL_Lagrange:
			    if(m_base[1]->Collocation())
			    {
				goto eQuadOrtho1;
			    }
			    break;
			case eFourier:
			    goto eQuadOrtho1;
			    break;
			default:
			    mat->SetMatForm(eSymmetric_Positive_Banded);
			    mat->SetBwidth(m_base[0]->GetBasisOrder()*m_base[1]->GetBasisOrder());
			    break;
			eQuadOrtho1:
			    {
				switch(m_base[0]->GetBasisType())
				{
				case eOrtho_A: case eLegendre:
				    if(m_base[0]->ExactIprodInt())
				    {
					goto eQuadOrtho2;
				    }
				    break;
				case eGLL_Lagrange:
				    if(m_base[0]->Collocation())
				    {
					goto eQuadOrtho2;
				    }
				    break;
				case eFourier:
				    goto eQuadOrtho2;
				    break;
				default:
				    mat->SetMatForm(eSymmetric_Positive_Banded);
				    mat->SetBwidth(m_base[0]->GetBasisOrder());
				    break;
				eQuadOrtho2:
				    {
					mat->SetMatForm(eSymmetric_Positive_Banded);
					mat->SetBwidth(1);
					break;
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
  
    } //end namespace			
}//end namespace

/** 
* $Log: StdQuadExp.cpp,v $
* Revision 1.7  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.6  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/07/02 17:16:18  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.2  2006/06/01 14:13:36  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.38  2006/04/25 20:23:34  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.37  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.36  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.35  2006/03/05 22:11:03  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.34  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.33  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
**/ 





