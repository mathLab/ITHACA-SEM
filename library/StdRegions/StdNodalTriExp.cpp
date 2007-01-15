///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.cpp
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
// Description: Nodal triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTriExp.h>

namespace Nektar
{
    namespace StdRegions
    {

        StdMatrix StdNodalTriExp::s_elmtmats;

        StdNodalTriExp::StdNodalTriExp(const BasisKey &Ba, const BasisKey &Bb, 
				       NodalBasisType Ntype):
	    StdTriExp(Ba,Bb)
        {
            m_nbtype = Ntype;

            ASSERTL0(m_base[0]->GetBasisOrder() == m_base[1]->GetBasisOrder(),
                "Nodal basis initiated with different orders in the a "
		     "and b directions");
        }

        StdNodalTriExp::StdNodalTriExp(const BasisKey &Ba,  const BasisKey &Bb,
		       NodalBasisType Ntype, double *coeffs,  double *phys):
	    StdTriExp(Ba,Bb,coeffs,phys)
        {    
	    m_nbtype = Ntype;

            ASSERTL0(m_base[0]->GetBasisOrder() == m_base[1]->GetBasisOrder(),
		     "Nodal basis initiated with different "
		     "orders in the a and b directions");
        }
	
        StdNodalTriExp::StdNodalTriExp(const StdNodalTriExp &T):
	    StdTriExp(T)
        {
	    m_nbtype = T.m_nbtype;
        }

        // Destructor
        StdNodalTriExp::~StdNodalTriExp()
        { 
        }

        void StdNodalTriExp::GenNBasisTransMatrix(double * outarray)
        {
            int             i,j;
            int             tot_order = GetNcoeffs();
            const double*   r; 
            const double*   s; 
            const double*   t;
            double          c[2];

            NBasisManagerSingleton::Instance().GetNodePoints(m_nbtype,
                m_base[0]->GetBasisOrder(),r,s,t);

            for(i = 0; i < tot_order; ++i)
            {
                // fill physical space with mode i
                StdTriExp::FillMode(i,m_phys);

                // interpolate mode i to the Nodal points 'j' and store in outarray
                for(j = 0; j < tot_order; ++j)
                {
                    c[0] = r[j];
                    c[1] = s[j];
                    // define matrix in row major format to have rows of 
                    // all the different expansion bases defined at the nodal point 
                    outarray[j*tot_order+i] = StdTriExp::Evaluate(c);
                }
            }
        }

        void StdNodalTriExp::NodalToModal()
        {
            NodalToModal(m_coeffs); 
        }

        void StdNodalTriExp::NodalToModal(double *in_out_array)
        {
            StdMatContainer *M;

            M = GetNBasisTransMatrix();
            M->Solve(in_out_array,1);
        }


        void StdNodalTriExp::NodalToModalTranspose()
        {
            NodalToModal(m_coeffs); 
        }

	// Operate with transpose of NodalToModal transformation
        void StdNodalTriExp::NodalToModalTranspose(double *in_out_array)
        {
            StdMatContainer *M;

            M = GetNBasisTransMatrix();
            M->SolveTranspose(in_out_array,1);
        }


        StdMatContainer * StdNodalTriExp::GetNBasisTransMatrix() 
        {
            StdMatContainer * mat;
            mat = s_elmtmats.GetNBasisTrans(this);
            return mat;
        }

        void StdNodalTriExp::ModalToNodal()
        {
	    ModalToNodal(m_coeffs);
        }

        void StdNodalTriExp::ModalToNodal(double *in_out_array)
        {
            StdMatContainer *M;
            M = GetNBasisTransMatrix();
            M->Mxv(in_out_array,in_out_array);
        }


        //////////////////////////////
        /// Integration Methods
        //////////////////////////////


        void StdNodalTriExp::IProductWRTBase(const double * inarray, 
            double * outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                inarray,outarray);
        }

        /** \brief Calculate the inner product of inarray with respect to
        the basis B=base0[p]*base1[pq] and put into outarray:

        This function uses the StdTriExp routine and then 
        calls ModalToNodal to transform to Nodal basis
        **/

        void StdNodalTriExp:: IProductWRTBase(const double *base0, 
            const double *base1, 
            const double *inarray,
            double *outarray)
        {
            // Take inner product with respect to Orthgonal basis using
            // StdTri routine

            StdTriExp::IProductWRTBase(base0,base1,inarray,outarray);

	    NodalToModalTranspose(outarray);
        }


        /** Fill outarray with nodal mode 'mode' of expansion and put in m_phys
        */

        void StdNodalTriExp::FillMode(const int mode, double *outarray)
        {

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs,m_coeffs,1);

            m_coeffs[mode] = 1.0;
	    BwdTrans(outarray);
        }

        StdMatContainer * StdNodalTriExp::GetMassMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalMass(this);
            return tmp;
        }

        StdMatContainer * StdNodalTriExp::GetLapMatrix() 
        {
            StdMatContainer * tmp;
            tmp = s_elmtmats.GetLocalLap(this);
            return tmp;
        }

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        void StdNodalTriExp::Deriv(const double *inarray, double *outarray_d0, 
            double *outarray_d1)
        {
            StdTriExp::Deriv(inarray, outarray_d0, outarray_d1);
        }

        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////

        // Currently convert nodal values into tranformed values and
        // backward transform

        void StdNodalTriExp::BwdTrans(double * outarray)
        {
            BstShrDArray tmp  = GetDoubleTmpSpace(m_ncoeffs);

            // save nodal values
            Blas::Dcopy(m_ncoeffs,m_coeffs,1,tmp.get(),1);
            NodalToModal();
            StdTriExp::BwdTrans(outarray);
            Blas::Dcopy(m_ncoeffs,tmp.get(),1,m_coeffs,1);
        }

        void StdNodalTriExp::FwdTrans(const double * inarray)
        {
            StdMatContainer *M;
	    
            IProductWRTBase(inarray,m_coeffs);
            M = GetMassMatrix();
	    M->ShowMatrixStructure(stdout);
            M->Solve(m_coeffs,1);
        }

        double StdNodalTriExp::Evaluate(const double * coords)
        {
            return StdTriExp::Evaluate(coords);
        }


        void  StdNodalTriExp::MapTo(const int edge_ncoeffs, 
				    const BasisType Btype,
				    const int eid, 
				    const EdgeOrientation eorient, 
				    StdExpMap &Map)
        {
	    
            int i;
            int *dir, order0,order1;
            BstShrIArray wsp; 

            ASSERTL2(eid>=0&&eid<=2,"eid must be between 0 and 2");

            ASSERTL2(Btype == m_base[0]->GetBasisType(),
                "Expansion type of edge and StdQuadExp are different");

            // make sure have correct memory storage
            if(edge_ncoeffs != Map.GetLen())
            {
                Map.SetMapMemory(edge_ncoeffs);
            }

            order0 = m_base[0]->GetBasisOrder();
            order1 = m_base[1]->GetBasisOrder();

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
                dir[1] = 0; 
                dir[0] = 1;
		
                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    dir[i] = edge_ncoeffs-i+1;
                }
            }

            // Set up Mapping details
            switch (eid)
            {
            case 0:
		Map[dir[0]] = 0;
		Map[dir[1]] = 1;
		    
		
		for(i = 2; i < edge_ncoeffs;  ++i)
		{
		    Map[dir[i]] = i+1; 
		}
                break;
            case 1:
		Map[dir[0]] = 1;
		Map[dir[1]] = 2;
		
                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = i+edge_ncoeffs-1; 
                }
                break;
            case 2:
		Map[dir[0]] = 0;
		Map[dir[1]] = 2;

                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = i+2*edge_ncoeffs-3; 
                }
                break;
            }
        }

	void StdNodalTriExp::MapTo_ModalFormat(const int edge_ncoeffs, 
						const BasisType Btype, 
						const int eid, 
						const EdgeOrientation eorient,
						StdExpMap &Map)
	{	    
	    MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
	}

	
	void StdNodalTriExp::SetInvInfo(StdMatContainer *mat, MatrixType Mform)
	{
	    mat->SetLda(m_ncoeffs);
	    mat->SetMatForm(eSymmetric_Positive);
	    
	    if(GeoFacType() == eRegular)
	    {
		switch(Mform)
		{
		case eMassMatrix:
		    // Nothing additional to be done 
		    break;
		case eLapMatrix:
		    mat->SetMatForm(eSymmetric);	
		    break;
                case eNBasisTrans:
                    mat->SetMatForm(eGeneral_Full);	
                    break;
		default:
		    ASSERTL0(false, "MatrixType not known");
		    break;
		    
		}
	    }
	}
	
    } // end of namespace
} // end of namespace

/** 
* $Log: StdNodalTriExp.cpp,v $
* Revision 1.4  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.3  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.19  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.18  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.17  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.16  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
**/ 

