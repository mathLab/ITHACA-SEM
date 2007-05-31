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
	
        StdTriExp::StdTriExp() // default constructor of StdExpansion is directly called. 
        {
        } //default constructor


        StdTriExp::StdTriExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb):
	    StdExpansion2D(Ba.GetNumModes()*(Ba.GetNumModes()+1)/2+
			   Ba.GetNumModes()*(Bb.GetNumModes()-Ba.GetNumModes()), 
			   Ba,Bb)
        {    

            if(Ba.GetNumModes() >  Bb.GetNumModes())
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
        NekDouble StdTriExp::Integral(const ConstArray<OneD, NekDouble>& inarray)
        {
            ConstArray<OneD, NekDouble> w0,z1,w1;
            int    i,nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> w1_tmp = Array<OneD, NekDouble>(nquad1);

	    w0 = ExpPointsProperties(0)->GetW();
	    ExpPointsProperties(1)->GetZW(z1,w1);

            switch(m_base[1]->GetPointsType())
	    {
	    case LibUtilities::eGaussLobattoLegendre: // Legendre inner product 
		for(i = 0; i < nquad1; ++i)
                {
                    w1_tmp[i] = 0.5*(1-z1[i])*w1[i];
                }
                break;
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (0,1) Jacobi Inner product 
                Vmath::Smul(nquad1,0.5,(NekDouble *)w1.get(),1,w1_tmp.get(),1);      
                break;
            }
	    
            return StdExpansion2D::Integral(inarray,w0,w1_tmp);
        }


	void StdTriExp::IProductWRTBase(const ConstArray<OneD, NekDouble>& inarray, 
					Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),inarray,
                outarray);
        }
    
	void StdTriExp:: IProductWRTBase(const ConstArray<OneD, NekDouble>& base0, 
					 const ConstArray<OneD, NekDouble>& base1, 
					 const ConstArray<OneD, NekDouble>& inarray, 
					 Array<OneD, NekDouble> & outarray)
        {
            int    i,mode;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
            int    order1 = m_base[1]->GetNumModes();
            ConstArray<OneD, NekDouble> z1,w0,w1;
            Array<OneD, NekDouble> tmp  = Array<OneD, NekDouble>(nquad0*nquad1);
            Array<OneD, NekDouble> tmp1 = Array<OneD, NekDouble>(nquad0*nquad1);

	    w0 = ExpPointsProperties(0)->GetW();
	    ExpPointsProperties(1)->GetZW(z1,w1);

            ASSERTL2((m_base[1]->GetBasisType() == LibUtilities::eOrtho_B)||
                (m_base[1]->GetBasisType() == LibUtilities::eModified_B), 
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
                Vmath::Vmul(nquad0,(NekDouble*)&inarray[0]+i*nquad0,1,
			    w0.get(),1, &tmp[0]+i*nquad0,1);
            }

            switch(m_base[1]->GetPointsType())
            {
            case LibUtilities::eGaussLobattoLegendre: // Legendre inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i], &tmp[0]+i*nquad0,1);
                }
                break;
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*w1[i], &tmp[0]+i*nquad0,1);      
                }
                break;
            }
	    
            // Inner product with respect to 'a' direction 
            Blas::Dgemm('T','N',nquad1,order0,nquad0,1.0,&tmp[0],nquad0,
			base0.get(),nquad0,0.0,&tmp1[0],nquad1);

            // Inner product with respect to 'b' direction 
            for(mode=i=0; i < order0; ++i)
            {
                Blas::Dgemv('T',nquad1,order1-i,1.0, base1.get()+mode*nquad1,
			    nquad1,&tmp1[0]+i*nquad1,1, 0.0, 
			    &outarray[0] + mode,1);
                mode += order1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                outarray[1] += Blas::Ddot(nquad1,base1.get()+nquad1,1,
					  &tmp1[0]+nquad1,1);
            }
        }
    
	void StdTriExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int    i,m;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   order0 = m_base[0]->GetNumModes();
            int   order1 = m_base[1]->GetNumModes();
            ConstArray<OneD, NekDouble> base0 = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> base1 = m_base[1]->GetBdata();
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
            if((mode == 1)&&(m_base[0]->GetBasisType() == LibUtilities::eModified_A))
            {
                Vmath::Fill(nquad0*nquad1,1.0,&outarray[0],1);
            }
            else
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad0,(NekDouble *)(base0.get()+mode0*nquad0),
				 1,&outarray[0]+i*nquad0,1);
                }
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode*nquad1),1,&outarray[0]+i,
                    nquad0,&outarray[0]+i,nquad0);
            }
        }
    
	//-----------------------------
	// Differentiation Methods
	//-----------------------------
	
	void StdTriExp::PhysDeriv(const ConstArray<OneD, NekDouble>& inarray, 
				  Array<OneD, NekDouble> &out_d0, 
				  Array<OneD, NekDouble> &out_d1,
                                  Array<OneD, NekDouble> &out_d2)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD, NekDouble> z0,z1;
            Array<OneD, NekDouble> d0;
            Array<OneD, NekDouble> wsp1  = Array<OneD, NekDouble>(nquad0*nquad1);
            NekDouble *gfac = wsp1.get();

	    z0 = ExpPointsProperties(0)->GetZ();
	    z1 = ExpPointsProperties(1)->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac[i] = 2.0/(1-z1[i]);
            }
	    
            if(out_d1.num_elements() > 0)// if no d1 required do not need to calculate both deriv
            {
                PhysTensorDeriv(inarray, out_d0, out_d1);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],&out_d0[0]+i*nquad0,1);
                }
            }
            else
            {
                if(out_d0.num_elements() > 0)// need other local callopsed derivative for d1 
                {
                    d0 = Array<OneD, NekDouble>(nquad0*nquad1);
                }

                PhysTensorDeriv(inarray, d0, out_d1);

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac[i],&d0[0]+i*nquad0,1);
                }

                // set up geometric factor: (1_z0)/(1-z1)
                for(i = 0; i < nquad0; ++i)
                {
                    gfac[i] = 0.5*(1+z0[i]);
                }

                for(i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,gfac,1,&d0[0]+i*nquad0,1,&out_d1[0]+i*nquad0,1,
				 &out_d1[0]+i*nquad0,1);
                }	
            }
        }
	

        ///////////////////////////////
        // Evaluation Methods
        ///////////////////////////////

        void StdTriExp::BwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
				 Array<OneD, NekDouble> &outarray)
        {
            int           i,mode;
            int           nquad0 = m_base[0]->GetNumPoints();
            int           nquad1 = m_base[1]->GetNumPoints();
            int           order0 = m_base[0]->GetNumModes();
            int           order1 = m_base[1]->GetNumModes();
            ConstArray<OneD, NekDouble> base0  = m_base[0]->GetBdata();
            ConstArray<OneD, NekDouble> base1  = m_base[1]->GetBdata();
            Array<OneD, NekDouble> tmp  = Array<OneD, NekDouble>(order0*nquad1);


            ASSERTL2((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)||
		     (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
		     "Basis[1] is not of general tensor type");

            for(i = mode = 0; i < order0; ++i)
            {
                Blas::Dgemv('N', nquad1,order1-i,1.0,base1.get()+mode*nquad1,
			    nquad1,&inarray[0]+mode,1,0.0,&tmp[0]+i*nquad1,1);
                mode += order1-i;
            }

            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                Blas::Daxpy(nquad1,inarray[1],base1.get()+nquad1,1,
			    &tmp[0]+nquad1,1);
            }

            Blas::Dgemm('N','T', nquad0,nquad1,order0,1.0, base0.get(),nquad0, 
			&tmp[0], nquad1,0.0, &outarray[0], nquad0);
        }

        void StdTriExp::FwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
				 Array<OneD, NekDouble> &outarray)
        {
	    
            IProductWRTBase(inarray,outarray);

	    // get Mass matrix
	    StdLinSysKey         masskey(eMassMatrix,DetShapeType(),*this);
	    DNekLinSysSharedPtr  matsys = m_stdLinSysManager[masskey];
	    
	    // solve inverse of system
	    DNekVec   v(m_ncoeffs,outarray,eWrapper);
	    matsys->Solve(v,v);
	}

        NekDouble StdTriExp::PhysEvaluate(const ConstArray<OneD, NekDouble>& coords)
        {
            Array<OneD, NekDouble> coll = Array<OneD, NekDouble>(2);

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

            return  StdExpansion2D::PhysEvaluate2D(coll); 
        }
	

        void  StdTriExp::MapTo(const int edge_ncoeffs, 
			       const LibUtilities::BasisType Btype,
			       const int eid, 
			       const EdgeOrientation eorient, 
			       StdExpMap &Map)
        {
	    
            int i;
            int *dir, order0,order1;
            Array<OneD, int> wsp; 

            ASSERTL2(eid>=0&&eid<=2,"eid must be between 0 and 2");
            ASSERTL2(Btype == LibUtilities::eModified_A,"Mapping only set up "
                "for Modified_A edges");
            ASSERTL2(Btype == m_base[0]->GetBasisType(),
                "Expansion type of edge and StdQuadExp are different");

            // make sure haved correct memory storage
            if(edge_ncoeffs != Map.GetLen())
            {
                Map.SetMapMemory(edge_ncoeffs);
            }

            order0 = m_base[0]->GetNumModes();
            order1 = m_base[1]->GetNumModes();

            wsp = Array<OneD, int>(edge_ncoeffs);
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
                    dir[i] = i;
                }
            }

            // Set up Mapping details
            switch (eid)
            {
            case 0:
		{
                    int cnt = 0;
		    
                    for(i = 0; i < edge_ncoeffs; cnt+=order1-i, ++i)
                    {
                        Map[dir[i]] = cnt; 
                    }
                }
                break;
            case 1:
                Map[dir[0]] = order1;
                Map[dir[1]] = 1;
		
                for(i = 2; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = order1+i-1; 
                }
                break;
            case 2:
                for(i = 0; i < edge_ncoeffs; ++i)
                {
                    Map[dir[i]] = i; 
                }
                break;
            }
        }

	// currently same as MapTo 

	void StdTriExp::MapTo_ModalFormat(const int edge_ncoeffs, 
					   const LibUtilities::BasisType Btype, 
					   const int eid, 
					   const EdgeOrientation eorient,
					   StdExpMap &Map)
	{	    
	    MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
	}


        void StdTriExp::WriteToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD, NekDouble> z0,z1;

	    z0 = ExpPointsProperties(0)->GetZ();
	    z1 = ExpPointsProperties(1)->GetZ();

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
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  cnt = 0;
            Array<OneD, NekDouble> wsp  = Array<OneD, NekDouble>(order0*order1);

            NekDouble *mat = wsp.get(); 

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

	void StdTriExp::GetCoords(Array<OneD, NekDouble> &coords_0, 
				  Array<OneD, NekDouble> &coords_1)
	{
	    ConstArray<OneD, NekDouble> z0 = ExpPointsProperties(0)->GetZ();
	    ConstArray<OneD, NekDouble> z1 = ExpPointsProperties(1)->GetZ();
	    int nq0 = GetNumPoints(0);
	    int nq1 = GetNumPoints(1);
	    int i,j;

	    for(i = 0; i < nq1; ++i)
	    {
		for(j = 0; j < nq0; ++j)
		{
		    coords_0[i*nq0+j] = (1+z0[j])*(1-z1[i])/2.0 - 1.0;
		}
		Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
	    }
	}

    }//end namespace
}//end namespace


/** 
* $Log: StdTriExp.cpp,v $
* Revision 1.16  2007/05/22 02:01:50  bnelson
* Changed Array::size to Array::num_elements.
*
* Fixed some compiler errors in assertions.
*
* Revision 1.15  2007/05/15 05:18:24  bnelson
* Updated to use the new Array object.
*
* Revision 1.14  2007/04/10 14:00:46  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.13  2007/04/06 08:44:43  sherwin
* Update to make 2D regions work at StdRegions level
*
* Revision 1.12  2007/04/05 15:20:11  sherwin
* Updated 2D stuff to comply with SharedArray philosophy
*
* Revision 1.11  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.10  2007/01/17 16:36:58  pvos
* updating doxygen documentation
*
* Revision 1.9  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.8  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.7  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.6  2006/07/02 17:16:19  sherwin
*
* Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
* Revision 1.5  2006/06/13 18:05:02  sherwin
* Modifications to make MultiRegions demo ProjectLoc2D execute properly.
*
* Revision 1.4  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.3  2006/06/02 18:48:40  sherwin
* Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
*
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


