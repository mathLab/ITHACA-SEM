///////////////////////////////////////////////////////////////////////////////
//
// File HexExp.cpp
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
// Description: Methods for Hex expansion in local regoins
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/HexExp.h>

namespace Nektar
{
  namespace LocalRegions 
  {

	HexExp::HexExp(const LibUtilities::BasisKey &Ba, 
                       const LibUtilities::BasisKey &Bb, 
		       const LibUtilities::BasisKey &Bc, 
                       const SpatialDomains::HexGeomSharedPtr &geom):
            StdRegions::StdHexExp(Ba,Bb,Bc),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom = geom;
            
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&HexExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&HexExp::CreateStaticCondMatrix, this, _1));
            }
            
            GenMetricInfo();
        }

	HexExp::HexExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb, 
		       const LibUtilities::BasisKey &Bc ):
            StdRegions::StdHexExp(Ba,Bb,Bc),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {

            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&HexExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&HexExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors. 
            m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(); 
            int coordim = 3;
            Array<OneD,NekDouble> ndata = Array<OneD,NekDouble>(coordim*coordim*coordim,0.0); 
            ndata[0] = ndata[4] = 1.0; //TODO must check
            m_metricinfo->ResetGmat(ndata,1,3,coordim); //TODO must check
            m_metricinfo->ResetJac(1,ndata); //TODO must check
        }

 	HexExp::HexExp(const HexExp &T):
            StdRegions::StdHexExp(T),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom       = T.m_geom;
            m_metricinfo = T.m_metricinfo;
        }

	// by default the StdHexExp destructor will be called      
	HexExp::~HexExp()
	{
	}

	//TODO: implement
	void HexExp::GenMetricInfo()
        {
            SpatialDomains::GeomFactorsSharedPtr Xgfac;
            
            Xgfac = m_geom->GetGeomFactors();
            
            if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {
                int nq0 = m_base[0]->GetNumPoints();
                int nq1 = m_base[1]->GetNumPoints();
                int nq = nq0*nq1;
                int coordim = m_geom->GetCoordim();
                int expdim = 2;
                SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;

                LibUtilities::BasisSharedPtr CBasis0;
                LibUtilities::BasisSharedPtr CBasis1;
     
                ConstArray<OneD,NekDouble> ojac = Xgfac->GetJac();
                ConstArray<TwoD,NekDouble> ogmat = Xgfac->GetGmat();
                Array<OneD,NekDouble> njac(nq);
                Array<OneD,NekDouble> ngmat(2*coordim*nq);
                
//                 CBasis0 = m_geom->GetBasis(0,0); // this assumes all goembasis are same
//                 CBasis1 = m_geom->GetBasis(0,1);
                int Cnq0 = CBasis0->GetNumPoints();
                int Cnq1 = CBasis1->GetNumPoints();

                
                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                    AllocateSharedPtr(gtype,expdim,coordim); 
                
                //basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
                   !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    int i;   

                    // interpolate Jacobian        
                    Interp2D(CBasis0->GetBasisKey(),
                             CBasis1->GetBasisKey(),
                             &ojac[0],
                             m_base[0]->GetBasisKey(),
                             m_base[1]->GetBasisKey(), 
                             &njac[0]);
                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    Array<OneD,NekDouble> dxdxi(nq);
                    for(i = 0; i < 2*coordim; ++i)
                    {
                        Vmath::Vmul(nq,&ojac[0],1,&ogmat[i][0],1,&dxdxi[0],1);
                        Interp2D(CBasis0->GetBasisKey(),
                                 CBasis1->GetBasisKey(), 
                                 &dxdxi[0], 
                                 m_base[0]->GetBasisKey(),
                                 m_base[1]->GetBasisKey(),
                                 &ngmat[0] + i*nq);
                        Vmath::Vdiv(nq,&ngmat[0]+i*nq,1,&njac[0],1,&ngmat[0]+i*nq,1);
                    }
                    
                    m_metricinfo->ResetGmat(ngmat,nq,2,coordim); 

                    // interpolate normals
                    Array<TwoD,NekDouble> newnorm = Array<TwoD,NekDouble>(4,coordim*max(nq0,nq1));    
                    ConstArray<TwoD,NekDouble> normals = Xgfac->GetNormals();
                    
                    for(i = 0; i < coordim; ++i)
                    {
                        Interp1D(CBasis0->GetBasisKey(),&(normals[0])[i*Cnq0],
                                 m_base[0]->GetBasisKey(),&(newnorm[0])[i*nq0]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[1])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[1])[i*nq1]);
                        
                        Interp1D(CBasis0->GetBasisKey(),&(normals[2])[i*Cnq0],
                                 m_base[0]->GetBasisKey(),&(newnorm[2])[i*nq0]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[3])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[3])[i*nq1]);
                    }
                    
                    m_metricinfo->ResetNormals(newnorm);

                    //                    NEKERROR(ErrorUtil::ewarning,
                    //"Need to check/debug routine for deformed elements");
                }
                else // Same data can be used 
                {                   
                    // Copy Jacobian
                    Blas::Dcopy(nq,&ojac[0],1,&njac[0],1);                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    Blas::Dcopy(2*coordim*nq,&ogmat[0][0],1,ngmat.data(),1);                    
                    m_metricinfo->ResetGmat(ngmat,nq,2,coordim);                 

                    m_metricinfo->ResetNormals(Xgfac->GetNormals());

                    NEKERROR(ErrorUtil::ewarning,
                             "Need to check/debug routine for deformed elements");
                }
                
            }
        }


 	NekDouble HexExp::Integral(const ConstArray<OneD,NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            NekDouble returnVal;
            Array<OneD,NekDouble> tmp   = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            
            // multiply inarray with Jacobian
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble) jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            
            // call StdHexExp version;
            returnVal = StdHexExp::Integral(tmp);
            
            return  returnVal;
        }
 
        void HexExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(), m_base[1]->GetBdata(), m_base[2]->GetBdata(), inarray,outarray);
        }

        void HexExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
                                     const ConstArray<OneD,NekDouble> &base1,
				     const ConstArray<OneD,NekDouble> &base2,  
                                     const ConstArray<OneD,NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray )
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
	    int    nquad2 = m_base[2]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            StdHexExp::IProductWRTBase(base0, base1, base2, tmp, outarray);
        }

	void HexExp::FwdTrans(const ConstArray<OneD,NekDouble> & inarray, Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);
                
                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetShapeType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }

	NekDouble HexExp::PhysEvaluate(const ConstArray<OneD,NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);
            
        ASSERTL0(m_geom,"m_geom not defined");
      //  m_geom->GetLocCoords(coord,Lcoord);
        
        return StdHexExp::PhysEvaluate(Lcoord);
        }

	// get the coordinates "coords" at the local coordinates "Lcoords"
        

	//TODO: implement FillGeom 
        void HexExp::GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, Array<OneD,NekDouble> &coords)
        {
            int  i;
            
            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 && 
                     Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0 &&
                     Lcoords[2] >= -1.0 && Lcoords[2]  <=1.0,
                     "Local coordinates are not in region [-1,1]");
            
//             m_geom->FillGeom();
//             
//             for(i = 0; i < m_geom->GetCoordDim(); ++i)
//             {
//                 coords[i] = m_geom->GetCoord(i,Lcoords);
//             }
        }

 	DNekMatSharedPtr HexExp::GetStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
	        LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
                HexExpSharedPtr tmp = MemoryManager<HexExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);
                
                return tmp->StdHexExp::GetStdMatrix(mkey);                
            }
            else
            {
                return StdHexExp::GetStdMatrix(mkey);
            }
        }

        DNekBlkMatSharedPtr HexExp::GetStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
	        LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
                HexExpSharedPtr tmp = MemoryManager<HexExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);
                
                return tmp->StdHexExp::GetStdStaticCondMatrix(mkey);                
            }
            else
            {
                return StdHexExp::GetStdStaticCondMatrix(mkey);
            }
        }

	//TODO: implement
	DNekScalBlkMatSharedPtr HexExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;
	    return returnval;
	}

	//TODO: implement
	DNekScalMatSharedPtr HexExp::CreateMatrix(const MatrixKey &mkey)
        {
         	DNekScalMatSharedPtr returnval;
		return returnval;
        }


  }//end of namespace
}//end of namespace

/** 
 *    $Log: HexExp.cpp,v $
 *    Revision 1.2  2007/07/20 00:45:50  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.1  2006/05/04 18:58:45  kirby
 *    *** empty log message ***
 *
 *    Revision 1.8  2006/03/12 07:43:31  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
