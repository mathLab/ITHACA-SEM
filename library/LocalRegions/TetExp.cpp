///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/TetExp.cpp,v $ 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/TetExp.h>

namespace Nektar
{
  namespace LocalRegions 
  {

    	TetExp::TetExp(const LibUtilities::BasisKey &Ba,
                   const LibUtilities::BasisKey &Bb,
		   const LibUtilities::BasisKey &Bc,
                   const SpatialDomains::TetGeomSharedPtr &geom):
            StdRegions::StdTetExp(Ba,Bb,Bc),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom = geom;

            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&TetExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&TetExp::CreateStaticCondMatrix, this, _1));
            }

            GenMetricInfo();
        }

	TetExp::TetExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb,
		       const LibUtilities::BasisKey &Bc
		      ):
            StdRegions::StdTetExp(Ba,Bb,Bc),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {

           for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&TetExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&TetExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors. 
            m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(); 
            int coordim = 3;
            Array<OneD,NekDouble> ndata = Array<OneD,NekDouble>(coordim*coordim*coordim,0.0); 
            ndata[0] = ndata[4] = 1.0; //TODO must check
            m_metricinfo->ResetGmat(ndata,1,3,coordim);//TODO must check
            m_metricinfo->ResetJac(1,ndata); //TODO must check
        }


	TetExp::TetExp(const TetExp &T):
            StdRegions::StdTetExp(T),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom       = T.m_geom;
            m_metricinfo = T.m_metricinfo;
        }


        TetExp::~TetExp()
        {
        }

	//TODO: implement
	void TetExp::GenMetricInfo()
        {
            SpatialDomains::GeomFactorsSharedPtr Xgfac;

            Xgfac = m_geom->GetGeomFactors();

            if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {
                int nq = GetTotPoints();   
                int coordim = m_geom->GetCoordim();
                int expdim = 2;
                SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;

                LibUtilities::BasisSharedPtr CBasis0;
                LibUtilities::BasisSharedPtr CBasis1;
//                 CBasis0 = m_geom->GetBasis(0,0); // assumes all goembasis are same
//                 CBasis1 = m_geom->GetBasis(0,1);
                int Cnq0 = CBasis0->GetNumPoints();
                int Cnq1 = CBasis1->GetNumPoints();

                ConstArray<OneD,NekDouble> ojac = Xgfac->GetJac();   
                ConstArray<TwoD,NekDouble> ogmat = Xgfac->GetGmat();
                Array<OneD,NekDouble> njac(nq);
                Array<OneD,NekDouble> ngmat(2*coordim*nq);

                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                    AllocateSharedPtr(gtype,expdim,coordim); 

                //basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
                   !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    int i;   
                    int nq0 = m_base[0]->GetNumPoints();
                    int nq1 = m_base[1]->GetNumPoints();

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

                        Interp1D(CBasis1->GetBasisKey(),&(normals[2])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[2])[i*nq1]);
                    }

                    m_metricinfo->ResetNormals(newnorm);

                    NEKERROR(ErrorUtil::ewarning,
                             "Need to check/debug routine for deformed elements");
                }
                else // Same data can be used 
                {                   
                    // Copy Jacobian
                    Blas::Dcopy(nq,&ojac[0],1,&njac[0],1);                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    ngmat = Array<OneD,NekDouble>(2*nq*coordim); 
                    Blas::Dcopy(2*coordim*nq,&ogmat[0][0],1,ngmat.data(),1);                    
                    m_metricinfo->ResetGmat(ngmat,nq,2,coordim);                 

                    m_metricinfo->ResetNormals(Xgfac->GetNormals());

                    NEKERROR(ErrorUtil::ewarning,
                             "Need to check/debug routine for deformed elements");
                }
            }
        } 

        //----------------------------
        // Integration Methods
        //----------------------------

	NekDouble TetExp::Integral(const ConstArray<OneD,NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
	    int    nquad2 = m_base[2]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            NekDouble retrunVal;
            Array<OneD,NekDouble> tmp   = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1, &tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble) jac[0], (NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            // call StdTetExp version;
            retrunVal = StdTetExp::Integral(tmp);

            return retrunVal; 
        }

	 void TetExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
                                      const ConstArray<OneD,NekDouble> &base1, 
				      const ConstArray<OneD,NekDouble> &base2, 
                                      const ConstArray<OneD,NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray)
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

            StdTetExp::IProductWRTBase(base0,base1,base2,tmp,outarray);
        }

	 void TetExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata(),
                            inarray,outarray);
        }


	void TetExp::FwdTrans(const ConstArray<OneD,NekDouble> & inarray,Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation())&&(m_base[2]->Collocation()))
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

	
	//TODO: implement
	void TetExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                               Array<OneD,NekDouble> &coords_1,
                               Array<OneD,NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
	    LibUtilities::BasisSharedPtr CBasis2;
            Array<OneD,NekDouble>  x;

            ASSERTL0(m_geom, "m_geom not define");

            // get physical points defined in Geom
           // m_geom->FillGeom(); //TODO: implement FillGeom()

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0, 
                         "output coords_2 is not defined");

		    //TODO: check GetBasis()
//                 CBasis0 = m_geom->GetBasis(2,0); 
//                 CBasis1 = m_geom->GetBasis(2,1);
//
//                 if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
//                    (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
//                 {
//                    x = m_geom->UpdatePhys(2); //TODO: check UpdatedPhys
//                     Blas::Dcopy(m_base[0]->GetNumPoints()*
//                                 m_base[1]->GetNumPoints(),
//                                 &x[0],1,&coords_2[0],1);
//                 }
//                 else // Interpolate to Expansion point distribution
//                 {
//                     Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
//                              m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_2[0]);
//                 }    
//             case 2:
//                 ASSERTL0(coords_1.num_elements(), 
//                          "output coords_1 is not defined");
//
//                 CBasis0 = m_geom->GetBasis(1,0); 
//                 CBasis1 = m_geom->GetBasis(1,1);
//
//                 if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
//                    (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
//                 {
//                     x = m_geom->UpdatePhys(1);
//                     Blas::Dcopy(m_base[0]->GetNumPoints()*
//                                 m_base[1]->GetNumPoints(),
//                                 &x[0],1,&coords_1[0],1);
//                 }
//                 else // Interpolate to Expansion point distribution
//                 {
//                     Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(1))[0],
//                              m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_1[0]);
//                 }
//             case 1:
//                 ASSERTL0(coords_0.num_elements(), 
//                          "output coords_0 is not defined");
//
//                 CBasis0 = m_geom->GetBasis(0,0); 
//                 CBasis1 = m_geom->GetBasis(0,1);
//
//                 if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
//                    (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
//                 {
//                     x = m_geom->UpdatePhys(0);
//                     Blas::Dcopy(m_base[0]->GetNumPoints()*
//                                 m_base[1]->GetNumPoints(),
//                                 &x[0],1,&coords_0[0],1);
//                 }
//                 else // Interpolate to Expansion point distribution
//                 {
//                     Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(0))[0],
//                              m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_0[0]);
//                 }
//                 break;
//             default:
//                 ASSERTL0(false,"Number of dimensions are greater than 2");
//                 break;
            }
        }  
  
	// get the coordinates "coords" at the local coordinates "Lcoords"
        void TetExp::GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, 
                              Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 && 
                     Lcoords[1] <= -1.0 && Lcoords[1] >= 1.0 &&
                     Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
                     "Local coordinates are not in region [-1,1]");

           // m_geom->FillGeom(); // TODO: implement FillGeom()

//             for(i = 0; i < m_geom->GetCoordDim(); ++i) //TODO: check GetCoordDim
//             {
//                 coords[i] = m_geom->GetCoord(i,Lcoords);
//             }
        }

	 void TetExp::WriteToFile(FILE *outfile)
        {
            int i,j,k;
            Array<OneD,NekDouble> coords[3];
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
	    int  nquad2 = m_base[2]->GetNumPoints();

            ASSERTL0(m_geom,"_geom not defined");  

	  //TODO: check GetCoordDim()
          // int  coordim   = m_geom->GetCoordDim();
	  int coordim;

            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            std::fprintf(outfile,"Variables = x");
            if(coordim == 2)
            {
              //  GetCoords(coords[0], coords[1]); //TODO: check GetCoords()
                fprintf(outfile,", y");
            }
            else if (coordim == 3)
            {
                GetCoords(coords[0],coords[1],coords[2]);
                fprintf(outfile,", y, z");
            }

            fprintf(outfile,", v\n");

            fprintf(outfile,"Zone, I=%d, J=%d, K=%d,F=Point\n",nquad0,nquad1, nquad2);

            for(i = 0; i < nquad0*nquad1*nquad2; ++i)
            {
                for(j = 0; j < coordim; ++j)
                {
		   for(k=0; k < coordim; ++k)
                   {
                     fprintf(outfile,"%lf ",coords[k][j]);
		   }
		   fprintf(outfile,"%lf \n",m_phys[j]);
                }
                fprintf(outfile,"%lf \n",m_phys[i]);
            }
        }

	void TetExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
        {
            int i,j,k;
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
	    int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD,NekDouble> coords[3];

            ASSERTL0(m_geom,"m_geom not defined");

            int     coordim  = m_geom->GetCoordim();

            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            GetCoords(coords[0],coords[1],coords[2]);

            if(dumpVar)
            { 
                outfile << "Variables = x";

                if(coordim == 2)
                {
                    outfile << ", y";
                }
                else if (coordim == 3)
                {
                    outfile << ", y, z";
                }
                outfile << ", v\n" << std::endl;
            }

            outfile << "Zone, I=" << nquad0 << ", J=" << 
                nquad1 << ", K=" << nquad2 << ", F=Point" << std::endl;

            for(i = 0; i < nquad0*nquad1*nquad2; ++i)
            {
                for(j = 0; j < coordim; ++j)
                {
		  for(k=0; k < coordim; ++k)
		  {
                     outfile << coords[k][j] << " ";
                  }
		    outfile << m_phys[j] << std::endl;
                }
                outfile << m_phys[i] << std::endl;
            }
        }

	DNekMatSharedPtr TetExp::GetStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
		LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();

                TetExpSharedPtr tmp = MemoryManager<TetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

                return tmp->StdTetExp::GetStdMatrix(mkey);
            }
            else
            {
                return StdTetExp::GetStdMatrix(mkey);
            }
        }

	 DNekBlkMatSharedPtr TetExp::GetStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
	        LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();

                TetExpSharedPtr tmp = MemoryManager<TetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

                return tmp->StdTetExp::GetStdStaticCondMatrix(mkey);                
            }
            else
            {
                return StdTetExp::GetStdStaticCondMatrix(mkey);
            }
        }

	NekDouble TetExp::PhysEvaluate(const ConstArray<OneD,NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(3);      

            ASSERTL0(m_geom,"m_geom not defined");
	
            //TODO: check GetLocCoords()
           // m_geom->GetLocCoords(coord, Lcoord, Lcoord);

            return StdTetExp::PhysEvaluate(Lcoord);
        }


	DNekScalMatSharedPtr TetExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());                      
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetShapeType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();
                        
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                        
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    { 
		      //TODO: make sure 3D Laplacian is set up for Hex in three-dimensional in Standard Region.
                     //   ASSERTL1(m_geom->GetCoordDim() == 2,"Standard Region Laplacian is only set up for Quads in two-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetShapeType(), *this);  
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetShapeType(), *this);  
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetShapeType(), *this);  
                        
                        DNekMatSharedPtr lap00 = GetStdMatrix(*lap00key.GetStdMatKey());
                        DNekMatSharedPtr lap01 = GetStdMatrix(*lap01key.GetStdMatKey());
                        DNekMatSharedPtr lap11 = GetStdMatrix(*lap11key.GetStdMatKey());
                        
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        ConstArray<TwoD,NekDouble> gmat = m_metricinfo->GetGmat();
                        
                        int rows = lap00->GetRows();
                        int cols = lap00->GetColumns();
                        
                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        
                        (*lap) = (gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0]) * (*lap00) + 
                            (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0]) * (*lap01 + Transpose(*lap01)) +
                            (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0]) * (*lap11);
                        
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetScaleFactor();
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetShapeType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetShapeType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }


	DNekScalBlkMatSharedPtr TetExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            unsigned int exp_size[] = {nbdry,nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHelmholtz: // special case since Helmholtz not defined in StdRegions

                // use Deformed case for both regular and deformed geometries 
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            default:
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                }
                else
                {
                    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(*(mkey.GetStdMatKey()));                    
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                break;
            UseLocRegionsMatrix:
                {
                    int i,j;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

                    const ConstArray<OneD,int> bmap = GetBoundaryMap();
                    const ConstArray<OneD,int> imap = GetInteriorMap();

                    for(i = 0; i < nbdry; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*A)(i,j) = mat(bmap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*B)(i,j) = mat(bmap[i],imap[j]);
                        }
                    }                    

                    for(i = 0; i < nint; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*C)(i,j) = mat(imap[i],bmap[j]);
                        }

                        for(j = 0; j < nint; ++j)
                        {
                            (*D)(i,j) = mat(imap[i],imap[j]);
                        }
                    }

                    // Calculate static condensed system 
                    if(nint)
                    {
                        D->Invert();
                        (*B) = (*B)*(*D);
                        (*A) = (*A) - (*B)*(*C);
                    }
                    
                    DNekScalMatSharedPtr     Atmp;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,A));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,C));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,D));

                }
            }
            
            return returnval;
        }



  
  }//end of namespace
}//end of namespace

/** 
 *    $Log: TetExp.cpp,v $
 *    Revision 1.2  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.9  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/



