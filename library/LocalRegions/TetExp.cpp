///////////////////////////////////////////////////////////////////////////////
//
// File TetExp.cpp
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
#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/TetExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        TetExp::TetExp( const LibUtilities::BasisKey &Ba,
                        const LibUtilities::BasisKey &Bb,
                        const LibUtilities::BasisKey &Bc,
                        const SpatialDomains::TetGeomSharedPtr &geom
                        ):
            StdRegions::StdTetExp(Ba,Bb,Bc),
            m_geom(geom),
            m_metricinfo(),
            m_matrixManager(std::string("TetExpMatrix")),
            m_staticCondMatrixManager(std::string("TetExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoExpansionType,*this),
                                                boost::bind(&TetExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoExpansionType,*this),
                                                          boost::bind(&TetExp::CreateStaticCondMatrix, this, _1));
            }

            GenMetricInfo();
        }

        TetExp::TetExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb,
                       const LibUtilities::BasisKey &Bc
                       ):
            StdRegions::StdTetExp(Ba,Bb,Bc),
            m_geom(),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("TetExpMatrix")),
            m_staticCondMatrixManager(std::string("TetExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoExpansionType,*this),
                                                boost::bind(&TetExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoExpansionType,*this),
                                                          boost::bind(&TetExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors.
            int coordim = 3;
            Array<OneD,NekDouble> ndata = Array<OneD,NekDouble>(coordim*coordim*coordim,0.0);
            ndata[0] = ndata[4] = 1.0; //TODO must check
            m_metricinfo->ResetGmat(ndata, 1, 3, coordim);//TODO must check
            m_metricinfo->ResetJac(1, ndata); //TODO must check

        }


        TetExp::TetExp(const TetExp &T):
            StdRegions::StdTetExp(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(std::string("TetExpMatrix")),
            m_staticCondMatrixManager(std::string("TetExpStaticCondMatrix"))
        {
        } 

        TetExp::~TetExp()
        {
        }

        //TODO: check following computations and function
        void TetExp::GenMetricInfo()
        {
            m_metricinfo = m_geom->GetGeomFactors(m_base);
//             SpatialDomains::GeomFactorsSharedPtr Xgfac;

//             Xgfac = m_geom->GetGeomFactors();

//             if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
//             {
//                 m_metricinfo = Xgfac;
//             }
//             else
//             {
//                 int nq = GetTotPoints();   
//                 int coordim = m_geom->GetCoordim();
//                 int expdim = 3;
//                 SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;

//                 LibUtilities::BasisSharedPtr CBasis0;
//                 LibUtilities::BasisSharedPtr CBasis1;
//                 LibUtilities::BasisSharedPtr CBasis2;
//                 CBasis0 = m_geom->GetBasis(0,0); // assumes all goembasis are same
//                 CBasis1 = m_geom->GetBasis(0,1);
//                 CBasis2 = m_geom->GetBasis(0,2);
//                 int Cnq0 = CBasis0->GetNumPoints();
//                 int Cnq1 = CBasis1->GetNumPoints();
//                 int Cnq2 = CBasis2->GetNumPoints();

//                 Array<OneD, const NekDouble> ojac = Xgfac->GetJac();
//                 Array<TwoD, const NekDouble> ogmat = Xgfac->GetGmat();
//                 Array<OneD,NekDouble> njac(nq);
//                 Array<OneD,NekDouble> ndata(3*coordim*nq);//TODO: check ndata

//                 m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr(gtype,expdim,coordim);

//                 //basis are different distributions
//                 if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
//                    !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))||
//                    !(m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
//                 {
//                     int i;   
//                     int nq0 = m_base[0]->GetNumPoints();
//                     int nq1 = m_base[1]->GetNumPoints();
//                     int nq2 = m_base[2]->GetNumPoints();

//                     // interpolate Jacobian        
//                     LibUtilities::Interp3D(CBasis0->GetPointsKey(),
//                              CBasis1->GetPointsKey(),
// 		             CBasis2->GetPointsKey(),
//                              &ojac[0],
//                              m_base[0]->GetPointsKey(),
//                              m_base[1]->GetPointsKey(),
//                              m_base[2]->GetPointsKey(),
//                              &njac[0]);

//                     m_metricinfo->ResetJac(nq,njac);

//                     // interpolate Geometric data
//                     Array<OneD,NekDouble> dxdxi(nq);
//                     for(i = 0; i < 2*coordim; ++i) //TODO : find out why  2*coordim
//                     {
//                         Vmath::Vmul(nq,&ojac[0],1,&ogmat[i][0],1,&dxdxi[0],1);
//                         LibUtilities::Interp2D(CBasis0->GetPointsKey(),
//                                  CBasis1->GetPointsKey(), 
//                                  &dxdxi[0], 
//                                  m_base[0]->GetPointsKey(),
//                                  m_base[1]->GetPointsKey(),
//                                  &ndata[0] + i*nq);
//                         Vmath::Vdiv(nq,&ndata[0]+i*nq,1,&njac[0],1,&ndata[0]+i*nq,1);
//                     }
//                     m_metricinfo->ResetGmat(ndata,nq,3,coordim); 

//                     NEKERROR(ErrorUtil::ewarning,
//                              "Need to check/debug routine for deformed elements");
//                 }
//                 else // Same data can be used 
//                 {                   
//                     // Copy Jacobian
//                     Blas::Dcopy(nq,&ojac[0],1,&njac[0],1);
//                     m_metricinfo->ResetJac(nq,njac);

//                     // interpolate Geometric data
//                     ndata = Array<OneD,NekDouble>(3*nq*coordim); //TODO: check this function
//                     Blas::Dcopy(3*coordim*nq, &ogmat[0][0], 1, ndata.data(), 1); //TODO: check this function

//                     m_metricinfo->ResetGmat(ndata,nq,3,coordim);

//                     NEKERROR(ErrorUtil::ewarning,
//                              "Need to check/debug routine for deformed elements");
//                 }
//             }
        }

        //----------------------------
        // Integration Methods
        //----------------------------

        /** \brief Integrate the physical point list \a inarray over region
            and return the value

            Inputs:\n

            - \a inarray: definition of function to be returned at quadrature point
            of expansion.

            Outputs:\n

            - returns \f$\int^1_{-1}\int^1_{-1} \int^1_{-1} u(\eta_1, \eta_2, \eta_3) J[i,j,k] d
            \eta_1 d \eta_2 d \eta_3 \f$ where \f$inarray[i,j,k] = u(\eta_{1i},\eta_{2j},\eta_{3k})
            \f$ and \f$ J[i,j,k] \f$ is the Jacobian evaluated at the quadrature point.
        */
        NekDouble TetExp::Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
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

        void TetExp::IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                     const Array<OneD, const NekDouble>& base1, 
                                     const Array<OneD, const NekDouble>& base2, 
                                     const Array<OneD, const NekDouble>& inarray, 
                                     Array<OneD, NekDouble> & outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
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

        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////
        
        /** 
            \brief Calculate the deritive of the physical points 
        **/
        void TetExp::PhysDeriv(const Array<OneD, const NekDouble> & inarray,
                               Array<OneD,NekDouble> &out_d0,
                               Array<OneD,NekDouble> &out_d1,
                               Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nquad2 = m_base[2]->GetNumPoints();
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);
            Array<OneD,NekDouble> Diff2 = Array<OneD,NekDouble>(nquad0*nquad1*nquad2);

            StdTetExp::PhysDeriv(inarray, Diff0, Diff1, Diff2);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[0][0],1,&Diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[1][0],1,&Diff1[0],1, &out_d0[0], 1,&out_d0[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[2][0],1,&Diff2[0],1, &out_d0[0], 1,&out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[3][0],1,&Diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[4][0],1,&Diff1[0],1, &out_d1[0], 1,&out_d1[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[5][0],1,&Diff2[0],1, &out_d1[0], 1,&out_d1[0],1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[6][0],1,&Diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[7][0],1,&Diff1[0],1, &out_d2[0], 1, &out_d2[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[8][0],1,&Diff2[0],1, &out_d2[0], 1, &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[1][0],&Diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[2][0],&Diff2[0],1, &out_d0[0], 1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[3][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[4][0],&Diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[5][0],&Diff2[0],1, &out_d1[0], 1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1*nquad2,gmat[6][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[7][0],&Diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[8][0],&Diff2[0],1, &out_d2[0], 1);
                }
            }
        }


        /** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs  
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - \a (this)->m_coeffs: updated array of expansion coefficients.             
        */
        
        void TetExp::FwdTrans(const Array<OneD, const NekDouble> & inarray,Array<OneD,NekDouble> &outarray)
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
                                              DetExpansionType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];

                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);

                out = (*matsys)*in;
            }
        }


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
            //             m_geom->FillGeom();  //TODO: implement
            
            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements(), "output coords_2 is not defined");
                CBasis0 = m_geom->GetBasis(2,0); 
                CBasis1 = m_geom->GetBasis(2,1);
                CBasis2 = m_geom->GetBasis(2,2);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_2, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(2))[0],
                             m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), &coords_2[0]);
                }    
            case 2:
                ASSERTL0(coords_1.num_elements(), "output coords_1 is not defined");
                
                CBasis0 = m_geom->GetBasis(1,0); 
                CBasis1 = m_geom->GetBasis(1,1);
                CBasis2 = m_geom->GetBasis(1,2);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_1, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(1))[0],
                             m_base[0]->GetPointsKey(), m_base[1]->GetPointsKey(), m_base[2]->GetPointsKey(), &coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements(), "output coords_0 is not defined");
                
                CBasis0 = m_geom->GetBasis(0,0); 
                CBasis1 = m_geom->GetBasis(0,1);
                CBasis2 = m_geom->GetBasis(0,2);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey()))&&
                   (m_base[2]->GetBasisKey().SamePoints(CBasis2->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),
                                x, 1, coords_0, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp3D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), CBasis2->GetPointsKey(), &(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),m_base[2]->GetPointsKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 3");
                break;
            }
        }

	
      
        // get the coordinates "coords" at the local coordinates "Lcoords"
        void TetExp::GetCoord(const Array<OneD, const NekDouble> &Lcoords, Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 && 
                     Lcoords[1] <= -1.0 && Lcoords[1] >= 1.0 &&
                     Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
                     "Local coordinates are not in region [-1,1]");

            // m_geom->FillGeom(); // TODO: implement FillGeom()

            for(i = 0; i < m_geom->GetCoordDim(); ++i) 
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        void TetExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
        {
            if(format==eTecplot)
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
                    outfile << ", "<< var << std::endl << std::endl;
                }
                
                outfile << "Zone, I=" << nquad0 << ", J=" << nquad1 << ", K=" << nquad2 << ", F=Point" << std::endl;
                
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
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }
      
        NekDouble TetExp::PhysEvaluate(const Array<OneD, const NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(3);
            
            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);
            
            return StdTetExp::PhysEvaluate(Lcoord);
        }


        DNekMatSharedPtr TetExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
            StdRegions::StdTetExpSharedPtr tmp = MemoryManager<StdTetExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);

            return tmp->GetStdMatrix(mkey);
        }


        DNekMatSharedPtr TetExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
            case StdRegions::eHybridDGHelmBndLam:
                returnval = Expansion3D::GenMatrix(mkey);
                break;
            default:
                returnval = StdTetExp::GenMatrix(mkey);
            }
            
            return returnval;            
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
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),
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
                        // TODO: make sure 3D Laplacian is set up for StdHexExp.
                        // ASSERTL1(m_geom->GetCoordDim() == 2,"Standard Region Laplacian is only set up for Quads in two-dimensional");
                        ASSERTL1(m_geom->GetCoordDim() == 3,"Standard Region Laplacian is only set up for Hex in three-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);

                        DNekMatSharedPtr lap00 = GetStdMatrix(*lap00key.GetStdMatKey());
                        DNekMatSharedPtr lap01 = GetStdMatrix(*lap01key.GetStdMatKey());
                        DNekMatSharedPtr lap11 = GetStdMatrix(*lap11key.GetStdMatKey());

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

                        int rows = lap00->GetRows();
                        int cols = lap00->GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0]) * (*lap00) +
                            (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0]) * (*lap01 + Transpose(*lap01)) +
                            (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0]) * (*lap11);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac, lap);
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstant(0);
                    MatrixKey masskey(StdRegions::eMass, mkey.GetExpansionType(), *this);
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian, mkey.GetExpansionType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows, cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;

                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one, helm);
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

            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            unsigned int exp_size[] = {nbdry, nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks, nblks, exp_size, exp_size); //Really need a constructor which takes Arrays
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
                    
                    //TODO: check below
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

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    GetBoundaryMap(bmap);
                    GetInteriorMap(imap); 

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
 *    Revision 1.19  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.18  2009/01/21 16:59:57  pvos
 *    Added additional geometric factors to improve efficiency
 *
 *    Revision 1.17  2008/09/09 15:05:09  sherwin
 *    Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
 *
 *    Revision 1.16  2008/08/14 22:12:57  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.15  2008/07/09 11:44:49  sherwin
 *    Replaced GetScaleFactor call with GetConstant(0)
 *
 *    Revision 1.14  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.13  2008/06/14 01:20:53  ehan
 *    Clean up the codes
 *
 *    Revision 1.12  2008/06/06 23:25:21  ehan
 *    Added doxygen documentation
 *
 *    Revision 1.11  2008/06/05 20:18:47  ehan
 *    Fixed undefined function GetGtype() in the ASSERTL2().
 *
 *    Revision 1.10  2008/06/02 23:35:26  ehan
 *    Fixed warning : no new line at end of file
 *
 *    Revision 1.9  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.8  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.7  2008/05/29 01:02:13  bnelson
 *    Added precompiled header support.
 *
 *    Revision 1.6  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.5  2008/03/17 10:36:17  pvos
 *    Clean up of the code
 *
 *    Revision 1.4  2008/02/16 05:52:49  ehan
 *    Added PhysDeriv and virtual functions.
 *
 *    Revision 1.3  2008/02/05 00:40:57  ehan
 *    Added initial tetrahedral expansion.
 *
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
