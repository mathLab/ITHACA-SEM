///////////////////////////////////////////////////////////////////////////////
//
// File PrismExp.cpp
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
// Description:  PrismExp routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/PrismExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        PrismExp::PrismExp(const LibUtilities::BasisKey &Ba,
                           const LibUtilities::BasisKey &Bb,
                           const LibUtilities::BasisKey &Bc,
                           const SpatialDomains::PrismGeomSharedPtr &geom):
            StdExpansion  (StdRegions::StdPrismData::getNumberOfCoefficients(
                               Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()),
                           3, Ba, Bb, Bc),
            Expansion     (),
            StdExpansion3D(StdRegions::StdPrismData::getNumberOfCoefficients(
                               Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes()), 
                           Ba, Bb, Bc),
            Expansion3D   (),
            StdPrismExp   (Ba, Bb, Bc),
            m_geom        (geom),
            m_metricinfo  (m_geom->GetGeomFactors(m_base)),
            m_matrixManager(
                    boost::bind(&PrismExp::CreateMatrix, this, _1),
                    std::string("PrismExpMatrix")),
            m_staticCondMatrixManager(
                    boost::bind(&PrismExp::CreateStaticCondMatrix, this, _1),
                    std::string("PrismExpStaticCondMatrix"))
        {
        }

        PrismExp::PrismExp(const PrismExp &T):
            StdExpansion(T),
            Expansion(T),
            StdExpansion3D(T),
            Expansion3D(T),
            StdRegions::StdPrismExp(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(T.m_matrixManager),
            m_staticCondMatrixManager(T.m_staticCondMatrixManager)
        {
        } 

        PrismExp::~PrismExp()
        {
        }


        //-------------------------------
        // Integration Methods
        //-------------------------------
        
        /** 
         * \brief Integrate the physical point list \a inarray over prismatic
         * region and return the value.
         *
         * Inputs:\n
         *
         * - \a inarray: definition of function to be returned at quadrature
         * point of expansion.
         *
         * Outputs:\n
         * 
         * - returns \f$\int^1_{-1}\int^1_{-1}\int^1_{-1} u(\bar \eta_1,
         *  \xi_2, \xi_3) J[i,j,k] d \bar \eta_1 d \xi_2 d \xi_3 \f$ \n \f$ =
         *  \sum_{i=0}^{Q_1 - 1} \sum_{j=0}^{Q_2 - 1} \sum_{k=0}^{Q_3 - 1}
         *  u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0})w_{i}^{0,0}
         *  w_{j}^{0,0} \hat w_{k}^{1,0} \f$ \n where \f$ inarray[i,j, k] =
         *  u(\bar \eta_{1i}^{0,0}, \xi_{2j}^{0,0},\xi_{3k}^{1,0}) \f$, \n
         *  \f$\hat w_{i}^{1,0} = \frac {w_{j}^{1,0}} {2} \f$ \n and \f$
         *  J[i,j,k] \f$ is the Jacobian evaluated at the quadrature point.
        */
        NekDouble PrismExp::v_Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,       NekDouble> tmp(nquad0*nquad1*nquad2);

            // Multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,(NekDouble)jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            
            // Call StdPrismExp version.
            return StdPrismExp::v_Integral(tmp);
        }


        //----------------------------
        // Differentiation Methods
        //----------------------------
        void PrismExp::v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD,       NekDouble>& out_d0,
                                         Array<OneD,       NekDouble>& out_d1,
                                         Array<OneD,       NekDouble>& out_d2)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,       NekDouble> diff0(nquad0*nquad1*nquad2);
            Array<OneD,       NekDouble> diff1(nquad0*nquad1*nquad2);
            Array<OneD,       NekDouble> diff2(nquad0*nquad1*nquad2);

            StdPrismExp::v_PhysDeriv(inarray, diff0, diff1, diff2);
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[0][0],1,&diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[1][0],1,&diff1[0],1, &out_d0[0], 1,&out_d0[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[2][0],1,&diff2[0],1, &out_d0[0], 1,&out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[3][0],1,&diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[4][0],1,&diff1[0],1, &out_d1[0], 1,&out_d1[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[5][0],1,&diff2[0],1, &out_d1[0], 1,&out_d1[0],1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1*nquad2,&gmat[6][0],1,&diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[7][0],1,&diff1[0],1, &out_d2[0], 1, &out_d2[0],1);
                    Vmath::Vvtvp (nquad0*nquad1*nquad2,&gmat[8][0],1,&diff2[0],1, &out_d2[0], 1, &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul (nquad0*nquad1*nquad2,gmat[0][0],&diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[1][0],&diff1[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[2][0],&diff2[0],1, &out_d0[0], 1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Smul (nquad0*nquad1*nquad2,gmat[3][0],&diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[4][0],&diff1[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[5][0],&diff2[0],1, &out_d1[0], 1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Smul (nquad0*nquad1*nquad2,gmat[6][0],&diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[7][0],&diff1[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1*nquad2,gmat[8][0],&diff2[0],1, &out_d2[0], 1);
                }
            }
        }

        //---------------------------------------
        // Transforms
        //---------------------------------------
        
        /** 
         * \brief Forward transform from physical quadrature space stored in
         * \a inarray and evaluate the expansion coefficients and store in \a
         * (this)->m_coeffs
         * 
         * Inputs:\n
         *
         * - \a inarray: array of physical quadrature points to be transformed
         *
         * Outputs:\n
         *
         * - (this)->_coeffs: updated array of expansion coefficients.
         */
        void PrismExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                        Array<OneD,       NekDouble>& outarray)
        {
            if(m_base[0]->Collocation() && 
               m_base[1]->Collocation() && 
               m_base[2]->Collocation())
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
            }
            else
            {
                v_IProductWRTBase(inarray, outarray);

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


        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /**
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0*base1*base2 and put into outarray:
         *            
         * \f$ \begin{array}{rcl} I_{pqr} = (\phi_{pqr}, u)_{\delta} & = &
         * \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \sum_{k=0}^{nq_2} \psi_{p}^{a}
         * (\bar \eta_{1i}) \psi_{q}^{a} (\xi_{2j}) \psi_{pr}^{b} (\xi_{3k})
         * w_i w_j w_k u(\bar \eta_{1,i} \xi_{2,j} \xi_{3,k}) J_{i,j,k}\\ & =
         * & \sum_{i=0}^{nq_0} \psi_p^a(\bar \eta_{1,i}) \sum_{j=0}^{nq_1}
         * \psi_{q}^a(\xi_{2,j}) \sum_{k=0}^{nq_2} \psi_{pr}^b u(\bar
         * \eta_{1i},\xi_{2j},\xi_{3k}) J_{i,j,k} \end{array} \f$ \n
         *  
         * where
         *  
         * \f$ \phi_{pqr} (\xi_1 , \xi_2 , \xi_3) = \psi_p^a (\bar \eta_1)
         * \psi_{q}^a (\xi_2) \psi_{pr}^b (\xi_3) \f$ \n
         * 
         * which can be implemented as \n \f$f_{pr} (\xi_{3k}) =
         * \sum_{k=0}^{nq_3} \psi_{pr}^b u(\bar \eta_{1i},\xi_{2j},\xi_{3k})
         * J_{i,j,k} = {\bf B_3 U} \f$ \n \f$ g_{q} (\xi_{3k}) =
         * \sum_{j=0}^{nq_1} \psi_{q}^a (\xi_{2j}) f_{pr} (\xi_{3k}) = {\bf
         * B_2 F} \f$ \n \f$ (\phi_{pqr}, u)_{\delta} = \sum_{k=0}^{nq_0}
         * \psi_{p}^a (\xi_{3k}) g_{q} (\xi_{3k}) = {\bf B_1 G} \f$
         */
        void PrismExp::v_IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                         const Array<OneD, const NekDouble>& base1, 
                                         const Array<OneD, const NekDouble>& base2, 
                                         const Array<OneD, const NekDouble>& inarray, 
                                               Array<OneD,       NekDouble>& outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            int nquad2 = m_base[2]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,       NekDouble> tmp(nquad0*nquad1*nquad2);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1*nquad2,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1*nquad2,jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }

            StdPrismExp::v_IProductWRTBase(base0,base1,base2,tmp,outarray,1);
        }

        void PrismExp::v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
                                               Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase(m_base[0]->GetBdata(),
                              m_base[1]->GetBdata(),
                              m_base[2]->GetBdata(),
                              inarray, outarray);
        }
        
        
        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        void PrismExp::v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                   Array<OneD, NekDouble> &coords_1,
                                   Array<OneD, NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            LibUtilities::BasisSharedPtr CBasis2;
            Array<OneD,NekDouble>  x;
            
            ASSERTL0(m_geom, "m_geom not defined");
            
            // get physical points defined in Geom
            m_geom->FillGeom();  

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
                        //Blas::Dcopy(GetTotPoints(), x, 1, coords_2, 1);
                        Blas::Dcopy(m_base[0]->GetNumPoints()*m_base[1]->GetNumPoints()*m_base[2]->GetNumPoints(),x, 1, coords_2, 1);
                    }
                    else // LibUtilities::Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(),
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(2))[0],
                                               m_base[0]->GetPointsKey(), 
                                               m_base[1]->GetPointsKey(), 
                                               m_base[2]->GetPointsKey(), 
                                               &coords_2[0]);
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
                        Blas::Dcopy(GetTotPoints(), x, 1, coords_1, 1);
                    }
                    else // Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(), 
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(1))[0],
                                               m_base[0]->GetPointsKey(), 
                                               m_base[1]->GetPointsKey(), 
                                               m_base[2]->GetPointsKey(), 
                                               &coords_1[0]);
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
                        Blas::Dcopy(GetTotPoints(), x, 1, coords_0, 1);
                    }
                    else // Interpolate to Expansion point distribution
                    {
                        LibUtilities::Interp3D(CBasis0->GetPointsKey(), 
                                               CBasis1->GetPointsKey(), 
                                               CBasis2->GetPointsKey(), 
                                               &(m_geom->UpdatePhys(0))[0],
                                               m_base[0]->GetPointsKey(),
                                               m_base[1]->GetPointsKey(),
                                               m_base[2]->GetPointsKey(),
                                               &coords_0[0]);
                    }
                    break;
                default:
                    ASSERTL0(false,"Number of dimensions are greater than 3");
                    break;
            }
        }

        /**
         * @brief Get the coordinates #coords at the local coordinates
         * #Lcoords.
         */
        void PrismExp::v_GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                                        Array<OneD,       NekDouble>& coords)
        {
            int i;
            
            ASSERTL1(Lcoords[0] <= -1.0 && Lcoords[0] >= 1.0 && 
                     Lcoords[1] <= -1.0 && Lcoords[1] >= 1.0 &&
                     Lcoords[2] <= -1.0 && Lcoords[2] >= 1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();
            
            for(i = 0; i < m_geom->GetCoordDim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }

        NekDouble PrismExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord)
        {
            return v_PhysEvaluate(coord, m_phys);
        }

        NekDouble PrismExp::v_PhysEvaluate(const Array<OneD, const NekDouble>& coord, 
                                           const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> Lcoord(3);

            ASSERTL0(m_geom,"m_geom not defined");
	
            m_geom->GetLocCoords(coord, Lcoord);

            return StdPrismExp::PhysEvaluate(Lcoord, physvals);
        }


        //---------------------------------------
        // Helper functions
        //---------------------------------------

        const SpatialDomains::GeomFactorsSharedPtr& PrismExp::v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

        const SpatialDomains::GeometrySharedPtr PrismExp::v_GetGeom() const
        {
            return m_geom;
        }

        const SpatialDomains::Geometry3DSharedPtr& PrismExp::v_GetGeom3D() const
        {
            return m_geom;
        }
        
        int PrismExp::v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        void PrismExp::v_WriteToFile(std::ofstream &outfile, 
                                     OutputFormat   format, 
                                     const bool     dumpVar, 
                                     std::string    var)
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
	
        
        //---------------------------------------
        // Matrix creation functions
        //---------------------------------------
        
        DNekMatSharedPtr PrismExp::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                    returnval = StdPrismExp::GenMatrix(mkey);
            }
            
            return returnval;            
        }
        
        DNekMatSharedPtr PrismExp::v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            LibUtilities::BasisKey bkey2 = m_base[2]->GetBasisKey();
            StdRegions::StdPrismExpSharedPtr tmp = 
                MemoryManager<StdPrismExp>::AllocateSharedPtr(bkey0, bkey1, bkey2);
            
            return tmp->GetStdMatrix(mkey); 
        }

        DNekScalMatSharedPtr& PrismExp::v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        DNekScalBlkMatSharedPtr& PrismExp::v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }

        DNekScalMatSharedPtr PrismExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,
                     "Geometric information is not set up");
            
            switch(mkey.GetMatrixType())
            {
                case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {   
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
                    }
                    break;
                }
                
                case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),*this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble fac = 1.0/(m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(mkey);
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);                        
                    }
                    break;
                }
                
                case StdRegions::eWeakDeriv0:
                case StdRegions::eWeakDeriv1:
                case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);
                        
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
                        int dir;

                        switch(mkey.GetMatrixType())
                        {
                            case StdRegions::eWeakDeriv0:
                                dir = 0;
                                break;
                            case StdRegions::eWeakDeriv1:
                                dir = 1;
                                break;
                            case StdRegions::eWeakDeriv2:
                                dir = 2;
                                break;
                        }

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);  
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetExpansionType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(deriv0key);
                        DNekMat &deriv1 = *GetStdMatrix(deriv1key);
                        
                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = gmat[2*dir][0]*deriv0 + gmat[2*dir+1][0]*deriv1;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,WeakDeriv);
                    }
                    break;
                }
                
                case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(mkey);

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        //ASSERTL1(m_geom->GetCoordDim() == 2,"Standard Region Laplacian is only set up for Quads in two-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap02key(StdRegions::eLaplacian02,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap12key(StdRegions::eLaplacian12,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap22key(StdRegions::eLaplacian22,
                                           mkey.GetExpansionType(), *this);

                        DNekMat &lap00
                                    = *GetStdMatrix(lap00key);
                        DNekMat &lap01
                                    = *GetStdMatrix(lap01key);
                        DNekMat &lap02
                                    = *GetStdMatrix(lap02key);
                        DNekMat &lap11
                                    = *GetStdMatrix(lap11key);
                        DNekMat &lap12
                                    = *GetStdMatrix(lap12key);
                        DNekMat &lap22
                                    = *GetStdMatrix(lap22key);

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat
                                                    = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>
                                                ::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0] + gmat[3][0]*gmat[3][0]
                                        + gmat[6][0]*gmat[6][0])*lap00
                               + (gmat[1][0]*gmat[1][0] + gmat[4][0]*gmat[4][0]
                                        + gmat[7][0]*gmat[7][0])*lap11
                               + (gmat[2][0]*gmat[2][0] + gmat[5][0]*gmat[5][0]
                                        + gmat[8][0]*gmat[8][0])*lap22
                               + (gmat[0][0]*gmat[1][0] + gmat[3][0]*gmat[4][0]
                                        + gmat[6][0]*gmat[7][0])
                                 *(lap01 + Transpose(lap01))
                               + (gmat[0][0]*gmat[2][0] + gmat[3][0]*gmat[5][0]
                                        + gmat[6][0]*gmat[8][0])
                                 *(lap02 + Transpose(lap02))
                               + (gmat[1][0]*gmat[2][0] + gmat[4][0]*gmat[5][0]
                                        + gmat[7][0]*gmat[8][0])
                                 *(lap12 + Transpose(lap12));

                        returnval = MemoryManager<DNekScalMat>
                                                ::AllocateSharedPtr(jac,lap);
                        /*
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);  
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);  

                        DNekMat &lap00 = *GetStdMatrix(lap00key);
                        DNekMat &lap01 = *GetStdMatrix(lap01key);
                        DNekMat &lap11 = *GetStdMatrix(lap11key);

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0]+gmat[2][0]*gmat[2][0])*lap00 + 
                            (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0])*(lap01 + Transpose(lap01)) +
                            (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0])*lap11;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,lap);
                        */
                    }
                    break;
                }
                
                case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                    break;
                }

                case StdRegions::eHybridDGHelmholtz:
                case StdRegions::eHybridDGLamToU:
                case StdRegions::eHybridDGLamToQ0:
                case StdRegions::eHybridDGLamToQ1:
                case StdRegions::eHybridDGHelmBndLam:
                {
                    NekDouble one    = 1.0;
                    
                    DNekMatSharedPtr mat = GenMatrix(mkey);
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    
                    break;
                }
                
                case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;
                    
                    MatrixKey hkey(StdRegions::eHybridDGHelmholtz, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
//                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
//                                                  DetExpansionType(),*this,
//                                                  mkey.GetConstant(0),
//                                                  mkey.GetConstant(1));
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    break;
                }
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }
                
            return returnval;
        }

        DNekScalBlkMatSharedPtr PrismExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() != SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            unsigned int exp_size[] = {nbdry,nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eLaplacian: 
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
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(mkey);
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
