///////////////////////////////////////////////////////////////////////////////
//
// File Expansion3D.cpp
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
// Description: File for Expansion3D routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <SpatialDomains/Geometry3D.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        Expansion3D::Expansion3D(){}
        
        void Expansion3D::v_AddHDGHelmholtzTraceTerms(
            const NekDouble                                tau,
            const Array<OneD, const NekDouble>            &inarray, 
            Array<OneD,StdRegions::StdExpansionSharedPtr> &FaceExp,  
            const StdRegions::VarCoeffMap                 &dirForcing,
            Array<OneD,NekDouble>                         &outarray)
        {
            
            ASSERTL0(&inarray[0] != &outarray[0],
                     "Input and output arrays use the same memory");
            
            int f,cnt;
            int order_f;
            int nfaces = GetNfaces();
            Array<OneD, const NekDouble> tmp;
            
            cnt = 0;
            for(f = 0; f < nfaces; ++f)
            {
                order_f = FaceExp[f]->GetNcoeffs();                    
                Vmath::Vcopy(order_f,tmp = inarray + cnt, 1, FaceExp[f]->UpdateCoeffs(), 1);
                FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(), FaceExp[f]->UpdatePhys());
                AddHDGHelmholtzFaceTerms(tau, f, FaceExp[f], dirForcing, outarray);
                cnt += order_f;
            }
        }
        
        //  evaluate additional terms in HDG face. Note that this assumes that
        // edges are unpacked into local cartesian order. 
        void Expansion3D::v_AddHDGHelmholtzFaceTerms(
            const NekDouble                    tau,
            const int                          face,
            StdRegions::StdExpansionSharedPtr  FaceExp,
            const StdRegions::VarCoeffMap     &varcoeffs,
            Array<OneD,NekDouble>             &outarray)
        {
            int i,j,n;
            int nquad_f = FaceExp->GetNumPoints(0)*FaceExp->GetNumPoints(1); 
            int order_f = FaceExp->GetNcoeffs();
            int coordim = GetCoordim();
            int ncoeffs = GetNcoeffs();

            Array<OneD, NekDouble> inval   (nquad_f);
            Array<OneD, NekDouble> outcoeff(order_f);
            Array<OneD, NekDouble> tmpcoeff(ncoeffs);

            const Array<OneD, const Array<OneD, NekDouble> > &normals
                = GetFaceNormal(face);

            DNekScalMat &invMass = *GetLocMatrix(StdRegions::eInvMass);
            
            DNekVec Coeffs(ncoeffs,outarray,eWrapper);
            DNekVec Tmpcoeff(ncoeffs,tmpcoeff,eWrapper);

            StdRegions::IndexMapKey ikey(
                StdRegions::eFaceToElement, DetExpansionType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetFaceOrient(face));
            StdRegions::IndexMapValuesSharedPtr map = 
                StdExpansion::GetIndexMap(ikey);

            StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                   StdRegions::eWeakDeriv1,
                                                   StdRegions::eWeakDeriv2};

            // @TODO Variable coefficients
            /*
            StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
                                                    StdRegions::eVarCoeffD11,
                                                    StdRegions::eVarCoeffD22};
            Array<OneD, NekDouble> varcoeff_work(nquad_f);
            StdRegions::VarCoeffMap::const_iterator x;
            ///// @TODO: What direction to use here??
            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
            {
                GetPhysFaceVarCoeffsFromElement(face,FaceExp,x->second,varcoeff_work);
                Vmath::Vmul(nquad_f,varcoeff_work,1,FaceExp->GetPhys(),1,FaceExp->UpdatePhys(),1);
            }
            */

            //================================================================
            // Add F = \tau <phi_i,in_phys>
            // Fill face and take inner product
            FaceExp->IProductWRTBase(FaceExp->GetPhys(),
                                     FaceExp->UpdateCoeffs());

            // add data to out array
            const Array<OneD, const NekDouble> &faceCoeffs = 
                FaceExp->GetCoeffs();

            for(i = 0; i < order_f; ++i)
            {
                outarray[(*map)[i].index] += (*map)[i].sign*tau*faceCoeffs[i];
            }
            //================================================================

            NekDouble scale = invMass.Scale();
            const NekDouble *data = invMass.GetRawPtr();

            //===============================================================
            // Add -\sum_i D_i^T M^{-1} G_i + E_i M^{-1} G_i = 
            //                         \sum_i D_i M^{-1} G_i term

            // Three independent direction
            for(n = 0; n < coordim; ++n)
            {
                Vmath::Vmul(nquad_f,normals[n],1,FaceExp->GetPhys(),1,inval,1);
                
                if (m_negatedNormals[face])
                {
                    Vmath::Neg(nquad_f, inval, 1);
                }

                // @TODO Multiply by variable coefficients
                // @TODO: Document this (probably not needed)
                /*
                StdRegions::VarCoeffMap::const_iterator x;
                if ((x = varcoeffs.find(VarCoeff[n])) != varcoeffs.end())
                {
                    GetPhysEdgeVarCoeffsFromElement(edge,FaceExp,x->second,varcoeff_work);
                    Vmath::Vmul(nquad_f,varcoeff_work,1,FaceExp->GetPhys(),1,FaceExp->UpdatePhys(),1);
                }
                */
                
                FaceExp->IProductWRTBase(inval,outcoeff);
                
                // M^{-1} G
                for(i = 0; i < ncoeffs; ++i)
                {
                    tmpcoeff[i] = 0;
                    for(j = 0; j < order_f; ++j)
                    {
                        tmpcoeff[i] += scale*data[i+(*map)[j].index*ncoeffs]*(*map)[j].sign*outcoeff[j];
                    }
                }
                
                DNekScalMat &Dmat = *GetLocMatrix(DerivType[n]);
                Coeffs = Coeffs  + Dmat*Tmpcoeff;       
                
                /*
                if(varcoeffs.find(VarCoeff[n]) != varcoeffs.end())
                {
                    MatrixKey mkey(DerivType[n], DetExpansionType(), *this, StdRegions::NullConstFactorMap, varcoeffs);
                    DNekScalMat &Dmat = *GetLocMatrix(mkey);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;                 
                }

                else
                {
                    DNekScalMat &Dmat = *GetLocMatrix(DerivType[n]);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;       
                }
                */
            }
        }

        /**
         * Computes the C matrix entries due to the presence of the identity
         * matrix in Eqn. 32.
         */
        void Expansion3D::v_AddNormTraceInt(const int dir,
                                          Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,StdRegions::StdExpansionSharedPtr> &FaceExp,
                                          Array<OneD,NekDouble> &outarray,
                                          const StdRegions::VarCoeffMap &varcoeffs)
        {
            int i,f,cnt;
            int order_f,nquad_f;
            int nfaces = GetNfaces();
            int coordim = GetCoordim();

            cnt = 0;
            for(f = 0; f < nfaces; ++f)
            {
                order_f = FaceExp[f]->GetNcoeffs();
                nquad_f = FaceExp[f]->GetNumPoints(0)*FaceExp[f]->GetNumPoints(1);

                const Array<OneD, const Array<OneD, NekDouble> > normals = GetFaceNormal(f);
                
                for(i = 0; i < order_f; ++i)
                {
                    FaceExp[f]->SetCoeff(i,inarray[i+cnt]);
                }
                cnt += order_f;
                
                FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(),
                                     FaceExp[f]->UpdatePhys());
                
                // Multiply by variable coefficient
                /// @TODO: Document this
//                StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
//                                                        StdRegions::eVarCoeffD11,
//                                                        StdRegions::eVarCoeffD22};
//                StdRegions::VarCoeffMap::const_iterator x;
//                Array<OneD, NekDouble> varcoeff_work(nquad_e);

//                if ((x = varcoeffs.find(VarCoeff[dir])) != varcoeffs.end())
//                {
//                    GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
//                    Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
//                }

                Vmath::Vmul(nquad_f,normals[dir],1,
                            FaceExp[f]->GetPhys(),1,
                            FaceExp[f]->UpdatePhys(),1);
                
                if (m_negatedNormals[f])
                {
                    Vmath::Neg(nquad_f, FaceExp[f]->UpdatePhys(), 1);
                }

                AddFaceBoundaryInt(f,FaceExp[f],outarray,varcoeffs);
            }
        }


        /**
         * For a given face add the \tilde{F}_1j contributions
         */
        void Expansion3D::v_AddFaceBoundaryInt(const int face,
                                              StdRegions::StdExpansionSharedPtr &FaceExp,
                                              Array <OneD,NekDouble > &outarray,
                                              const StdRegions::VarCoeffMap &varcoeffs)
        {
            int i;
            int order_f = FaceExp->GetNcoeffs();
            int nquad_f = FaceExp->GetNumPoints(0)*FaceExp->GetNumPoints(1);
            Array<OneD, NekDouble> coeff(order_f);

            StdRegions::IndexMapKey ikey(
                StdRegions::eFaceToElement, DetExpansionType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetFaceOrient(face));
            StdRegions::IndexMapValuesSharedPtr map = 
                StdExpansion::GetIndexMap(ikey);

//            StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
//                                                    StdRegions::eVarCoeffD11,
//                                                    StdRegions::eVarCoeffD22};
//            StdRegions::VarCoeffMap::const_iterator x;
//            Array<OneD, NekDouble> varcoeff_work(nquad_e);
//
///// @TODO Variable coeffs
//            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
//            {
//                GetPhysEdgeVarCoeffsFromElement(edge,EdgeExp,x->second,varcoeff_work);
//                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp->GetPhys(),1,EdgeExp->UpdatePhys(),1);
//            }

            FaceExp->IProductWRTBase(FaceExp->GetPhys(),coeff);

            // add data to out array
            for(i = 0; i < order_f; ++i)
            {
                outarray[(*map)[i].index] += (*map)[i].sign*coeff[i];
            }
        }

        /**
         * @brief Align face orientation with the geometry orientation.
         */	
        void Expansion3D::SetFaceToGeomOrientation(
            const int face, Array<OneD, NekDouble> &inout)
        {
            int j,k;
            int nface = GetFaceNcoeffs(face);
            Array<OneD, NekDouble> f_in(nface);
            Vmath::Vcopy(nface,&inout[0],1,&f_in[0],1);
            
            // retreiving face to element map for standard face orientation and
            // for actual face orientation
            StdRegions::IndexMapKey ikey1(
                StdRegions::eFaceToElement, DetExpansionType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, StdRegions::eDir1FwdDir1_Dir2FwdDir2);
            StdRegions::IndexMapValuesSharedPtr map1 = 
                StdExpansion::GetIndexMap(ikey1);
            StdRegions::IndexMapKey ikey2(
                StdRegions::eFaceToElement, DetExpansionType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetFaceOrient(face));
            StdRegions::IndexMapValuesSharedPtr map2 = 
                StdExpansion::GetIndexMap(ikey2);
            
            ASSERTL1((*map1).num_elements() == (*map2).num_elements(),
                     "There is an error with the GetFaceToElementMap");
            
            for(j = 0; j < (*map1).num_elements(); ++j)
            {
                // j = index in the standard orientation
                for(k = 0; k < (*map2).num_elements(); ++k)
                {
                    // k = index in the actual orientation
                    if((*map1)[j].index == (*map2)[k].index && k != j)
                    {
                        inout[k] = f_in[j];
                        //checking if sign is changing
                        if((*map1)[j].sign != (*map2)[k].sign)
                            inout[k] *= -1.0;
                        break;
                    }
                }
            }
        }
        
        /**
         * @brief Align trace orientation with the geometry orientation.
         */	
        void Expansion3D::SetTraceToGeomOrientation(Array<OneD, NekDouble> &inout)
        {
            int i,cnt = 0;
            int nfaces = GetNfaces();

            Array<OneD, NekDouble> f_tmp;
            
            for(i = 0; i < nfaces; ++i)
            {
                SetFaceToGeomOrientation(i, f_tmp = inout + cnt);
                cnt += GetFaceNcoeffs(i);
            }
        }

        /**
         * Computes matrices needed for the HDG formulation. References to
         * equations relate to the following paper (with a suitable changes in
         * formulation to adapt to 3D):
         *   R. M. Kirby, S. J. Sherwin, B. Cockburn, To CG or to HDG: A
         *   Comparative Study, J. Sci. Comp P1-30
         *   DOI 10.1007/s10915-011-9501-7
         *   NOTE: VARIABLE COEFFICIENTS CASE IS NOT IMPLEMENTED
         */
        DNekMatSharedPtr Expansion3D::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            //Variable coefficients are not implemented/////////
            ASSERTL1(!mkey.HasVarCoeff(StdRegions::eVarCoeffD00),
                     "Matrix construction is not implemented for variable "
                     "coefficients at the moment");
            ////////////////////////////////////////////////////
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
                // (Z^e)^{-1} (Eqn. 33, P22)
                case StdRegions::eHybridDGHelmholtz:
                {
                    ASSERTL1(IsBoundaryInteriorExpansion(),
                             "HybridDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");
                    
                    int       i,j,k;
                    NekDouble lambdaval = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    NekDouble tau       = mkey.GetConstFactor(StdRegions::eFactorTau);
                    int       ncoeffs   = GetNcoeffs();
                    int       nfaces    = GetNfaces();

                    Array<OneD,unsigned int> fmap;
                    Array<OneD,int> sign;
                    ExpansionSharedPtr FaceExp;
                    StdRegions::StdExpansionSharedPtr FaceExp2;

                    int order_f, coordim = GetCoordim();
                    DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;
                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    StdRegions::VarCoeffType Coeffs[3] = {StdRegions::eVarCoeffD00,
                                                            StdRegions::eVarCoeffD11,
                                                            StdRegions::eVarCoeffD22};

                    for(i=0;  i < coordim; ++i)
                    {
                        DNekScalMat &Dmat = *GetLocMatrix(DerivType[i]);
                        Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        
                        /*
                        if(mkey.HasVarCoeff(Coeffs[i]))
                        {
                            MatrixKey DmatkeyL(DerivType[i], DetExpansionType(), *this, 
                                               StdRegions::NullConstFactorMap,
                                               mkey.GetVarCoeffAsMap(Coeffs[i]));
                            MatrixKey DmatkeyR(DerivType[i], DetExpansionType(), *this);

                            DNekScalMat &DmatL = *GetLocMatrix(DmatkeyL);
                            DNekScalMat &DmatR = *GetLocMatrix(DmatkeyR);
                            Mat = Mat + DmatL*invMass*Transpose(DmatR);
                        }
                        else
                        {
                            DNekScalMat &Dmat = *GetLocMatrix(DerivType[i]);
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }
                        */
                    }

                    // Add Mass Matrix Contribution for Helmholtz problem
                    DNekScalMat  &Mass = *GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;

                    // Add tau*E_l using elemental mass matrices on each edge
                    for(i = 0; i < nfaces; ++i)
                    {
                        FaceExp = GetFaceExp(i);
                        order_f = FaceExp->GetNcoeffs();  
                        StdRegions::IndexMapKey ikey(
                            StdRegions::eFaceToElement, DetExpansionType(), 
                            GetBasisNumModes(0), GetBasisNumModes(1), 
                            GetBasisNumModes(2), i, GetFaceOrient(i));
                        StdRegions::IndexMapValuesSharedPtr map = 
                            StdExpansion::GetIndexMap(ikey);

                        // @TODO: Document
                        /*
                        StdRegions::VarCoeffMap edgeVarCoeffs;
                        if (mkey.HasVarCoeff(StdRegions::eVarCoeffD00))
                        {
                            Array<OneD, NekDouble> mu(nq);
                            GetPhysEdgeVarCoeffsFromElement(
                                i, EdgeExp2, 
                                mkey.GetVarCoeff(StdRegions::eVarCoeffD00), mu);
                            edgeVarCoeffs[StdRegions::eVarCoeffMass] = mu;
                        }
                        DNekScalMat &eMass = *EdgeExp->GetLocMatrix(
                            StdRegions::eMass, 
                            StdRegions::NullConstFactorMap, edgeVarCoeffs);
                        */

                        DNekScalMat &eMass = *FaceExp->GetLocMatrix(StdRegions::eMass);

                        for(j = 0; j < order_f; ++j)
                        {
                            for(k = 0; k < order_f; ++k)
                            {
                                Mat((*map)[j].index,(*map)[k].index) += 
                                    tau*(*map)[j].sign*(*map)[k].sign*eMass(j,k);
                            }
                        }
                    }
                    break;
                }

                // U^e (P22)
                case StdRegions::eHybridDGLamToU:
                {
                    int i,j,k;
                    int nbndry = NumDGBndryCoeffs();
                    int ncoeffs = GetNcoeffs();
                    int nfaces  = GetNfaces();
                    NekDouble lambdaval = mkey.GetConstFactor(StdRegions::eFactorLambda);
                    NekDouble tau       = mkey.GetConstFactor(StdRegions::eFactorTau);
                    
                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    StdRegions::StdExpansionSharedPtr FaceExp;
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Umat = *returnval;
                    
                    // Z^e matrix
                    MatrixKey newkey(StdRegions::eInvHybridDGHelmholtz, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat  &invHmat = *GetLocMatrix(newkey);

                    Array<OneD,unsigned int> fmap;
                    Array<OneD,int> sign;
                    
                    //alternative way to add boundary terms contribution
                    int bndry_cnt = 0;
                    for(i = 0; i < nfaces; ++i)
                    {
                        FaceExp = GetFaceExp(i);//temporary, need to rewrite AddHDGHelmholtzFaceTerms
                        int nface = GetFaceNcoeffs(i);
                        Array<OneD, NekDouble> face_lambda(nface);
			
                        const Array<OneD, const Array<OneD, NekDouble> > normals
                            = GetFaceNormal(i);

                        //cout << endl << "face #" << i;
                        //cout << endl << "nquad_f " << FaceExp[i]->GetNumPoints(0)*FaceExp[i]->GetNumPoints(1);
                        //cout << endl << "normals[0] " <<  normals[0].num_elements();
                        //cout << endl << "normals[1] " <<  normals[1].num_elements();
                        //cout << endl << "normals[2] " <<  normals[2].num_elements();

                        for(j = 0; j < nface; ++j)
                        {
                            Vmath::Zero(nface,&face_lambda[0],1);
                            Vmath::Zero(ncoeffs,&f[0],1);
                            face_lambda[j] = 1.0;

                            SetFaceToGeomOrientation(i, face_lambda);
							
                            Vmath::Vcopy(nface, face_lambda, 1, FaceExp->UpdateCoeffs(), 1);

                            FaceExp->BwdTrans(FaceExp->GetCoeffs(), FaceExp->UpdatePhys());

                            AddHDGHelmholtzFaceTerms(tau, i, FaceExp, mkey.GetVarCoeffs(), f);
							
                            Ulam = invHmat*F; // generate Ulam from lambda

                            // fill column of matrix
                            for(k = 0; k < ncoeffs; ++k)
                            {
                                Umat(k,bndry_cnt) = Ulam[k]; 
                            }
							
                            ++bndry_cnt;
                        }
                    }

                    //// Set up face expansions from local geom info
                    //for(i = 0; i < nfaces; ++i)
                    //{
                    //    FaceExp[i] = GetFaceExp(i);
                    //}
                    //
                    //// for each degree of freedom of the lambda space
                    //// calculate Umat entry 
                    //// Generate Lambda to U_lambda matrix 
                    //for(j = 0; j < nbndry; ++j)
                    //{
                    //    // standard basis vectors e_j
                    //    Vmath::Zero(nbndry,&lambda[0],1);
                    //    Vmath::Zero(ncoeffs,&f[0],1);
                    //    lambda[j] = 1.0;
                    //    
					//	//cout << Lambda;
                    //    SetTraceToGeomOrientation(lambda);
					//	//cout << Lambda << endl;
                    //    
                    //    // Compute F = [I   D_1 M^{-1}   D_2 M^{-1}] C e_j
                    //    AddHDGHelmholtzTraceTerms(tau, lambda, FaceExp, mkey.GetVarCoeffs(), f);
                    //    
                    //    // Compute U^e_j
                    //    Ulam = invHmat*F; // generate Ulam from lambda
                    //    
                    //    // fill column of matrix
                    //    for(k = 0; k < ncoeffs; ++k)
                    //    {
                    //        Umat(k,j) = Ulam[k]; 
                    //    }
                    //}
                }
                break;
            // Q_0, Q_1, Q_2 matrices (P23)
            // Each are a product of a row of Eqn 32 with the C matrix.
            // Rather than explicitly computing all of Eqn 32, we note each
            // row is almost a multiple of U^e, so use that as our starting
            // point.
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
                {
                    int i,j,k,dir;
                    int nbndry = NumDGBndryCoeffs();
                    //int nquad  = GetNumPoints(0);
                    int ncoeffs = GetNcoeffs();
                    int nfaces  = GetNfaces();

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    Array<OneD,StdRegions::StdExpansionSharedPtr>  FaceExp(nfaces);
                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    // Lambda to U matrix
                    MatrixKey lamToUkey(StdRegions::eHybridDGLamToU, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &lamToU = *GetLocMatrix(lamToUkey);

                    // Inverse mass matrix 
                    DNekScalMat &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    
                    for(i = 0; i < nfaces; ++i)
                    {
                        FaceExp[i] = GetFaceExp(i);
                    }

                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv0);
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv1);
                        break;
                    case StdRegions::eHybridDGLamToQ2:
                        dir = 2;
                        Dmat = GetLocMatrix(StdRegions::eWeakDeriv2);
                        break;
                    default:
                        ASSERTL0(false,"Direction not known");
                        break;
                    }
                
                    // for each degree of freedom of the lambda space
                    // calculate Qmat entry 
                    // Generate Lambda to Q_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        lambda[j] = 1.0;
                        
                        // for lambda[j] = 1 this is the solution to ulam
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Ulam[k] = lamToU(k,j);
                        }
                        
                        // -D^T ulam
                        Vmath::Neg(ncoeffs,&ulam[0],1);
                        F = Transpose(*Dmat)*Ulam; 
                        
                        SetTraceToGeomOrientation(lambda);
                        
                        // Add the C terms resulting from the I's on the
                        // diagonals of Eqn 32
                        AddNormTraceInt(dir,lambda,FaceExp,f,mkey.GetVarCoeffs());
                        
                        // finally multiply by inverse mass matrix
                        Ulam = invMass*F; 
                        
                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*ncoeffs,1);
                    }
                }
                break;            
            // Matrix K (P23)
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int i,j,f,cnt;
                    int order_f, nquad_f;
                    int nbndry  = NumDGBndryCoeffs();
                    int coordim = GetCoordim();
                    int nfaces  = GetNfaces();
                    NekDouble tau = mkey.GetConstFactor(StdRegions::eFactorTau);

                    Array<OneD,NekDouble>       work, varcoeff_work;
                    Array<OneD,const Array<OneD, NekDouble> > normals; 
                    Array<OneD,StdRegions::StdExpansionSharedPtr>  FaceExp(nfaces);
                    Array<OneD, NekDouble> lam(nbndry); 
                    
                    Array<OneD,unsigned int>    fmap;
                    Array<OneD, int>            sign;
                    StdRegions::Orientation facedir;
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    DNekScalMatSharedPtr LamToQ[3];
                    
                    // Matrix to map Lambda to U
                    MatrixKey LamToUkey(StdRegions::eHybridDGLamToU, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LamToU = *GetLocMatrix(LamToUkey);

                    // Matrix to map Lambda to Q0
                    MatrixKey LamToQ0key(StdRegions::eHybridDGLamToQ0, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[0] = GetLocMatrix(LamToQ0key);
 
                    // Matrix to map Lambda to Q1
                    MatrixKey LamToQ1key(StdRegions::eHybridDGLamToQ1, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[1] = GetLocMatrix(LamToQ1key);

                    // Matrix to map Lambda to Q2
                    MatrixKey LamToQ2key(StdRegions::eHybridDGLamToQ2, DetExpansionType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[2] = GetLocMatrix(LamToQ2key);

                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < nfaces; ++i)
                    {
                        FaceExp[i] = GetFaceExp(i);
                    }

                    // Set up matrix derived from <mu, Q_lam.n - \tau (U_lam - Lam) > 
                    for(i = 0; i < nbndry; ++i)
                    {
                        cnt = 0;
                        
                        Vmath::Zero(nbndry,lam,1);
                        lam[i] = 1.0;
                        SetTraceToGeomOrientation(lam);

                        for(f = 0; f < nfaces; ++f)
                        {
                            order_f = FaceExp[f]->GetNcoeffs();  
                            nquad_f = FaceExp[f]->GetNumPoints(0)*FaceExp[f]->GetNumPoints(1);    
                            normals = GetFaceNormal(f);
                            facedir = GetFaceOrient(f);
                            
                            work = Array<OneD,NekDouble>(nquad_f);
                            varcoeff_work = Array<OneD, NekDouble>(nquad_f);

                            StdRegions::IndexMapKey ikey(
                                StdRegions::eFaceToElement, DetExpansionType(), 
                                GetBasisNumModes(0), GetBasisNumModes(1), 
                                GetBasisNumModes(2), f, GetFaceOrient(f));
                            StdRegions::IndexMapValuesSharedPtr map = 
                                StdExpansion::GetIndexMap(ikey);

                            // @TODO Variable coefficients
                            /*
                            StdRegions::VarCoeffType VarCoeff[3] = {StdRegions::eVarCoeffD00,
                                                                    StdRegions::eVarCoeffD11,
                                                                    StdRegions::eVarCoeffD22};
                            const StdRegions::VarCoeffMap &varcoeffs = mkey.GetVarCoeffs();
                            StdRegions::VarCoeffMap::const_iterator x;
                            */

                            // Q0 * n0 (BQ_0 terms)
                            for(j = 0; j < order_f; ++j)
                            {
                                FaceExp[f]->SetCoeff(j,(*map)[j].sign*(*LamToQ[0])((*map)[j].index,i));
                            }
                            
                            FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(),
                                                 FaceExp[f]->UpdatePhys());

                            // @TODO Variable coefficients
                            // Multiply by variable coefficient
                            /*
                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */
          
                            Vmath::Vmul(nquad_f,normals[0],1,FaceExp[f]->GetPhys(),1,work,1);
                            
                            // Q1 * n1 (BQ_1 terms)
                            for(j = 0; j < order_f; ++j)
                            {
                                FaceExp[f]->SetCoeff(j,(*map)[j].sign*(*LamToQ[1])((*map)[j].index,i));
                            }
                            
                            FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(),
                                                 FaceExp[f]->UpdatePhys());

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[1])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Vvtvp(nquad_f,normals[1],1,
                                         FaceExp[f]->GetPhys(),1,
                                         work,1,work,1);
                            
                            // Q2 * n2 (BQ_2 terms)
                            for(j = 0; j < order_f; ++j)
                            {
                                FaceExp[f]->SetCoeff(j,(*map)[j].sign*(*LamToQ[2])((*map)[j].index,i));
                            }
                            
                            FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(),
                                                 FaceExp[f]->UpdatePhys());

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[2])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Vvtvp(nquad_f,normals[2],1,
                                         FaceExp[f]->GetPhys(),1,
                                         work,1,work,1);
                            
                            if (m_negatedNormals[f])
                            {
                                Vmath::Neg(nquad_f, work, 1);
                            }

                            // - tau (ulam - lam)
                            // Corresponds to the G and BU terms.
                            for(j = 0; j < order_f; ++j)
                            {
                                FaceExp[f]->SetCoeff(j,(*map)[j].sign*LamToU((*map)[j].index,i) - lam[cnt+j]);
                            }
                            
                            FaceExp[f]->BwdTrans(FaceExp[f]->GetCoeffs(),
                                                 FaceExp[f]->UpdatePhys());

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,FaceExp[f],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_f,varcoeff_work,1,FaceExp[f]->GetPhys(),1,FaceExp[f]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Svtvp(nquad_f,-tau,FaceExp[f]->GetPhys(),1,
                                         work,1,work,1);

                            // @TODO Add variable coefficients
                            FaceExp[f]->IProductWRTBase(work,FaceExp[f]->UpdateCoeffs());
                            
                            SetFaceToGeomOrientation(f, FaceExp[f]->UpdateCoeffs());

                            for(j = 0; j < order_f; ++j)
                            {
                                BndMat(cnt+j,i) = FaceExp[f]->GetCoeff(j);
                            }
                            
                            cnt += order_f;
                        }
                    }
                }
                break;
            default:
                ASSERTL0(false,"This matrix type cannot be generated from this class");
                break;
            }
            
            return returnval;
        }

        void Expansion3D::SetFaceExp(const int face, Expansion2DSharedPtr &f)
        {
            int nFaces = GetNfaces();
            ASSERTL1(face < nFaces, "Face is out of range.");
            if (m_faceExp.size() < nFaces)
            {
                m_faceExp.resize(nFaces);
            }
            m_faceExp[face] = f;
        }
        
        Expansion2DSharedPtr Expansion3D::GetFaceExp(const int face)
        {
            return m_faceExp[face].lock();
        }
        
        void Expansion3D::v_AddFaceNormBoundaryInt(
            const int                            face,
            StdRegions::StdExpansionSharedPtr   &FaceExp,
            const Array<OneD, const NekDouble>  &Fn,
                  Array<OneD,       NekDouble>  &outarray)
        {
            int i;
            
            StdRegions::IndexMapKey ikey(
                StdRegions::eFaceToElement, DetExpansionType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetFaceOrient(face));
            StdRegions::IndexMapValuesSharedPtr map = 
                StdExpansion::GetIndexMap(ikey);

            int order_e = (*map).num_elements(); // Order of the element
            int n_coeffs = FaceExp->GetCoeffs().num_elements(); // Order of the trace

            if(n_coeffs != order_e) // Going to orthogonal space
            {
                ASSERTL0(false, "Variable order not supported in 3D.");
            }
            else
            {
                FaceExp->IProductWRTBase(Fn,FaceExp->UpdateCoeffs());
                
                LocalRegions::Expansion2DSharedPtr locExp = 
                    boost::dynamic_pointer_cast<
                        LocalRegions::Expansion2D>(FaceExp);
                
                /*
                 * Coming into this routine, the velocity V will have been
                 * multiplied by the trace normals to give the input vector
                 * Vn. By convention, these normals are inwards facing for
                 * elements which have FaceExp as their right-adjacent face.
                 * This conditional statement therefore determines whether the
                 * normals must be negated, since the integral being performed
                 * here requires an outwards facing normal.
                 */ 
                if (locExp->GetRightAdjacentElementFace() != -1)
                {
                    if (locExp->GetRightAdjacentElementExp()->GetGeom3D()->GetGlobalID() 
                        == GetGeom3D()->GetGlobalID())
                    {
                        Vmath::Neg(order_e,FaceExp->UpdateCoeffs(),1);
                    }
                }
            }
            
            for(i = 0; i < order_e; ++i)
            {
                outarray[(*map)[i].index] += (*map)[i].sign*FaceExp->GetCoeff(i);
            }
        }

        void Expansion3D::v_AddRobinMassMatrix(
            const int face, 
            const Array<OneD, const NekDouble > &primCoeffs, 
            DNekMatSharedPtr &inoutmat)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");
            ASSERTL1(inoutmat->GetRows() == inoutmat->GetColumns(),
                     "Assuming that input matrix was square");
            
            int i,j;
            int id1,id2;
            Expansion2DSharedPtr faceExp = m_faceExp[face].lock();
            int order_f = faceExp->GetNcoeffs();
         
            Array<OneD, unsigned int> map;
            Array<OneD,          int> sign;
            
            StdRegions::VarCoeffMap varcoeffs;
            varcoeffs[StdRegions::eVarCoeffPrimative] = primCoeffs;

            StdRegions::ExpansionType expType = 
                faceExp->DetExpansionType();

            LocalRegions::MatrixKey mkey(
                StdRegions::eMass,
                expType, 
                *faceExp, 
                StdRegions::NullConstFactorMap, 
                varcoeffs);

            DNekScalMat &facemat = *faceExp->GetLocMatrix(mkey);

            // Now need to identify a map which takes the local face
            // mass matrix to the matrix stored in inoutmat;
            // This can currently be deduced from the size of the matrix
            
            // - if inoutmat.m_rows() == v_NCoeffs() it is a full
            //   matrix system
            
            // - if inoutmat.m_rows() == v_NumBndCoeffs() it is a
            //  boundary CG system

            // - if inoutmat.m_rows() == v_NumDGBndCoeffs() it is a
            //  trace DG system; still needs implementing.
            int rows = inoutmat->GetRows();

            if (rows == GetNcoeffs())
            {
                GetFaceToElementMap(face,GetFaceOrient(face),map,sign);
            }
            else if(rows == NumBndryCoeffs())
            {
                int nbndry = NumBndryCoeffs();
                Array<OneD,unsigned int> bmap(nbndry);

                GetFaceToElementMap(face,GetFaceOrient(face),map,sign);
                GetBoundaryMap(bmap);
                
                for(i = 0; i < order_f; ++i)
                {
                    for(j = 0; j < nbndry; ++j)
                    {
                        if(map[i] == bmap[j])
                        {
                            map[i] = j;
                            break;
                        }
                    }
                    ASSERTL1(j != nbndry,"Did not find number in map");
                }
            }
            else if (rows == NumDGBndryCoeffs())
            {
                // possibly this should be a separate method
                int cnt = 0; 
                map  = Array<OneD, unsigned int> (order_f);
                sign = Array<OneD,          int> (order_f,1);
                
                StdRegions::IndexMapKey ikey1(
                    StdRegions::eFaceToElement, DetExpansionType(), 
                    GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                    face, GetFaceOrient(face));
                StdRegions::IndexMapValuesSharedPtr map1 = 
                    StdExpansion::GetIndexMap(ikey1);
                StdRegions::IndexMapKey ikey2(
                    StdRegions::eFaceToElement, DetExpansionType(), 
                    GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                    face, StdRegions::eDir1FwdDir1_Dir2FwdDir2);
                StdRegions::IndexMapValuesSharedPtr map2 = 
                    StdExpansion::GetIndexMap(ikey2);
                
                ASSERTL1((*map1).num_elements() == (*map2).num_elements(),
                         "There is an error with the GetFaceToElementMap");
                
                for (i = 0; i < face; ++i)
                {
                    cnt += GetFaceNcoeffs(i);
                }

                for(i = 0; i < (*map1).num_elements(); ++i)
                {
                    int idx = -1;

                    for(j = 0; j < (*map2).num_elements(); ++j)
                    {
                        if((*map1)[i].index == (*map2)[j].index)
                        {
                            idx = j;
                            break;
                        }
                    }
                    
                    ASSERTL2(idx >= 0, "Index not found");
                    map [i] = idx + cnt;
                    sign[i] = (*map2)[idx].sign;
                }
            }
            else
            {
                ASSERTL0(false,"Could not identify matrix type from dimension");
            }

            for(i = 0; i < order_f; ++i)
            {
                id1 = map[i];
                for(j = 0; j < order_f; ++j)
                {
                    id2 = map[j];
                    (*inoutmat)(id1,id2) += facemat(i,j)*sign[i]*sign[j];
                }
            }
        }
    } //end of namespace
} //end of namespace
