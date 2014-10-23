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

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <SpatialDomains/Geometry3D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        //  evaluate additional terms in HDG face. Note that this assumes that
        // edges are unpacked into local cartesian order. 
        void Expansion3D::AddHDGHelmholtzFaceTerms(
            const NekDouble                    tau,
            const int                          face,
            Array<OneD, NekDouble>            &facePhys,
            const StdRegions::VarCoeffMap     &varcoeffs,
            Array<OneD, NekDouble>            &outarray)
        {
            ExpansionSharedPtr FaceExp = GetFaceExp(face);
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
                StdRegions::eFaceToElement, DetShapeType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetForient(face));
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
            FaceExp->IProductWRTBase(facePhys, outcoeff);

            for(i = 0; i < order_f; ++i)
            {
                outarray[(*map)[i].index] += (*map)[i].sign*tau*outcoeff[i];
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
                Vmath::Vmul(nquad_f, normals[n], 1, facePhys, 1, inval, 1);

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
                
                FaceExp->IProductWRTBase(inval, outcoeff);

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
        void Expansion3D::AddNormTraceInt(
            const int                        dir,
            Array<OneD, const NekDouble>    &inarray,
            Array<OneD, ExpansionSharedPtr> &FaceExp,
            Array<OneD, NekDouble>          &outarray,
            const StdRegions::VarCoeffMap   &varcoeffs)
        {
            int i,f,cnt;
            int order_f,nquad_f;
            int nfaces = GetNfaces();

            cnt = 0;
            for(f = 0; f < nfaces; ++f)
            {
                order_f = FaceExp[f]->GetNcoeffs();
                nquad_f = FaceExp[f]->GetNumPoints(0)*FaceExp[f]->GetNumPoints(1);

                const Array<OneD, const Array<OneD, NekDouble> > &normals = GetFaceNormal(f);
                Array<OneD, NekDouble> faceCoeffs(order_f);
                Array<OneD, NekDouble> facePhys  (nquad_f);

                for(i = 0; i < order_f; ++i)
                {
                    faceCoeffs[i] = inarray[i+cnt];
                }
                cnt += order_f;

                FaceExp[f]->BwdTrans(faceCoeffs, facePhys);

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

                Vmath::Vmul(nquad_f, normals[dir], 1, facePhys, 1, facePhys, 1);
                
                if (m_negatedNormals[f])
                {
                    Vmath::Neg(nquad_f, facePhys, 1);
                }

                AddFaceBoundaryInt(f, FaceExp[f], facePhys, outarray, varcoeffs);
            }
        }

        //shorter version of the above (coefficients are already set for faces)
        void Expansion3D::AddNormTraceInt(
            const int                             dir,
            Array<OneD, ExpansionSharedPtr>      &FaceExp,
            Array<OneD, Array<OneD, NekDouble> > &faceCoeffs,
            Array<OneD, NekDouble>               &outarray)
        {
            int f, cnt;
            int order_f, nquad_f;
            int nfaces = GetNfaces();

            cnt = 0;
            for(f = 0; f < nfaces; ++f)
            {
                order_f = FaceExp[f]->GetNcoeffs();
                nquad_f = FaceExp[f]->GetNumPoints(0)*FaceExp[f]->GetNumPoints(1);

                const Array<OneD, const Array<OneD, NekDouble> > &normals = GetFaceNormal(f);
                Array<OneD, NekDouble> facePhys(nquad_f);
                
                cnt += order_f;
                
                FaceExp[f]->BwdTrans(faceCoeffs[f], facePhys);
                
                Vmath::Vmul(nquad_f, normals[dir], 1, facePhys, 1, facePhys, 1);
                
                if (m_negatedNormals[f])
                {
                    Vmath::Neg(nquad_f, facePhys, 1);
                }

                AddFaceBoundaryInt(f, FaceExp[f], facePhys, outarray);
            }
        }

        /**
         * For a given face add the \tilde{F}_1j contributions
         */
        void Expansion3D::AddFaceBoundaryInt(
            const int                      face,
            ExpansionSharedPtr            &FaceExp,
            Array<OneD, NekDouble>        &facePhys,
            Array<OneD, NekDouble>        &outarray,
            const StdRegions::VarCoeffMap &varcoeffs)
        {
            int i;
            int order_f = FaceExp->GetNcoeffs();
            Array<OneD, NekDouble> coeff(order_f);

            StdRegions::IndexMapKey ikey(
                StdRegions::eFaceToElement, DetShapeType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetForient(face));
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

            FaceExp->IProductWRTBase(facePhys, coeff);

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
                StdRegions::eFaceToElement, DetShapeType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, StdRegions::eDir1FwdDir1_Dir2FwdDir2);
            StdRegions::IndexMapValuesSharedPtr map1 = 
                StdExpansion::GetIndexMap(ikey1);
            StdRegions::IndexMapKey ikey2(
                StdRegions::eFaceToElement, DetShapeType(), 
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetForient(face));
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
                    ExpansionSharedPtr FaceExp2;

                    int order_f, coordim = GetCoordim();
                    DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;
                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    // StdRegions::VarCoeffType Coeffs[3] = {StdRegions::eVarCoeffD00,
                    //                                       StdRegions::eVarCoeffD11,
                    //                                       StdRegions::eVarCoeffD22};

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
                            StdRegions::eFaceToElement, DetShapeType(), 
                            GetBasisNumModes(0), GetBasisNumModes(1), 
                            GetBasisNumModes(2), i, GetForient(i));
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
                    int       i,j,k;
                    int       nbndry  = NumDGBndryCoeffs();
                    int       ncoeffs = GetNcoeffs();
                    int       nfaces  = GetNfaces();
                    NekDouble tau     = mkey.GetConstFactor(StdRegions::eFactorTau);
                    
                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    ExpansionSharedPtr FaceExp;
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Umat = *returnval;
                    
                    // Z^e matrix
                    MatrixKey newkey(StdRegions::eInvHybridDGHelmholtz, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
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

                        for(j = 0; j < nface; ++j)
                        {
                            Vmath::Zero(nface,&face_lambda[0],1);
                            Vmath::Zero(ncoeffs,&f[0],1);
                            face_lambda[j] = 1.0;

                            SetFaceToGeomOrientation(i, face_lambda);

                            Array<OneD, NekDouble> tmp(FaceExp->GetTotPoints());
                            FaceExp->BwdTrans(face_lambda, tmp);
                            AddHDGHelmholtzFaceTerms(tau, i, tmp, mkey.GetVarCoeffs(), f);
							
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
                    Array<OneD, ExpansionSharedPtr>  FaceExp(nfaces);
                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    // Lambda to U matrix
                    MatrixKey lamToUkey(StdRegions::eHybridDGLamToU, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
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
                    int nfaces  = GetNfaces();
                    NekDouble tau = mkey.GetConstFactor(StdRegions::eFactorTau);

                    Array<OneD,NekDouble>       work, varcoeff_work;
                    Array<OneD,const Array<OneD, NekDouble> > normals; 
                    Array<OneD, ExpansionSharedPtr>  FaceExp(nfaces);
                    Array<OneD, NekDouble> lam(nbndry); 
                    
                    Array<OneD,unsigned int>    fmap;
                    Array<OneD, int>            sign;
                    
                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    DNekScalMatSharedPtr LamToQ[3];
                    
                    // Matrix to map Lambda to U
                    MatrixKey LamToUkey(StdRegions::eHybridDGLamToU, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat &LamToU = *GetLocMatrix(LamToUkey);

                    // Matrix to map Lambda to Q0
                    MatrixKey LamToQ0key(StdRegions::eHybridDGLamToQ0, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[0] = GetLocMatrix(LamToQ0key);
 
                    // Matrix to map Lambda to Q1
                    MatrixKey LamToQ1key(StdRegions::eHybridDGLamToQ1, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    LamToQ[1] = GetLocMatrix(LamToQ1key);

                    // Matrix to map Lambda to Q2
                    MatrixKey LamToQ2key(StdRegions::eHybridDGLamToQ2, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
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
                            
                            work = Array<OneD,NekDouble>(nquad_f);
                            varcoeff_work = Array<OneD, NekDouble>(nquad_f);

                            StdRegions::IndexMapKey ikey(
                                StdRegions::eFaceToElement, DetShapeType(), 
                                GetBasisNumModes(0), GetBasisNumModes(1), 
                                GetBasisNumModes(2), f, GetForient(f));
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
                            Array<OneD, NekDouble> faceCoeffs(order_f);
                            Array<OneD, NekDouble> facePhys  (nquad_f);
                            for(j = 0; j < order_f; ++j)
                            {
                                faceCoeffs[j] = (*map)[j].sign*(*LamToQ[0])((*map)[j].index,i);
                            }
                            
                            FaceExp[f]->BwdTrans(faceCoeffs, facePhys);

                            // @TODO Variable coefficients
                            // Multiply by variable coefficient
                            /*
                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */
          
                            Vmath::Vmul(nquad_f, normals[0], 1, facePhys, 1, work, 1);
                            
                            // Q1 * n1 (BQ_1 terms)
                            for(j = 0; j < order_f; ++j)
                            {
                                faceCoeffs[j] = (*map)[j].sign*(*LamToQ[1])((*map)[j].index,i);
                            }
                            
                            FaceExp[f]->BwdTrans(faceCoeffs, facePhys);

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[1])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Vvtvp(nquad_f, normals[1], 1, facePhys, 1, work, 1, work, 1);
                            
                            // Q2 * n2 (BQ_2 terms)
                            for(j = 0; j < order_f; ++j)
                            {
                                faceCoeffs[j] = (*map)[j].sign*(*LamToQ[2])((*map)[j].index,i);
                            }
                            
                            FaceExp[f]->BwdTrans(faceCoeffs, facePhys);

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[2])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,EdgeExp[e],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_e,varcoeff_work,1,EdgeExp[e]->GetPhys(),1,EdgeExp[e]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Vvtvp(nquad_f, normals[2], 1, facePhys, 1,
                                                  work,       1, work,     1);

                            if (m_negatedNormals[f])
                            {
                                Vmath::Neg(nquad_f, work, 1);
                            }

                            // - tau (ulam - lam)
                            // Corresponds to the G and BU terms.
                            for(j = 0; j < order_f; ++j)
                            {
                                faceCoeffs[j] = (*map)[j].sign*LamToU((*map)[j].index,i) - lam[cnt+j];
                            }
                            
                            FaceExp[f]->BwdTrans(faceCoeffs, facePhys);

                            // @TODO Variable coefficients
                            // Multiply by variable coefficients
                            /*
                            if ((x = varcoeffs.find(VarCoeff[0])) != varcoeffs.end())
                            {
                                GetPhysEdgeVarCoeffsFromElement(e,FaceExp[f],x->second,varcoeff_work);
                                Vmath::Vmul(nquad_f,varcoeff_work,1,FaceExp[f]->GetPhys(),1,FaceExp[f]->UpdatePhys(),1);
                            }
                            */

                            Vmath::Svtvp(nquad_f, -tau, facePhys, 1,
                                         work, 1, work, 1);

                            // @TODO Add variable coefficients
                            FaceExp[f]->IProductWRTBase(work, faceCoeffs);
                            
                            SetFaceToGeomOrientation(f, faceCoeffs);

                            for(j = 0; j < order_f; ++j)
                            {
                                BndMat(cnt+j,i) = faceCoeffs[j];
                            }

                            cnt += order_f;
                        }
                    }
                }
                break;
            //HDG postprocessing
            case StdRegions::eInvLaplacianWithUnityMean:
                {
                    MatrixKey lapkey(StdRegions::eLaplacian, DetShapeType(), *this, mkey.GetConstFactors(), mkey.GetVarCoeffs());
                    DNekScalMat  &LapMat = *GetLocMatrix(lapkey);
                    
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(LapMat.GetRows(),LapMat.GetColumns());
                    DNekMatSharedPtr lmat = returnval;
                    
                    (*lmat) = LapMat;

                    // replace first column with inner product wrt 1
                    int nq = GetTotPoints();
                    Array<OneD, NekDouble> tmp(nq);
                    Array<OneD, NekDouble> outarray(m_ncoeffs);
                    Vmath::Fill(nq,1.0,tmp,1);
                    IProductWRTBase(tmp, outarray);

                    Vmath::Vcopy(m_ncoeffs,&outarray[0],1,
                                 &(lmat->GetPtr())[0],1);

                    //cout << endl << *lmat << endl;
                    lmat->Invert();
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
            const int                           face,
            const ExpansionSharedPtr           &FaceExp,
            const Array<OneD, const NekDouble> &Fn,
                  Array<OneD,       NekDouble> &outarray)
        {
            int i, j;
            
            /*
             * Coming into this routine, the velocity V will have been
             * multiplied by the trace normals to give the input vector Vn. By
             * convention, these normals are inwards facing for elements which
             * have FaceExp as their right-adjacent face.  This conditional
             * statement therefore determines whether the normals must be
             * negated, since the integral being performed here requires an
             * outwards facing normal.
             */
            if (m_requireNeg.size() == 0)
            {
                m_requireNeg.resize(GetNfaces());
                
                for (i = 0; i < GetNfaces(); ++i)
                {
                    m_requireNeg[i] = false;
                    if (m_negatedNormals[i])
                    {
                        m_requireNeg[i] = true;
                        continue;
                    }
                    
                    Expansion2DSharedPtr faceExp = m_faceExp[i].lock();

                    if (faceExp->GetRightAdjacentElementExp())
                    {
                        if (faceExp->GetRightAdjacentElementExp()->GetGeom3D()
                            ->GetGlobalID() == GetGeom3D()->GetGlobalID())
                        {
                            m_requireNeg[i] = true;
                        }
                    }
                }
            }

            StdRegions::IndexMapKey ikey(
                StdRegions::eFaceToElement, DetShapeType(),
                GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                face, GetForient(face));
            StdRegions::IndexMapValuesSharedPtr map =
            StdExpansion::GetIndexMap(ikey);

            int order_e  = (*map).num_elements(); // Order of the element
            int n_coeffs = FaceExp->GetNcoeffs();

            Array<OneD, NekDouble> faceCoeffs(n_coeffs);

            if (n_coeffs != order_e) // Going to orthogonal space
            {
                Array<OneD, NekDouble> coeff(n_coeffs);
                Array<OneD, NekDouble> array(n_coeffs);

                FaceExp->FwdTrans(Fn, faceCoeffs);

                int NumModesElementMax  = FaceExp->GetBasis(0)->GetNumModes();
                int NumModesElementMin  = m_base[0]->GetNumModes();

                FaceExp->ReduceOrderCoeffs(NumModesElementMin,
                                           faceCoeffs,
                                           faceCoeffs);

                StdRegions::StdMatrixKey masskey(
                    StdRegions::eMass, FaceExp->DetShapeType(), *FaceExp);
                FaceExp->MassMatrixOp(
                                      faceCoeffs,faceCoeffs,masskey);

                // Reorder coefficients for the lower degree face.
                int offset1 = 0, offset2 = 0;

                if (FaceExp->DetShapeType() == LibUtilities::eQuadrilateral)
                {
                    for (i = 0; i < NumModesElementMin; ++i)
                    {
                        for (j = 0; j < NumModesElementMin; ++j)
                        {
                            faceCoeffs[offset1+j] =
                            faceCoeffs[offset2+j];
                        }
                        offset1 += NumModesElementMin;
                        offset2 += NumModesElementMax;
                    }

                    // Extract lower degree modes. TODO: Check this is correct.
                    for (i = NumModesElementMin; i < NumModesElementMax; ++i)
                    {
                        for (j = NumModesElementMin; j < NumModesElementMax; ++j)
                        {
                            faceCoeffs[i*NumModesElementMax+j] = 0.0;
                        }
                    }
                }

                if (FaceExp->DetShapeType() == LibUtilities::eTriangle)
                {

                    // Reorder coefficients for the lower degree face.
                    int offset1 = 0, offset2 = 0;

                    for (i = 0; i < NumModesElementMin; ++i)
                    {
                        for (j = 0; j < NumModesElementMin-i; ++j)
                        {
                            faceCoeffs[offset1+j] =
                            faceCoeffs[offset2+j];
                        }
                        offset1 += NumModesElementMin-i;
                        offset2 += NumModesElementMax-i;
                    }
                }

            }
            else
            {
                FaceExp->IProductWRTBase(Fn, faceCoeffs);
            }

            if (m_requireNeg[face])
            {
                for (i = 0; i < order_e; ++i)
                {
                    outarray[(*map)[i].index] -= (*map)[i].sign * faceCoeffs[i];
                }
            }
            else
            {
                for (i = 0; i < order_e; ++i)
                {
                    outarray[(*map)[i].index] += (*map)[i].sign * faceCoeffs[i];
                }
            }
        }


        /**
         * @brief Evaluate coefficients of weak deriviative in the direction dir
         * given the input coefficicents incoeffs and the imposed boundary
         * values in EdgeExp (which will have its phys space updated).
         */
        void Expansion3D::v_DGDeriv(
            int                                   dir,
            const Array<OneD, const NekDouble>   &incoeffs,
            Array<OneD, ExpansionSharedPtr>      &FaceExp,
            Array<OneD, Array<OneD, NekDouble> > &faceCoeffs,
            Array<OneD, NekDouble>               &out_d)
        {
            int ncoeffs = GetNcoeffs();
            StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                   StdRegions::eWeakDeriv1,
                                                   StdRegions::eWeakDeriv2};

            DNekScalMat &InvMass = *GetLocMatrix(StdRegions::eInvMass);
            DNekScalMat &Dmat    = *GetLocMatrix(DerivType[dir]);

            Array<OneD, NekDouble> coeffs = incoeffs;
            DNekVec     Coeffs  (ncoeffs,coeffs, eWrapper);

            Coeffs = Transpose(Dmat)*Coeffs;
            Vmath::Neg(ncoeffs, coeffs,1);

            // Add the boundary integral including the relevant part of
            // the normal
            AddNormTraceInt(dir, FaceExp, faceCoeffs, coeffs);

            DNekVec Out_d (ncoeffs,out_d,eWrapper);

            Out_d  = InvMass*Coeffs;
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
            varcoeffs[StdRegions::eVarCoeffMass] = primCoeffs;

            LibUtilities::ShapeType shapeType = 
                faceExp->DetShapeType();

            LocalRegions::MatrixKey mkey(
                StdRegions::eMass,
                shapeType, 
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
                GetFaceToElementMap(face,GetForient(face),map,sign);
            }
            else if(rows == NumBndryCoeffs())
            {
                int nbndry = NumBndryCoeffs();
                Array<OneD,unsigned int> bmap(nbndry);

                GetFaceToElementMap(face,GetForient(face),map,sign);
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
                    StdRegions::eFaceToElement, DetShapeType(), 
                    GetBasisNumModes(0), GetBasisNumModes(1), GetBasisNumModes(2),
                    face, GetForient(face));
                StdRegions::IndexMapValuesSharedPtr map1 = 
                    StdExpansion::GetIndexMap(ikey1);
                StdRegions::IndexMapKey ikey2(
                    StdRegions::eFaceToElement, DetShapeType(), 
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

        DNekMatSharedPtr Expansion3D::v_BuildVertexMatrix(
            const DNekScalMatSharedPtr &r_bnd)
        {
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_vertexmatrix;

            int nVerts, vid1, vid2, vMap1, vMap2;
            NekDouble VertexValue;

            nVerts = GetNverts();

            m_vertexmatrix =
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nVerts, nVerts, 0.0, storage);
            DNekMat &VertexMat = (*m_vertexmatrix);

            for (vid1 = 0; vid1 < nVerts; ++vid1)
            {
                vMap1 = GetVertexMap(vid1,true);

                for (vid2 = 0; vid2 < nVerts; ++vid2)
                {
                    vMap2 = GetVertexMap(vid2,true);
                    VertexValue = (*r_bnd)(vMap1, vMap2);
                    VertexMat.SetValue(vid1, vid2, VertexValue);
                }
            }

            return m_vertexmatrix;
        }

        DNekMatSharedPtr Expansion3D::v_BuildTransformationMatrix(
            const DNekScalMatSharedPtr &r_bnd, 
            const StdRegions::MatrixType matrixType)
        {
            int nVerts, nEdges;
            int eid, fid, vid, n, i;

            int nBndCoeffs = NumBndryCoeffs();

            const SpatialDomains::Geometry3DSharedPtr &geom = GetGeom3D();

            // Get geometric information about this element
            nVerts = GetNverts();
            nEdges = GetNedges();

            /*************************************/
            /* Vetex-edge & vertex-face matrices */
            /*************************************/

            /**
             * The matrix component of \f$\mathbf{R}\f$ is given by \f[
             * \mathbf{R^{T}_{v}}=
             * -\mathbf{S}^{-1}_{ef,ef}\mathbf{S}^{T}_{v,ef}\f]
             *
             * For every vertex mode we extract the submatrices from statically
             * condensed matrix \f$\mathbf{S}\f$ corresponding to the coupling
             * between the attached edges and faces of a vertex
             * (\f$\mathbf{S_{ef,ef}}\f$). This matrix is then inverted and
             * multiplied by the submatrix representing the coupling between a
             * vertex and the attached edges and faces
             * (\f$\mathbf{S_{v,ef}}\f$).
             */

            int nmodes;
            int m;
            NekDouble VertexEdgeFaceValue;

            // The number of connected edges/faces is 3 (for all elements)
            int nConnectedEdges = 3;
            int nConnectedFaces = 3;

            // Location in the matrix
            Array<OneD, Array<OneD, unsigned int> >
                MatEdgeLocation(nConnectedEdges);
            Array<OneD, Array<OneD, unsigned int> >
                MatFaceLocation(nConnectedFaces);

            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_transformationmatrix;
            DNekMatSharedPtr m_transposedtransformationmatrix;

            m_transformationmatrix =
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nBndCoeffs, nBndCoeffs, 0.0, storage);
            m_transposedtransformationmatrix =
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nBndCoeffs, nBndCoeffs, 0.0, storage);

            DNekMat &R  = (*m_transformationmatrix);
            DNekMat &RT = (*m_transposedtransformationmatrix);

            // Build the vertex-edge/face transform matrix: This matrix is
            // constructed from the submatrices corresponding to the couping
            // between each vertex and the attached edges/faces
            for (vid = 0; vid < nVerts; ++vid)
            {
                // Row and column size of the vertex-edge/face matrix
                int efRow =
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 0)) +
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 1)) +
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 2)) +
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 0)) +
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 1)) +
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 2)) - 6;

                int nedgemodesconnected =
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 0)) +
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 1)) +
                    GetEdgeNcoeffs   (geom->GetVertexEdgeMap(vid, 2)) - 6;
                Array<OneD, unsigned int> edgemodearray(nedgemodesconnected);

                int nfacemodesconnected =
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 0)) +
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 1)) +
                    GetFaceIntNcoeffs(geom->GetVertexFaceMap(vid, 2));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);

                int offset = 0;

                // Create array of edge modes
                for (eid = 0; eid < nConnectedEdges; ++eid)
                {
                    MatEdgeLocation[eid] = GetEdgeInverseBoundaryMap(
                        geom->GetVertexEdgeMap(vid, eid));
                    nmodes = MatEdgeLocation[eid].num_elements();

                    if (nmodes)
                    {
                        Vmath::Vcopy(nmodes, &MatEdgeLocation[eid][0],
                                     1, &edgemodearray[offset], 1);
                    }

                    offset += nmodes;
                }

                offset = 0;

                // Create array of face modes
                for (fid = 0; fid < nConnectedFaces; ++fid)
                {
                    MatFaceLocation[fid] = GetFaceInverseBoundaryMap(
                        geom->GetVertexFaceMap(vid, fid));
                    nmodes = MatFaceLocation[fid].num_elements();

                    if (nmodes)
                    {
                        Vmath::Vcopy(nmodes, &MatFaceLocation[fid][0],
                                     1, &facemodearray[offset], 1);
                    }
                    offset += nmodes;
                }

                DNekMatSharedPtr m_vertexedgefacetransformmatrix =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        1, efRow, 0.0, storage);
                DNekMat &Sveft = (*m_vertexedgefacetransformmatrix);

                DNekMatSharedPtr m_vertexedgefacecoupling =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        1, efRow, 0.0, storage);
                DNekMat &Svef = (*m_vertexedgefacecoupling);

                // Vertex-edge coupling
                for (n = 0; n < nedgemodesconnected; ++n)
                {
                    // Matrix value for each coefficient location
                    VertexEdgeFaceValue = (*r_bnd)(GetVertexMap(vid),
                                                   edgemodearray[n]);

                    // Set the value in the vertex edge/face matrix
                    Svef.SetValue(0, n, VertexEdgeFaceValue);
                }

                // Vertex-face coupling
                for (n = 0; n < nfacemodesconnected; ++n)
                {
                    // Matrix value for each coefficient location
                    VertexEdgeFaceValue = (*r_bnd)(GetVertexMap(vid),
                                                   facemodearray[n]);

                    // Set the value in the vertex edge/face matrix
                    Svef.SetValue(0, n + nedgemodesconnected, VertexEdgeFaceValue);
                }

                /*
                 * Build the edge-face transform matrix: This matrix is
                 * constructed from the submatrices corresponding to the couping
                 * between the edges and faces on the attached faces/edges of a
                 * vertex.
                 */

                // Allocation of matrix to store edge/face-edge/face coupling
                DNekMatSharedPtr m_edgefacecoupling =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        efRow, efRow, 0.0, storage);
                DNekMat &Sefef = (*m_edgefacecoupling);

                NekDouble EdgeEdgeValue, FaceFaceValue;

                // Edge-edge coupling (S_{ee})
                for (m = 0; m < nedgemodesconnected; ++m)
                {
                    for (n = 0; n < nedgemodesconnected; ++n)
                    {
                        // Matrix value for each coefficient location
                        EdgeEdgeValue = (*r_bnd)(edgemodearray[n],
                                                 edgemodearray[m]);

                        // Set the value in the vertex edge/face matrix
                        Sefef.SetValue(n, m, EdgeEdgeValue);
                    }
                }

                // Face-face coupling (S_{ff})
                for (n = 0; n < nfacemodesconnected; ++n)
                {
                    for (m = 0; m < nfacemodesconnected; ++m)
                    {
                        // Matrix value for each coefficient location
                        FaceFaceValue = (*r_bnd)(facemodearray[n],
                                                 facemodearray[m]);
                        // Set the value in the vertex edge/face matrix
                        Sefef.SetValue(nedgemodesconnected + n,
                                       nedgemodesconnected + m, FaceFaceValue);
                    }
                }

                // Edge-face coupling (S_{ef} and trans(S_{ef}))
                for (n = 0; n < nedgemodesconnected; ++n)
                {
                    for (m = 0; m < nfacemodesconnected; ++m)
                    {
                        // Matrix value for each coefficient location
                        FaceFaceValue = (*r_bnd)(edgemodearray[n],
                                                 facemodearray[m]);

                        // Set the value in the vertex edge/face matrix
                        Sefef.SetValue(
                            n, nedgemodesconnected + m, FaceFaceValue);

                        // and transpose
                        Sefef.SetValue(
                            nedgemodesconnected + m, n, FaceFaceValue);
                    }
                }


                // Invert edge-face coupling matrix
                if (efRow)
                {
                    Sefef.Invert();

                    //R_{v}=-S_{v,ef}inv(S_{ef,ef})
                    Sveft = -Svef * Sefef;
                }

                // Populate R with R_{ve} components
                for (n = 0; n < edgemodearray.num_elements(); ++n)
                {
                    RT.SetValue(edgemodearray[n], GetVertexMap(vid),
                                Sveft(0, n));
                    R.SetValue(GetVertexMap(vid), edgemodearray[n],
                               Sveft(0, n));
                }

                // Populate R with R_{vf} components
                for (n = 0; n < facemodearray.num_elements(); ++n)
                {
                    RT.SetValue(facemodearray[n], GetVertexMap(vid),
                                Sveft(0, n + nedgemodesconnected));
                    R.SetValue(GetVertexMap(vid), facemodearray[n],
                               Sveft(0, n + nedgemodesconnected));
                }
            }

            /********************/
            /* edge-face matrix */
            /********************/

            /* 
             * The matrix component of \f$\mathbf{R}\f$ is given by \f[
             * \mathbf{R^{T}_{ef}}=-\mathbf{S}^{-1}_{ff}\mathbf{S}^{T}_{ef}\f]
             *
             * For each edge extract the submatrices from statically condensed
             * matrix \f$\mathbf{S}\f$ corresponding to inner products of modes
             * on the two attached faces within themselves as well as the
             * coupling matrix between the two faces
             * (\f$\mathbf{S}_{ff}\f$). This matrix of face coupling is then
             * inverted and multiplied by the submatrices of corresponding to
             * the coupling between the edge and attached faces
             * (\f$\mathbf{S}_{ef}\f$).
             */

            NekDouble EdgeFaceValue, FaceFaceValue;
            int efCol, efRow, nedgemodes;

            // Number of attached faces is always 2
            nConnectedFaces = 2;

            // Location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> >
                (nEdges);
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> >
                (nConnectedFaces);

            // Build the edge/face transform matrix: This matrix is constructed
            // from the submatrices corresponding to the couping between a
            // specific edge and the two attached faces.
            for (eid = 0; eid < nEdges; ++eid)
            {
                // Row and column size of the vertex-edge/face matrix
                efCol = GetFaceIntNcoeffs(geom->GetEdgeFaceMap(eid, 0)) +
                    GetFaceIntNcoeffs(geom->GetEdgeFaceMap(eid, 1));
                efRow = GetEdgeNcoeffs(eid) - 2;

                // Edge-face coupling matrix
                DNekMatSharedPtr m_efedgefacecoupling =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        efRow, efCol, 0.0, storage);
                DNekMat &Mef = (*m_efedgefacecoupling);

                // Face-face coupling matrix
                DNekMatSharedPtr m_effacefacecoupling =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        efCol, efCol, 0.0, storage);
                DNekMat &Meff = (*m_effacefacecoupling);

                // Edge-face transformation matrix
                DNekMatSharedPtr m_edgefacetransformmatrix =
                    MemoryManager<DNekMat>::AllocateSharedPtr(
                        efRow, efCol, 0.0, storage);
                DNekMat &Meft = (*m_edgefacetransformmatrix);

                int nfacemodesconnected =
                    GetFaceIntNcoeffs(geom->GetEdgeFaceMap(eid, 0)) +
                    GetFaceIntNcoeffs(geom->GetEdgeFaceMap(eid, 1));
                Array<OneD, unsigned int>
                    facemodearray(nfacemodesconnected);

                // Create array of edge modes
                Array<OneD, unsigned int> inedgearray
                    = GetEdgeInverseBoundaryMap(eid);
                nedgemodes = GetEdgeNcoeffs(eid) - 2;
                Array<OneD, unsigned int> edgemodearray(nedgemodes);

                if (nedgemodes)
                {
                    Vmath::Vcopy(nedgemodes, &inedgearray[0],
                                 1, &edgemodearray[0], 1);
                }

                int offset = 0;

                // Create array of face modes
                for (fid = 0; fid < nConnectedFaces; ++fid)
                {
                    MatFaceLocation[fid] = GetFaceInverseBoundaryMap(
                        geom->GetEdgeFaceMap(eid, fid));
                    nmodes = MatFaceLocation[fid].num_elements();

                    if (nmodes)
                    {
                        Vmath::Vcopy(nmodes, &MatFaceLocation[fid][0],
                                     1, &facemodearray[offset], 1);
                    }
                    offset += nmodes;
                }

                // Edge-face coupling
                for (n = 0; n < nedgemodes; ++n)
                {
                    for (m = 0; m < nfacemodesconnected; ++m)
                    {
                        // Matrix value for each coefficient location
                        EdgeFaceValue = (*r_bnd)(edgemodearray[n],
                                                 facemodearray[m]);

                        // Set the value in the edge/face matrix
                        Mef.SetValue(n, m, EdgeFaceValue);
                    }
                }

                // Face-face coupling
                for (n = 0; n < nfacemodesconnected; ++n)
                {
                    for (m = 0; m < nfacemodesconnected; ++m)
                    {
                        // Matrix value for each coefficient location
                        FaceFaceValue = (*r_bnd)(facemodearray[n],
                                                 facemodearray[m]);

                        // Set the value in the vertex edge/face matrix
                        Meff.SetValue(n, m, FaceFaceValue);
                    }
                }

                if (efCol)
                {
                    // Invert edge-face coupling matrix
                    Meff.Invert();

                    // trans(R_{ef})=-S_{ef}*(inv(S_{ff})
                    Meft = -Mef * Meff;
                }

                //Populate transformation matrix with Meft
                for (n = 0; n < Meft.GetRows(); ++n)
                {
                    for (m = 0; m < Meft.GetColumns(); ++m)
                    {
                        R.SetValue(edgemodearray[n], facemodearray[m],
                                   Meft(n, m));
                        RT.SetValue(facemodearray[m], edgemodearray[n],
                                    Meft(n, m));
                    }
                }
            }

            for (i = 0; i < R.GetRows(); ++i)
            {
                R.SetValue(i, i, 1.0);
                RT.SetValue(i, i, 1.0);
            }

            if ((matrixType == StdRegions::ePreconR)||
                (matrixType == StdRegions::ePreconRMass))
            {
                return m_transformationmatrix;
            }
            else if ((matrixType == StdRegions::ePreconRT)||
                     (matrixType == StdRegions::ePreconRTMass))
            {
                return m_transposedtransformationmatrix;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal, "unkown matrix type" );
                return NullDNekMatSharedPtr;
            }
	}

        /**
	 * \brief Build inverse and inverse transposed transformation matrix:
	 * \f$\mathbf{R^{-1}}\f$ and \f$\mathbf{R^{-T}}\f$
	 *
	 * \f\mathbf{R^{-T}}=[\left[\begin{array}{ccc} \mathbf{I} &
	 * -\mathbf{R}_{ef} & -\mathbf{R}_{ve}+\mathbf{R}_{ve}\mathbf{R}_{vf} \\
	 *  0 & \mathbf{I} & \mathbf{R}_{ef} \\
	 *  0 & 0 & \mathbf{I}} \end{array}\right]\f]
	 */
        DNekMatSharedPtr Expansion3D::v_BuildInverseTransformationMatrix(
            const DNekScalMatSharedPtr & m_transformationmatrix)
	{
            int i, j, n, eid = 0, fid = 0;
            int nCoeffs = NumBndryCoeffs();
            NekDouble MatrixValue;
            DNekScalMat &R = (*m_transformationmatrix);

            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;

            // Inverse transformation matrix
            DNekMatSharedPtr m_inversetransformationmatrix =
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nCoeffs, nCoeffs, 0.0, storage);
            DNekMat &InvR = (*m_inversetransformationmatrix);

            int nVerts = GetNverts();
            int nEdges = GetNedges();
            int nFaces = GetNfaces();

            int nedgemodes = 0;
            int nfacemodes = 0;
            int nedgemodestotal = 0;
            int nfacemodestotal = 0;

            for (eid = 0; eid < nEdges; ++eid)
            {
                nedgemodes = GetEdgeNcoeffs(eid) - 2;
                nedgemodestotal += nedgemodes;
            }

            for (fid = 0; fid < nFaces; ++fid)
            {
                nfacemodes = GetFaceIntNcoeffs(fid);
                nfacemodestotal += nfacemodes;
            }

            Array<OneD, unsigned int>
                edgemodearray(nedgemodestotal);
            Array<OneD, unsigned int>
                facemodearray(nfacemodestotal);

            int offset = 0;

            // Create array of edge modes
            for (eid = 0; eid < nEdges; ++eid)
            {
                Array<OneD, unsigned int> edgearray
                    = GetEdgeInverseBoundaryMap(eid);
                nedgemodes = GetEdgeNcoeffs(eid) - 2;

                // Only copy if there are edge modes
                if (nedgemodes)
                {
                    Vmath::Vcopy(nedgemodes, &edgearray[0], 1,
                                 &edgemodearray[offset], 1);
                }

                offset += nedgemodes;
            }

            offset = 0;

            // Create array of face modes
            for (fid = 0; fid < nFaces; ++fid)
            {
                Array<OneD, unsigned int> facearray
                    = GetFaceInverseBoundaryMap(fid);
                nfacemodes = GetFaceIntNcoeffs(fid);

                // Only copy if there are face modes
                if (nfacemodes)
                {
                    Vmath::Vcopy(nfacemodes, &facearray[0], 1,
                                 &facemodearray[offset], 1);
                }

                offset += nfacemodes;
            }

            // Vertex-edge/face
            for (i = 0; i < nVerts; ++i)
            {
                for (j = 0; j < nedgemodestotal; ++j)
                {
                    InvR.SetValue(
                        GetVertexMap(i), edgemodearray[j],
                        -R(GetVertexMap(i), edgemodearray[j]));
                }
                for (j = 0; j < nfacemodestotal; ++j)
                {
                    InvR.SetValue(
                        GetVertexMap(i), facemodearray[j],
                        -R(GetVertexMap(i), facemodearray[j]));
                    for (n = 0; n < nedgemodestotal; ++n)
                    {
                        MatrixValue = InvR.GetValue(
                            GetVertexMap(i), facemodearray[j])
                            + R(GetVertexMap(i), edgemodearray[n])
                            * R(edgemodearray[n], facemodearray[j]);
                        InvR.SetValue(
                            GetVertexMap(i), facemodearray[j], MatrixValue);
                    }
                }
            }

            // Edge-face contributions
            for (i = 0; i < nedgemodestotal; ++i)
            {
                for (j = 0; j < nfacemodestotal; ++j)
                {
                    InvR.SetValue(
                        edgemodearray[i], facemodearray[j],
                        -R(edgemodearray[i], facemodearray[j]));
                }
            }

            for (i = 0; i < nCoeffs; ++i)
            {
                InvR.SetValue(i, i, 1.0);
            }

            return m_inversetransformationmatrix;
	}

        Array<OneD, unsigned int> Expansion3D::v_GetEdgeInverseBoundaryMap(
            int eid)
        {
            int n, j;
            int nEdgeCoeffs;
            int nBndCoeffs = NumBndryCoeffs();

            Array<OneD, unsigned int> bmap(nBndCoeffs);
            GetBoundaryMap(bmap);

            // Map from full system to statically condensed system (i.e reverse
            // GetBoundaryMap)
            map<int, int> invmap;
            for (j = 0; j < nBndCoeffs; ++j)
            {
                invmap[bmap[j]] = j;
            }

            // Number of interior edge coefficients
            nEdgeCoeffs = GetEdgeNcoeffs(eid) - 2;

            const SpatialDomains::Geometry3DSharedPtr &geom = GetGeom3D();

            Array<OneD, unsigned int> edgemaparray(nEdgeCoeffs);
            StdRegions::Orientation   eOrient  =
                geom->GetEorient(eid);
            Array<OneD, unsigned int> maparray =
                Array<OneD, unsigned int>(nEdgeCoeffs);
            Array<OneD, int> signarray         =
                Array<OneD, int>(nEdgeCoeffs, 1);

            // maparray is the location of the edge within the matrix
            GetEdgeInteriorMap(eid, eOrient, maparray, signarray);

            for (n = 0; n < nEdgeCoeffs; ++n)
            {
                edgemaparray[n] = invmap[maparray[n]];
            }

            return edgemaparray;
        }

        Array<OneD, unsigned int>
        Expansion3D::v_GetFaceInverseBoundaryMap(int fid, StdRegions::Orientation faceOrient)
        {
            int n, j;
            int nFaceCoeffs;

            int nBndCoeffs = NumBndryCoeffs();

            Array<OneD, unsigned int> bmap(nBndCoeffs);
            GetBoundaryMap(bmap);

            // Map from full system to statically condensed system (i.e reverse
            // GetBoundaryMap)
            map<int, int> reversemap;
            for (j = 0; j < bmap.num_elements(); ++j)
            {
                reversemap[bmap[j]] = j;
            }

            // Number of interior face coefficients
            nFaceCoeffs = GetFaceIntNcoeffs(fid);

            Array<OneD, unsigned int> facemaparray(nFaceCoeffs);
            StdRegions::Orientation   fOrient; 
            Array<OneD, unsigned int> maparray  =
                Array<OneD, unsigned int>(nFaceCoeffs);
            Array<OneD, int>          signarray =
                Array<OneD, int>(nFaceCoeffs, 1);

            if(faceOrient == StdRegions::eNoOrientation)
            {
                fOrient = GetForient(fid);
            }
            else
            {
                fOrient = faceOrient; 
            }
            
            // maparray is the location of the face within the matrix
            GetFaceInteriorMap(fid, fOrient, maparray, signarray);

            for (n = 0; n < nFaceCoeffs; ++n)
            {
                facemaparray[n] = reversemap[maparray[n]];
            }

            return facemaparray;
        }

        StdRegions::Orientation Expansion3D::v_GetForient(int face)
        {
            return m_geom->GetForient(face);
        }
    } //end of namespace
} //end of namespace
