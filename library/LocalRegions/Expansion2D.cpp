///////////////////////////////////////////////////////////////////////////////
//
// File Expansion2D.cpp
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
// Description: File for Expansion2D routines
//
///////////////////////////////////////////////////////////////////////////////
//#include <LocalRegions/SegExp.h>
#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        
        void Expansion2D::AddEdgeNormBoundaryInt(const int edge, 
                                                 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                 const Array<OneD, const NekDouble> &Fx,  
                                                 const Array<OneD, const NekDouble> &Fy,  
                                                 Array<OneD, NekDouble> &outarray)
        {
            int nquad_e = EdgeExp->GetNumPoints(0);
            Array<OneD, const NekDouble> normals = EdgeExp->GetPhysNormals();

            ASSERTL1(v_GetCoordim() == 2,"Routine only set up for two-dimensions");
            
            Vmath::Vmul(nquad_e,&(normals[0]),1,&Fx[0],1,
                        &(EdgeExp->UpdatePhys())[0],1);
            Vmath::Vvtvp(nquad_e,&(normals[nquad_e]),1,
                         &Fy[0],1,&(EdgeExp->GetPhys())[0],1,
                         &(EdgeExp->UpdatePhys())[0],1);

            AddEdgeNormBoundaryInt(edge, EdgeExp, EdgeExp->GetPhys(), outarray);
        }

        void Expansion2D::AddEdgeNormBoundaryInt(const int edge, 
						 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
						 const Array<OneD, const NekDouble> &Fn,  
						 Array<OneD, NekDouble> &outarray)
        {
	  int i;
	  Array<OneD,unsigned int>     map;
	  Array<OneD,int>              sign;
	  StdRegions::EdgeOrientation  edgedir = v_GetEorient(edge);
	  
	  v_GetEdgeToElementMap(edge,edgedir,map,sign);
	  int order_e = map.num_elements(); // Order of the element
	  int n_coeffs = (EdgeExp->GetCoeffs()).num_elements(); // Order of the trace
	  
	  if(n_coeffs!=order_e) // Going to orthogonal space
	  {
	      EdgeExp->FwdTrans(Fn,EdgeExp->UpdateCoeffs());
	      // negate backward edge values due to inwards normal definition
	      if(edgedir == StdRegions::eBackwards)
	      {
		  
		  Vmath::Neg(n_coeffs,EdgeExp->UpdateCoeffs(),1);
              }

	      Array<OneD, NekDouble> coeff(n_coeffs,0.0);
	      LibUtilities::BasisType btype = ((LibUtilities::BasisType) 1); //1-->Ortho_A
	      LibUtilities::BasisKey bkey_ortho(btype,EdgeExp->GetBasis(0)->GetNumModes(),EdgeExp->GetBasis(0)->GetPointsKey());
	      LibUtilities::BasisKey bkey(EdgeExp->GetBasis(0)->GetBasisType(),EdgeExp->GetBasis(0)->GetNumModes(),EdgeExp->GetBasis(0)->GetPointsKey());
	      LibUtilities::InterpCoeff1D(bkey,EdgeExp->GetCoeffs(),bkey_ortho,coeff);
	      // Cutting high frequencies
	      for(i = order_e; i < n_coeffs; i++)
		{
		  coeff[i] = 0.0;
		}	
	      LibUtilities::InterpCoeff1D(bkey_ortho,coeff,bkey,EdgeExp->UpdateCoeffs());
	      
	      StdRegions::StdMatrixKey masskey(StdRegions::eMass,StdRegions::eSegment,*EdgeExp);
	      EdgeExp->MassMatrixOp(EdgeExp->UpdateCoeffs(),EdgeExp->UpdateCoeffs(),masskey);
	    }
	  else
	    {
	      EdgeExp->IProductWRTBase(Fn,EdgeExp->UpdateCoeffs());

	      // negate backward edge values due to inwards normal definition
	      if(edgedir == StdRegions::eBackwards)
		{
		  Vmath::Neg(order_e,EdgeExp->UpdateCoeffs(),1);
		}
	    }
	  
	  // add data to outarray if forward edge normal is outwards
	  for(i = 0; i < order_e; ++i)
	    {
	      outarray[map[i]] += sign[i]*EdgeExp->GetCoeff(i);
	    }
        }

        void Expansion2D::AddEdgeNormBoundaryBiInt(const int edge, 
                                                 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                 const Array<OneD, const NekDouble> &Fwd,  
                                                 const Array<OneD, const NekDouble> &Bwd,  
                                                 Array<OneD, NekDouble> &outarray)
        {
            int i;
            int order_e = EdgeExp->GetNcoeffs();
            Array<OneD,unsigned int>     map;
            Array<OneD,int>              sign;
            StdRegions::EdgeOrientation  edgedir = v_GetEorient(edge);
            
            //    ASSERTL1(v_GetCoordim() == 2,"Routine only set up for two-dimensions");

            v_GetEdgeToElementMap(edge,edgedir,map,sign);

            if(edgedir == StdRegions::eForwards)
            {
                EdgeExp->IProductWRTBase(Fwd,EdgeExp->UpdateCoeffs());
            }

            else if(edgedir == StdRegions::eBackwards)
            {
                EdgeExp->IProductWRTBase(Bwd,EdgeExp->UpdateCoeffs());
            }

            // add data to outarray if forward edge normal is outwards
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*EdgeExp->GetCoeff(i);
            }
        }

        void Expansion2D::SetTraceToGeomOrientation(Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,  Array<OneD, NekDouble> &inout)
        {
            int i,cnt = 0;
            int nedges = v_GetNedges();
            Array<OneD, NekDouble> e_tmp;
            
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i]->SetCoeffsToOrientation(v_GetEorient(i),
                                                   e_tmp = inout + cnt, 
                                                   e_tmp = inout + cnt);
                cnt += v_GetEdgeNcoeffs(i);
            }
        }

        void Expansion2D::AddNormTraceInt(const int dir,
                                          Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                          Array<OneD,NekDouble> &outarray) 
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = v_GetNedges();

            Array<OneD,NekDouble> normals;

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);
                normals = EdgeExp[e]->GetPhysNormals();
                
                for(i = 0; i < order_e; ++i)
                {
                    EdgeExp[e]->SetCoeff(i,inarray[i+cnt]);
                }
                cnt += order_e;
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                     EdgeExp[e]->UpdatePhys());
                
                Vmath::Vmul(nquad_e,&(normals[dir*nquad_e]),1,
                            &(EdgeExp[e]->GetPhys())[0],1,
                            &(EdgeExp[e]->UpdatePhys())[0],1);

                // negate backwards normal
                if(v_GetEorient(e) == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                }
                
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }
        }

        void Expansion2D::AddNormTraceInt(const int dir,
                                          Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                          const Array<OneD, NekDouble> &directional,
                                          Array<OneD,NekDouble> &outarray)
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = v_GetNedges();
            int coordim = v_GetCoordim();

            Array<OneD,NekDouble> normals;

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);
                normals = EdgeExp[e]->GetPhysNormals();
                
                for(i = 0; i < order_e; ++i)
                {
                    EdgeExp[e]->SetCoeff(i,inarray[i+cnt]);
                }
                cnt += order_e;
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                     EdgeExp[e]->UpdatePhys());
                
                if(directional.num_elements()>0)
		{
                    Array<OneD, NekDouble> normalindir(nquad_e);
                    Getnormalindir(e,EdgeExp[e],normals,directional,normalindir);
                    Vmath::Vmul(nquad_e,&normalindir[0],1,
                                &(EdgeExp[e]->GetPhys())[0],1,
                                &(EdgeExp[e]->UpdatePhys())[0],1);
                }                

                else
                {
                    Vmath::Vmul(nquad_e,&(normals[dir*nquad_e]),1,
                                &(EdgeExp[e]->GetPhys())[0],1,
                                &(EdgeExp[e]->UpdatePhys())[0],1);
                }

                // negate backwards normal
                if(v_GetEorient(e) == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                }
                
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }
        }

        void Expansion2D::AddNormTraceInt(const int dir,
                                          Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                          Array<OneD,NekDouble> &outarray) 
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = v_GetNedges();

            Array<OneD,NekDouble> normals;

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);
                normals = EdgeExp[e]->GetPhysNormals();
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                     EdgeExp[e]->UpdatePhys());
                
                Vmath::Vmul(nquad_e,&(normals[dir*nquad_e]),1,
                            &(EdgeExp[e]->GetPhys())[0],1,
                            &(EdgeExp[e]->UpdatePhys())[0],1);

                // negate backwards normal
                if(v_GetEorient(e) == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                }
                
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }
        }


        void Expansion2D:: AddEdgeBoundaryInt(const int edge, 
                                              const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                              Array <OneD,NekDouble > &outarray)
        {
            
            int i;
            int order_e = EdgeExp->GetNcoeffs();                    
            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;
            Array<OneD, NekDouble> coeff(order_e);

            v_GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);

            EdgeExp->IProductWRTBase(EdgeExp->GetPhys(),coeff);

            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*coeff[i];
            }
        }

        
        // This method assumes that data in EdgeExp is orientated 
        // according to elemental counter clockwise format
        // AddHDGHelmholtzTraceTerms with directions
        void Expansion2D::AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                    const Array<OneD, const NekDouble> &inarray, 
                                                    Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,  
						    const Array<OneD, Array<OneD, NekDouble> > &dirForcing,
                                                    Array<OneD,NekDouble> &outarray,
                                                    const int matrixid)
        {
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");
            
            
            int e,cnt;
            int order_e;
            int nedges = v_GetNedges();
            Array<OneD, const NekDouble> tmp;
            
            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                
                order_e = EdgeExp[e]->GetNcoeffs();                    
                
                Vmath::Vcopy(order_e,tmp =inarray+cnt,1,EdgeExp[e]->UpdateCoeffs(),1);
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),EdgeExp[e]->UpdatePhys());
                
                AddHDGHelmholtzEdgeTerms(tau,e,EdgeExp,dirForcing,outarray,matrixid);
                
                cnt += order_e;
            }
        }
        
        //  evaluate additional terms in HDG edges. Not that this assumes that
        // edges are unpacked into local cartesian order. 
        void Expansion2D::AddHDGHelmholtzEdgeTerms(const NekDouble tau, 
                                                   const int edge,
                                                   Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp,
						   const Array<OneD, Array<OneD, NekDouble> > &dirForcing, 
                                                   Array <OneD,NekDouble > &outarray,
                                                   const int matrixid)
        {            
            int i,j,n;
            int nquad_e = EdgeExp[edge]->GetNumPoints(0); 
            int order_e = EdgeExp[edge]->GetNcoeffs();            
            int coordim = v_GetCoordim();
            int ncoeffs  = v_GetNcoeffs();

            Array<OneD, NekDouble> inval   (nquad_e);
            Array<OneD, NekDouble> outcoeff(order_e);
            Array<OneD, NekDouble> tmpcoeff(ncoeffs);
            Array<OneD, const NekDouble> normals = EdgeExp[edge]->GetPhysNormals();
            Array<OneD,unsigned int> emap;
            Array<OneD,int> sign;

            DNekScalMat  &invMass = *v_GetLocMatrix(StdRegions::eInvMass);
            
            StdRegions::EdgeOrientation edgedir = v_GetEorient(edge);

            DNekVec                Coeffs  (ncoeffs,outarray,eWrapper);
            DNekVec                Tmpcoeff(ncoeffs,tmpcoeff,eWrapper);
            
            v_GetEdgeToElementMap(edge,edgedir,emap,sign);

            //================================================================
            // Add F = \tau <phi_i,in_phys>
            // Fill edge and take inner product
            EdgeExp[edge]->IProductWRTBase(EdgeExp[edge]->GetPhys(),
                                           EdgeExp[edge]->UpdateCoeffs());
            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[emap[i]] += sign[i]*tau*EdgeExp[edge]->GetCoeff(i);
            }
            //================================================================

            //===============================================================
            // Add -\sum_i D_i^T M^{-1} G_i + E_i M^{-1} G_i = 
            //                         \sum_i D_i M^{-1} G_i term
            StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                   StdRegions::eWeakDeriv1,
                                                   StdRegions::eWeakDeriv2};

	    int nvarcoeffs = dirForcing.num_elements();

	    // Two independent direction
            for(n = 0; n < 2; ++n)
            {
                //G;
                if(nvarcoeffs>0)
		{
                    Array<OneD, NekDouble> normalindir(nquad_e);
                    Getnormalindir(edge,EdgeExp[edge],normals,dirForcing[n],normalindir);
                    Vmath::Vmul(nquad_e,&normalindir[0],1,&(EdgeExp[edge]->GetPhys())[0],1,&inval[0],1);
                }                
                
                else
                {
                    Vmath::Vmul(nquad_e,&normals[n*nquad_e],1,&(EdgeExp[edge]->GetPhys())[0],1,&inval[0],1);
                }
              
                // negate for backwards facing edge
                if(edgedir == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,inval,1);
                }
                
                EdgeExp[edge]->IProductWRTBase(inval,outcoeff);
                
                // M^{-1} G
                for(i = 0; i < ncoeffs; ++i)
                {
                    tmpcoeff[i] = 0;
                    for(j = 0; j < order_e; ++j)
                    {
                        tmpcoeff[i] += invMass(i,emap[j])*sign[j]*outcoeff[j];
                    }
                }

                if(nvarcoeffs>0)
                {
                    DNekScalMat &Dmat = *v_GetLocMatrix(StdRegions::eWeakDirectionalDeriv,dirForcing[n],matrixid+n*10000);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;                 
                }

                else
                {
                    DNekScalMat &Dmat = *v_GetLocMatrix(DerivType[n]);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;       
                }
            }
        }

      void Expansion2D::Getnormalindir(const int edge,
				       StdRegions::StdExpansion1DSharedPtr &EdgeExp_e,
				       const Array<OneD, const NekDouble> &normals, 
				       const Array<OneD, const NekDouble> &directional,
				       Array<OneD, NekDouble> &outarray)
      {
            int nquad_e = EdgeExp_e->GetNumPoints(0); 
            int coordim = v_GetCoordim();
            int nq = (directional.num_elements())/coordim;
            
            Array<OneD, NekDouble> dirForcing_e(nquad_e);
            Array<OneD, NekDouble> dirtemp(nq);

	    // Initialize of outarray as zero vectors
	    Vmath::Zero(nquad_e,outarray,1);
            
            for(int k=0; k < coordim; ++k)
            {
                // direction k for nth tangential basis is copied into dirtemp
                Vmath::Vcopy(nq, &directional[k*nq], 1, &dirtemp[0], 1);
                
                // get edge values of dirtemp
                v_GetEdgePhysVals(edge, EdgeExp_e, dirtemp, dirForcing_e);
                
                // new_normal = nx*dirForcingx + ny*dirForcingy + nz*dirForcingz
                Vmath::Vvtvp(nquad_e,&dirForcing_e[0],1,&normals[k*nquad_e],1,&outarray[0],1,&outarray[0],1);
            }
        }

        DNekMatSharedPtr Expansion2D::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
                {
                    
                    ASSERTL1(v_IsBoundaryInteriorExpansion(),
                             "HybridDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");
                    
                    int       i,j,k;
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    int       ncoeffs   = v_GetNcoeffs();
                    int       nedges    = v_GetNedges();
                    int       matrixid   = mkey.GetMatrixID();

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;
                    StdRegions::EdgeOrientation edgedir;
                    StdRegions::StdExpansion1DSharedPtr EdgeExp;

                    int order_e, coordim = v_GetCoordim();
                    DNekScalMat  &invMass = *v_GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};
                    DNekMat LocMat(ncoeffs,ncoeffs); 

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;

                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    int nvarcoeffs = mkey.GetNvariableCoefficients();

                    for(i=0;  i < 2; ++i)
                    {
                        if(nvarcoeffs>0)
                        {
                            DNekScalMat &Dmat = *v_GetLocMatrix(StdRegions::eWeakDirectionalDeriv,
                                                                mkey.GetVariableCoefficient(i),matrixid+i*10000);
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }

                        else
                        {
                            DNekScalMat &Dmat = *v_GetLocMatrix(DerivType[i]);
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }

                    }

                    // Add Mass Matrix Contribution
                    DNekScalMat  &Mass = *v_GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;                    

                    // Add tau*F_e using elemental mass matrices
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp = v_GetEdgeExp(i);
                        DNekScalMat &eMass = *EdgeExp->GetLocMatrix(StdRegions::eMass);
                        order_e = EdgeExp->GetNcoeffs();  
                        v_GetEdgeToElementMap(i,edgedir,emap,sign);

                        for(j = 0; j < order_e; ++j)
                        {
                            for(k = 0; k < order_e; ++k)
                            {
                                Mat(emap[j],emap[k]) = Mat(emap[j],emap[k]) + tau*sign[j]*sign[k]*eMass(j,k);
                            }
                        }
                    }
                }
                break;
            case StdRegions::eHybridDGLamToU:
                {
                    int i,j,k;
                    int nbndry = v_NumDGBndryCoeffs();
                    int ncoeffs = v_GetNcoeffs();
                    int nedges  = v_GetNedges();
                    int matrixid = mkey.GetMatrixID();
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    
                    Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Umat = *returnval;
                    
                    int nvarcoeffs = mkey.GetNvariableCoefficients();
                    Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);
                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }

                    // Helmholtz matrix
                    //  DNekScalMat  &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, lambdaval, tau);
		    DNekScalMat  &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, varcoeffs, matrixid, lambdaval, tau);

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;
                    
#if 0 
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = v_GetEdgeExp(i);
                      
                        DNekScalMat &MassMat = *(EdgeExp[i]->v_GetLocMatrix(StdRegions::eMass));
                        Array<OneD, const NekDouble> normals = EdgeExp[edge]->GetPhysNormals();
                        MatrixKey    mass_n0(StdRegions::eMass,*this,normals[0]);
                        DNekScalMat &MassMatN0 = *(EdgeExp[i]->v_GetLocMatrix(mass_n0));
                        MatrixKey    mass_n1(StdRegions::eMass,*this,normals[1]);
                        DNekScalMat &MassMatN0 = *(EdgeExp[i]->v_GetLocMatrix(mass_n1));

                        StdRegions::EdgeOrientation edgedir = v_GetEorient(edge);
                        v_GetEdgeToElementMap(edge,edgedir,emap,sign);
                        
                        // sum up appropriate terms in matrix 
                        

                        // D M^{-1}
                        
                    }

#else
                    for(i = 0; i < nedges; ++i)
                    {
                      EdgeExp[i] = v_GetEdgeExp(i);
                    }

                    // for each degree of freedom of the lambda space
                    // calculate Umat entry 
                    // Generate Lambda to U_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(ncoeffs,&f[0],1);
                        lambda[j] = 1.0;
                        
                        SetTraceToGeomOrientation(EdgeExp,lambda);
                        
                        //   AddHDGHelmholtzTraceTerms(tau,lambda,EdgeExp,f);
                        AddHDGHelmholtzTraceTerms(tau, lambda, EdgeExp, varcoeffs, f, matrixid);
                        
                        Ulam = invHmat*F; // generate Ulam from lambda
                        
                        // fill column of matrix
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Umat(k,j) = Ulam[k]; 
                        }
                    }
#endif

                }
                break;
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
                {
                    int i,j,k,dir;
                    int nbndry = v_NumDGBndryCoeffs();
                    int nquad  = v_GetNumPoints(0);
                    int ncoeffs = v_GetNcoeffs();
                    int nedges  = v_GetNedges();
                    int matrixid = mkey.GetMatrixID();

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    
                    Array<OneD,NekDouble> ulam(ncoeffs);
                    DNekVec Ulam(ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(ncoeffs);
                    DNekVec F(ncoeffs,f,eWrapper);
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    

                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    int nvarcoeffs = mkey.GetNvariableCoefficients();
                    Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);

                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }

                    // Helmholtz matrix
                    // DNekScalMat &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, lambdaval,tau);
                    
                    // Lambda to U matrix
                    // DNekScalMat &lamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval, tau);
                    
                    DNekScalMat  &invHmat = *v_GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, varcoeffs, matrixid, lambdaval, tau);
                    
                    // Lambda to U matrix
                    DNekScalMat &lamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, varcoeffs, matrixid, lambdaval, tau);

                    // Inverse mass matrix 
                    DNekScalMat &invMass = *v_GetLocMatrix(StdRegions::eInvMass);
                    
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = v_GetEdgeExp(i);
                    }

                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        if(nvarcoeffs>0)
                        {
                            Dmat = v_GetLocMatrix(StdRegions::eWeakDirectionalDeriv,varcoeffs[dir],matrixid+dir*10000);
                        }

                        else
                        {
                            Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv0); 
                        }
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        if(nvarcoeffs>0)
                        {
                            Dmat = v_GetLocMatrix(StdRegions::eWeakDirectionalDeriv,varcoeffs[dir],matrixid+dir*10000);
                        }

                        else
                        {
                            Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv1); 
                        }
                        break;
                    case StdRegions::eHybridDGLamToQ2:
                        dir = 2;
                        Dmat = v_GetLocMatrix(StdRegions::eWeakDeriv2); 
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
                        
                        SetTraceToGeomOrientation(EdgeExp,lambda);
                        
                        // + \tilde{G} \lambda
                        if(nvarcoeffs>0)
                        {
                            AddNormTraceInt(dir,lambda,EdgeExp,varcoeffs[dir],f); 
                        }
                        else
                        {
                            AddNormTraceInt(dir,lambda,EdgeExp,f); 
                        }
                        
                        // multiply by inverse mass matrix
                        Ulam = invMass*F; 
                        
                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*ncoeffs,1);
                    }
                }
                break;            
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int i,j,e,cnt;
                    int order_e, nquad_e;
                    int nbndry  = v_NumDGBndryCoeffs();
                    int coordim = v_GetCoordim();
                    int nedges  = v_GetNedges();
                    int matrixid = mkey.GetMatrixID();
                    
                    Array<OneD,NekDouble>       work;
                    Array<OneD,const NekDouble> normals; 
                    Array<OneD, NekDouble> normalindir;
                    Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    Array<OneD, NekDouble> lam(nbndry); 
                    
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    Array<OneD,unsigned int>    emap;
                    Array<OneD, int>            sign;
                    StdRegions::EdgeOrientation edgedir;
                    
                    int nvarcoeffs = mkey.GetNvariableCoefficients();
                    Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);
                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }

                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    // DNekScalMat &LamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval,tau); 
                    DNekScalMat &LamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, varcoeffs, matrixid, lambdaval, tau);                

                    // Matrix to map Lambda to Q0
                    // DNekScalMat &LamToQ0 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ0, lambdaval,tau);
                    DNekScalMat &LamToQ0 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ0, varcoeffs, matrixid, lambdaval,tau);
 
                    // Matrix to map Lambda to Q1
                    // DNekScalMat &LamToQ1 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ1, lambdaval,tau);
                    DNekScalMat &LamToQ1 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ1, varcoeffs, matrixid, lambdaval,tau);


                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = v_GetEdgeExp(i);
                    }

                    // Set up matrix derived from <mu, Q_lam.n - \tau (U_lam - Lam) > 
                    for(i = 0; i < nbndry; ++i)
                    {
                        cnt = 0;
                        
                        Vmath::Zero(nbndry,lam,1);
                        lam[i] = 1.0;
                        SetTraceToGeomOrientation(EdgeExp,lam);

                        for(e = 0; e < nedges; ++e)
                        {
                            order_e = EdgeExp[e]->GetNcoeffs();  
                            nquad_e = EdgeExp[e]->GetNumPoints(0);    
                            normals = EdgeExp[e]->GetPhysNormals();
                            edgedir = v_GetEorient(e);
                            
                            work = Array<OneD,NekDouble>(nquad_e);
                            
                            v_GetEdgeToElementMap(e,edgedir,emap,sign);

                            // Q0 * n0
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,sign[j]*LamToQ0(emap[j],i));
                            }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());
          

                            if(nvarcoeffs>0)
                            {
                                normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(0),normalindir);
                                Vmath::Vmul(nquad_e,normalindir,1,EdgeExp[e]->GetPhys(),1,work,1);
                            }

                            else
                            {
                                Vmath::Vmul(nquad_e,normals,1,EdgeExp[e]->GetPhys(),1,work,1);
                            }
                            
                            if(edgedir == StdRegions::eBackwards)
                            {
                                Vmath::Neg(nquad_e,work,1);
                            }

                            // Q1 * n1
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,sign[j]*LamToQ1(emap[j],i));
                            }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());

                            if(edgedir == StdRegions::eForwards)
                            {
                                if(nvarcoeffs>0)
                                {
                                    // normalindir = normals /cdot dirForcing[1]
                                    normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                    Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(1),normalindir);
                                    Vmath::Vvtvp(nquad_e,&normalindir[0],1,&EdgeExp[e]->GetPhys()[0],1,&work[0],1,&work[0],1);
                                }

                                else
                                {
                                    Vmath::Vvtvp(nquad_e,&normals[nquad_e],1,
                                                 &EdgeExp[e]->GetPhys()[0],1,
                                                 &work[0],1,&work[0],1);
                                }
                            }
                            else // subtrace values for negative normal
                            {
                                if(nvarcoeffs>0)
                                {
                                    // normalindir = normals /cdot dirForcing[1]
                                    normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                    Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(1),normalindir);
                                    Vmath::Vvtvm(nquad_e,&normalindir[0],1,&EdgeExp[e]->GetPhys()[0],1,&work[0],1,&work[0],1);
                                    Vmath::Neg(nquad_e,work,1);
                                }
                                
                                else
                                {
                                    Vmath::Vvtvm(nquad_e,&normals[nquad_e],1,
                                                 &EdgeExp[e]->GetPhys()[0],1,
                                                 &work[0],1,&work[0],1);
                                    Vmath::Neg(nquad_e,work,1);
                                }
                            }

                            // - tau (ulam - lam)
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,sign[j]*LamToU(emap[j],i) - lam[cnt+j]);
                            }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());
                            
                            Vmath::Svtvp(nquad_e,-tau,EdgeExp[e]->GetPhys(),1,
                                         work,1,work,1);

                            EdgeExp[e]->IProductWRTBase(work,EdgeExp[e]->UpdateCoeffs());
                            
                            EdgeExp[e]->SetCoeffsToOrientation(edgedir);
                            
                            for(j = 0; j < order_e; ++j)
                            {
                                BndMat(cnt+j,i) = EdgeExp[e]->GetCoeff(j);
                            }
                            
                            cnt += order_e;
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

      //Evaluate Coefficients of weak deriviative in the direction dir
      //given the input coefficicents incoeffs and the imposed
      //boundary values in EdgeExp (which will have its phys space updated);
      void Expansion2D::DGDeriv(int dir, 
                                const Array<OneD, const NekDouble>&incoeffs,
                                Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                Array<OneD, NekDouble> &out_d)
      {
        StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                               StdRegions::eWeakDeriv1,
                                               StdRegions::eWeakDeriv2};
        
          int ncoeffs = v_GetNcoeffs();
          int nedges  = v_GetNedges();
          
          DNekScalMat &InvMass = *v_GetLocMatrix(StdRegions::eInvMass);
          DNekScalMat &Dmat    = *v_GetLocMatrix(DerivType[dir]);
          
          Array<OneD, NekDouble> coeffs = incoeffs;
          DNekVec     Coeffs  (ncoeffs,coeffs, eWrapper);
          
          Coeffs = Transpose(Dmat)*Coeffs;
          Vmath::Neg(ncoeffs, coeffs,1);

          // Add the boundary integral including the relevant part of
          // the normal
          AddNormTraceInt(dir,EdgeExp,coeffs);
        
          DNekVec Out_d (ncoeffs,out_d,eWrapper);

          Out_d  = InvMass*Coeffs;
      }
    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion2D.cpp,v $
 *    Revision 1.15  2009/11/09 18:12:55  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.14  2009/11/09 15:43:51  sehunchun
 *    HDG2DManifold Solver with Variable coefficients
 *
 *    Revision 1.13  2009/11/06 21:43:56  sherwin
 *    DGDeriv function
 *
 *    Revision 1.12  2009/10/22 18:05:33  cbiotto
 *    Update for variable order expansion
 *
 *    Revision 1.11  2009/09/06 22:24:00  sherwin
 *    Updates for Navier-Stokes solver
 *
 *    Revision 1.10  2009/07/07 16:31:47  sehunchun
 *    Adding AddEdgeBoundaryBiInt to line integrate depending on Fwd and Bwd
 *
 *    Revision 1.9  2009/04/28 09:58:17  sherwin
 *    Updates to make HDG implementation consistent in 1D as to the 2D implementation
 *
 *    Revision 1.8  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.7  2009/04/20 16:12:28  sherwin
 *    Updates related to output format and optimising DG solver
 *
 *    Revision 1.6  2009/04/02 13:04:36  sherwin
 *    Modified Hybrid Matrix call to use matrix D M^{-1}D' formulation and removed operations based version
 *
 *    Revision 1.5  2008/10/04 19:34:09  sherwin
 *    Added an upwind method which takes the normal flux rather than than individual components
 *
 *    Revision 1.4  2008/08/27 16:35:13  pvos
 *    Small efficiency update
 *
 *    Revision 1.3  2008/08/20 09:16:39  sherwin
 *    Modified generation of HDG matrices so that they use Expansion1D, Expansion2D GenMatrix method rather than Expansion method. Have also removed methods which were generating edge expansions locally as this was too expensive
 *
 *    Revision 1.2  2008/08/18 08:30:36  sherwin
 *    Updates for HDG 1D work
 *
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
