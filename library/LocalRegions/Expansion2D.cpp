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
#include <SpatialDomains/Geometry.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        Expansion2D::Expansion2D() : 
            StdExpansion2D(), Expansion(), StdExpansion()
        {
        }
        
        void Expansion2D::v_AddEdgeNormBoundaryInt(const int edge,
                                                 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                 const Array<OneD, const NekDouble> &Fx,  
                                                 const Array<OneD, const NekDouble> &Fy,  
                                                 Array<OneD, NekDouble> &outarray)
        {
            ASSERTL1(GetCoordim() == 2,"Routine only set up for two-dimensions");

            const Array<OneD, const Array<OneD, NekDouble> > normals
                                    = GetEdgeNormal(edge);

            // We allow the case of mixed polynomial order by supporting only
            // those modes on the edge common to both adjoining elements. This
            // is enforced here by taking the minimum size and padding with
            // zeros.
            int nquad_e = min(EdgeExp->GetNumPoints(0), int(normals[0].num_elements()));
            int coordim = GetCoordim();

            Vmath::Zero(EdgeExp->GetNumPoints(0),EdgeExp->UpdatePhys(),1);
            Vmath::Vmul(nquad_e,normals[0],1,Fx,1,
                        EdgeExp->UpdatePhys(),1);
            Vmath::Vvtvp(nquad_e,normals[1],1,
                         Fy,1,EdgeExp->GetPhys(),1,
                         EdgeExp->UpdatePhys(),1);

            AddEdgeNormBoundaryInt(edge, EdgeExp, EdgeExp->GetPhys(), outarray);
        }

        void Expansion2D::v_AddEdgeNormBoundaryInt(const int edge,
						 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
						 const Array<OneD, const NekDouble> &Fn,  
						 Array<OneD, NekDouble> &outarray)
        {
	  int i;
	  Array<OneD,unsigned int>     map;
	  Array<OneD,int>              sign;
	  StdRegions::EdgeOrientation  edgedir = GetEorient(edge);
	  
	  GetEdgeToElementMap(edge,edgedir,map,sign);
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

        void Expansion2D::v_AddEdgeNormBoundaryBiInt(const int edge,
                                                 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                 const Array<OneD, const NekDouble> &Fwd,  
                                                 const Array<OneD, const NekDouble> &Bwd,  
                                                 Array<OneD, NekDouble> &outarray)
        {
            int i;
            int order_e = EdgeExp->GetNcoeffs();
            Array<OneD,unsigned int>     map;
            Array<OneD,int>              sign;
            StdRegions::EdgeOrientation  edgedir = GetEorient(edge);
            
            //    ASSERTL1(v_GetCoordim() == 2,"Routine only set up for two-dimensions");

            GetEdgeToElementMap(edge,edgedir,map,sign);

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
            int nedges = GetNedges();
            Array<OneD, NekDouble> e_tmp;
            
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i]->SetCoeffsToOrientation(GetEorient(i),
                                                   e_tmp = inout + cnt, 
                                                   e_tmp = inout + cnt);
                cnt += GetEdgeNcoeffs(i);
            }
        }

        void Expansion2D::v_AddNormTraceInt(const int dir,
                                          Array<OneD, const NekDouble> &inarray,
                                          Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                          Array<OneD,NekDouble> &outarray,
                                          const Array<OneD, NekDouble> &directional)
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = GetNedges();
            int coordim = GetCoordim();

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);

                const Array<OneD, const Array<OneD, NekDouble> > normals = GetEdgeNormal(e);
                
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
                    Vmath::Vmul(nquad_e,normalindir,1,
                                EdgeExp[e]->GetPhys(),1,
                                EdgeExp[e]->UpdatePhys(),1);
                }                

                else
                {
                    Vmath::Vmul(nquad_e,normals[dir],1,
                                EdgeExp[e]->GetPhys(),1,
                                EdgeExp[e]->UpdatePhys(),1);
                }

                // negate backwards normal
                if(GetEorient(e) == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,EdgeExp[e]->UpdatePhys(),1);
                }
                
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }
        }

        void Expansion2D::v_AddNormTraceInt(const int dir,
                                          Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                          Array<OneD,NekDouble> &outarray) 
        {
            int e,cnt;
            int order_e,nquad_e;
            int nedges = GetNedges();

            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);

                const Array<OneD, const Array<OneD, NekDouble> > normals = GetEdgeNormal(e);
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                     EdgeExp[e]->UpdatePhys());
                
                Vmath::Vmul(nquad_e,normals[dir],1,
                            EdgeExp[e]->GetPhys(),1,
                            EdgeExp[e]->UpdatePhys(),1);

                // negate backwards normal
                if(GetEorient(e) == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,EdgeExp[e]->UpdatePhys(),1);
                }
                
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }
        }


        void Expansion2D::v_AddEdgeBoundaryInt(const int edge,
                                              const StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                              Array <OneD,NekDouble > &outarray)
        {
            
            int i;
            int order_e = EdgeExp->GetNcoeffs();                    
            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;
            Array<OneD, NekDouble> coeff(order_e);

            GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);

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
        void Expansion2D::v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                    const Array<OneD, const NekDouble> &inarray, 
                                                    Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,  
						    const Array<OneD, Array<OneD, const NekDouble> > &dirForcing,
                                                    Array<OneD,NekDouble> &outarray)
        {
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");
            
            
            int e,cnt;
            int order_e;
            int nedges = GetNedges();
            Array<OneD, const NekDouble> tmp;
            
            cnt = 0;
            for(e = 0; e < nedges; ++e)
            {
                
                order_e = EdgeExp[e]->GetNcoeffs();                    
                
                Vmath::Vcopy(order_e,tmp =inarray+cnt,1,EdgeExp[e]->UpdateCoeffs(),1);
                
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),EdgeExp[e]->UpdatePhys());
                
                AddHDGHelmholtzEdgeTerms(tau,e,EdgeExp,dirForcing,outarray);
                
                cnt += order_e;
            }
        }
        
        //  evaluate additional terms in HDG edges. Not that this assumes that
        // edges are unpacked into local cartesian order. 
        void Expansion2D::v_AddHDGHelmholtzEdgeTerms(const NekDouble tau,
                                                   const int edge,
                                                   Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp,
						   const Array<OneD, Array<OneD, const  NekDouble> > &dirForcing, 
                                                   Array <OneD,NekDouble > &outarray)
        {            
            int i,j,n;
            int nquad_e = EdgeExp[edge]->GetNumPoints(0); 
            int order_e = EdgeExp[edge]->GetNcoeffs();            
            int coordim = GetCoordim();
            int ncoeffs  = GetNcoeffs();

            Array<OneD, NekDouble> inval   (nquad_e);
            Array<OneD, NekDouble> outcoeff(order_e);
            Array<OneD, NekDouble> tmpcoeff(ncoeffs);

            const Array<OneD, const Array<OneD, NekDouble> > normals
                                = GetEdgeNormal(edge);

            Array<OneD,unsigned int> emap;
            Array<OneD,int> sign;

            DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
            
            StdRegions::EdgeOrientation edgedir = GetEorient(edge);

            DNekVec                Coeffs  (ncoeffs,outarray,eWrapper);
            DNekVec                Tmpcoeff(ncoeffs,tmpcoeff,eWrapper);
            
            GetEdgeToElementMap(edge,edgedir,emap,sign);

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
            for(n = 0; n < coordim; ++n)
            {

                //G;
                  if(nvarcoeffs>0)
                  {
                      Array<OneD, NekDouble> normalindir(nquad_e);
                      Getnormalindir(edge,EdgeExp[edge],normals,dirForcing[n],normalindir);
                      Vmath::Vmul(nquad_e,normalindir,1,EdgeExp[edge]->GetPhys(),1,inval,1);
                   }                
                  
                  else
                  {
                      Vmath::Vmul(nquad_e,normals[n],1,EdgeExp[edge]->GetPhys(),1,inval,1);
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
                    DNekScalMat &Dmat = *GetLocMatrix(StdRegions::eWeakDirectionalDeriv,dirForcing[n]);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;                 
                }

                else
                {
                    DNekScalMat &Dmat = *GetLocMatrix(DerivType[n]);
                    Coeffs = Coeffs  + Dmat*Tmpcoeff;       
                }
            }
        }

      void Expansion2D::Getnormalindir(const int edge,
				       StdRegions::StdExpansion1DSharedPtr &EdgeExp_e,
				       const Array<OneD, const Array<OneD, NekDouble> > &normals, 
				       const Array<OneD, const NekDouble> &directional,
				       Array<OneD, NekDouble> &outarray)
      {
            int nquad_e = EdgeExp_e->GetNumPoints(0); 
            int coordim = GetCoordim();
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
                GetEdgePhysVals(edge, EdgeExp_e, dirtemp, dirForcing_e);
                
                // new_normal = nx*dirForcingx + ny*dirForcingy + nz*dirForcingz
                Vmath::Vvtvp(nquad_e,&dirForcing_e[0],1,&normals[k][0],1,&outarray[0],1,&outarray[0],1);
            }
        }

        DNekMatSharedPtr Expansion2D::v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
                {
                    ASSERTL1(IsBoundaryInteriorExpansion(),
                             "HybridDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");
                    
                    int       i,j,k;
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    int       ncoeffs   = GetNcoeffs();
                    int       nedges    = GetNedges();

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;
                    StdRegions::EdgeOrientation edgedir = StdRegions::eForwards;
                    StdRegions::StdExpansion1DSharedPtr EdgeExp;

                    int order_e, coordim = GetCoordim();
                    DNekScalMat  &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                                           StdRegions::eWeakDeriv1,
                                                           StdRegions::eWeakDeriv2};
                    DNekMat LocMat(ncoeffs,ncoeffs); 

                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(ncoeffs,ncoeffs);
                    DNekMat &Mat = *returnval;
                    Vmath::Zero(ncoeffs*ncoeffs,Mat.GetPtr(),1);

                    int nvarcoeffs = mkey.GetNvariableCoefficients();

                    for(i=0;  i < coordim; ++i)
                    {
                        if(nvarcoeffs>0)
                        {
                            DNekScalMat &Dmat = *GetLocMatrix(StdRegions::eWeakDirectionalDeriv,
                                                                mkey.GetVariableCoefficient(i));
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }

                        else
                        {
                            DNekScalMat &Dmat = *GetLocMatrix(DerivType[i]);
                            Mat = Mat + Dmat*invMass*Transpose(Dmat);
                        }

                    }

                    // Add Mass Matrix Contribution
                    DNekScalMat  &Mass = *GetLocMatrix(StdRegions::eMass);
                    Mat = Mat + lambdaval*Mass;                    

                    // Add tau*F_e using elemental mass matrices
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp = GetEdgeExp(i);
                        DNekScalMat &eMass = *EdgeExp->GetLocMatrix(StdRegions::eMass);
                        order_e = EdgeExp->GetNcoeffs();  
                        GetEdgeToElementMap(i,edgedir,emap,sign);

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
                    int nbndry = NumDGBndryCoeffs();
                    int ncoeffs = GetNcoeffs();
                    int nedges  = GetNedges();
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

                    Array<OneD, Array<OneD, const NekDouble> > varcoeffs(nvarcoeffs);
                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }

                    // Helmholtz matrix
		    DNekScalMat  &invHmat = *GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, varcoeffs, lambdaval, tau);

                    Array<OneD,unsigned int> emap;
                    Array<OneD,int> sign;
                    
                    for(i = 0; i < nedges; ++i)
                    {
                      EdgeExp[i] = GetEdgeExp(i);
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
                        
                        AddHDGHelmholtzTraceTerms(tau, lambda, EdgeExp, varcoeffs, f);
                        
                        Ulam = invHmat*F; // generate Ulam from lambda
                        
                        // fill column of matrix
                        for(k = 0; k < ncoeffs; ++k)
                        {
                            Umat(k,j) = Ulam[k]; 
                        }
                    }
                }
                break;
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
                {
                    int i,j,k,dir;
                    int nbndry = NumDGBndryCoeffs();
                    int nquad  = GetNumPoints(0);
                    int ncoeffs = GetNcoeffs();
                    int nedges  = GetNedges();

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
                    Array<OneD, Array<OneD,const NekDouble> > varcoeffs(nvarcoeffs);

                    if(nvarcoeffs>0)
                    {
                        for(int j=0; j<nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j);
                        }
                    }
                    
                    DNekScalMat  &invHmat = *GetLocMatrix(StdRegions::eInvHybridDGHelmholtz, varcoeffs, lambdaval, tau);
                    
                    // Lambda to U matrix
                    DNekScalMat &lamToU = *GetLocMatrix(StdRegions::eHybridDGLamToU, varcoeffs, lambdaval, tau);

                    // Inverse mass matrix 
                    DNekScalMat &invMass = *GetLocMatrix(StdRegions::eInvMass);
                    
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
                    }

                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case StdRegions::eHybridDGLamToQ0:
                        dir = 0;
                        if(nvarcoeffs>0)
                        {
                            Dmat = GetLocMatrix(StdRegions::eWeakDirectionalDeriv,varcoeffs[dir]);
                        }

                        else
                        {
                            Dmat = GetLocMatrix(StdRegions::eWeakDeriv0);
                        }
                        break;
                    case StdRegions::eHybridDGLamToQ1:
                        dir = 1;
                        if(nvarcoeffs>0)
                        {
                            Dmat = GetLocMatrix(StdRegions::eWeakDirectionalDeriv,varcoeffs[dir]);
                        }

                        else
                        {
                            Dmat = GetLocMatrix(StdRegions::eWeakDeriv1);
                        }
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
                        
                        SetTraceToGeomOrientation(EdgeExp,lambda);
                        
                        // + \tilde{G} \lambda
                        AddNormTraceInt(dir,lambda,EdgeExp,f,mkey.GetVariableCoefficient(dir)); 
                        
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
                    int nbndry  = NumDGBndryCoeffs();
                    int coordim = GetCoordim();
                    int nedges  = GetNedges();
                    
                    Array<OneD,NekDouble>       work;
                    Array<OneD,const Array<OneD, NekDouble> > normals; 
                    Array<OneD, NekDouble> normalindir;
                    Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    Array<OneD, NekDouble> lam(nbndry); 
                    
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    Array<OneD,unsigned int>    emap;
                    Array<OneD, int>            sign;
                    StdRegions::EdgeOrientation edgedir;
                    
                    int nvarcoeffs = mkey.GetNvariableCoefficients();
                    Array<OneD, Array<OneD,const NekDouble> > varcoeffs(nvarcoeffs);
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
                    
                    DNekScalMatSharedPtr LamToQ[3];
                    
                    // Matrix to map Lambda to U
                    DNekScalMat &LamToU = *GetLocMatrix(StdRegions::eHybridDGLamToU, varcoeffs, lambdaval, tau);

                    // Matrix to map Lambda to Q0
                    LamToQ[0] = GetLocMatrix(StdRegions::eHybridDGLamToQ0, varcoeffs, lambdaval,tau);
 
                    // Matrix to map Lambda to Q1
                    LamToQ[1] = GetLocMatrix(StdRegions::eHybridDGLamToQ1, varcoeffs, lambdaval,tau);

                    // Matrix to map Lambda to Q2 for 3D coordinates
                    if (coordim == 3)
                    {
                        LamToQ[2] = GetLocMatrix(StdRegions::eHybridDGLamToQ2, varcoeffs, lambdaval,tau);
                    }

                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
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

                            normals = GetEdgeNormal(e);
                            edgedir = GetEorient(e);
                            
                            work = Array<OneD,NekDouble>(nquad_e);
                            
                            GetEdgeToElementMap(e,edgedir,emap,sign);

                            // Q0 * n0
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,sign[j]*(*LamToQ[0])(emap[j],i));
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
                                Vmath::Vmul(nquad_e,normals[0],1,EdgeExp[e]->GetPhys(),1,work,1);
                            }
                            
                            if(edgedir == StdRegions::eBackwards)
                            {
                                Vmath::Neg(nquad_e,work,1);
                            }

                            // Q1 * n1
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,sign[j]*(*LamToQ[1])(emap[j],i));
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
                                    Vmath::Vvtvp(nquad_e,normalindir,1,EdgeExp[e]->GetPhys(),1,work,1,work,1);
                                }

                                else
                                {
                                    Vmath::Vvtvp(nquad_e,normals[1],1,
                                                 EdgeExp[e]->GetPhys(),1,
                                                 work,1,work,1);
                                }
                            }
                            else // subtrace values for negative normal
                            {
                                if(nvarcoeffs>0)
                                {
                                    // normalindir = normals /cdot dirForcing[1]
                                    normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                    Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(1),normalindir);
                                    Vmath::Vvtvm(nquad_e,normalindir,1,EdgeExp[e]->GetPhys(),1,work,1,work,1);
                                    Vmath::Neg(nquad_e,work,1);
                                }
                                
                                else
                                {
                                    Vmath::Vvtvm(nquad_e,normals[1],1,
                                                 EdgeExp[e]->GetPhys(),1,
                                                 work,1,work,1);
                                    Vmath::Neg(nquad_e,work,1);
                                }
                            }

                            // Q2 * n2
                            if (coordim == 3)
                            {
                                for(j = 0; j < order_e; ++j)
                                {
                                    EdgeExp[e]->SetCoeff(j,sign[j]*(*LamToQ[2])(emap[j],i));
                                }
                                
                                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                     EdgeExp[e]->UpdatePhys());
    
                                if(edgedir == StdRegions::eForwards)
                                {
                                    if(nvarcoeffs>0)
                                    {
                                        // normalindir = normals /cdot dirForcing[1]
                                        normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                        Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(2),normalindir);
                                        Vmath::Vvtvp(nquad_e,normalindir,1,EdgeExp[e]->GetPhys(),1,work,1,work,1);
                                    }
    
                                    else
                                    {
                                        Vmath::Vvtvp(nquad_e,normals[2],1,
                                                     EdgeExp[e]->GetPhys(),1,
                                                     work,1,work,1);
                                    }
                                }
                                else // subtrace values for negative normal
                                {
                                    if(nvarcoeffs>0)
                                    {
                                        // normalindir = normals /cdot dirForcing[1]
                                        normalindir = Array<OneD, NekDouble>(nquad_e,0.0);
                                        Getnormalindir(e,EdgeExp[e],normals,mkey.GetVariableCoefficient(2),normalindir);
                                        Vmath::Vvtvm(nquad_e,normalindir,1,EdgeExp[e]->GetPhys(),1,work,1,work,1);
                                        Vmath::Neg(nquad_e,work,1);
                                    }
                                    
                                    else
                                    {
                                        Vmath::Vvtvm(nquad_e,normals[2],1,
                                                     EdgeExp[e]->GetPhys(),1,
                                                     work,1,work,1);
                                        Vmath::Neg(nquad_e,work,1);
                                    }
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
      void Expansion2D::v_DGDeriv(int dir,
                                const Array<OneD, const NekDouble>&incoeffs,
                                Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,
                                Array<OneD, NekDouble> &out_d)
      {
        StdRegions::MatrixType DerivType[3] = {StdRegions::eWeakDeriv0,
                                               StdRegions::eWeakDeriv1,
                                               StdRegions::eWeakDeriv2};
        
          int ncoeffs = GetNcoeffs();
          int nedges  = GetNedges();

#if 1
          DNekScalMat &InvMass = *GetLocMatrix(StdRegions::eInvMass);
          DNekScalMat &Dmat    = *GetLocMatrix(DerivType[dir]);
          
          Array<OneD, NekDouble> coeffs = incoeffs;
          DNekVec     Coeffs  (ncoeffs,coeffs, eWrapper);
          
          Coeffs = Transpose(Dmat)*Coeffs;
          Vmath::Neg(ncoeffs, coeffs,1);

          // Add the boundary integral including the relevant part of
          // the normal
          AddNormTraceInt(dir,EdgeExp,coeffs);
        
          DNekVec Out_d (ncoeffs,out_d,eWrapper);

          Out_d  = InvMass*Coeffs;
#else
          StdRegions::MatrixType LamToQType[3] = {StdRegions::eHybridDGLamToQ0,
                                                  StdRegions::eHybridDGLamToQ1,
                                                  StdRegions::eHybridDGLamToQ2};
          
          DNekScalMat &LamToQ = *GetLocMatrix(LamToQType[dir],0.0,1.0);
          int nedgetot = 0;
          for(i = 0; i < EdgeExp.num_elements(); ++i)
          {
              nedgetot += EdgeExp[i]->GetNcoeffs();
          }
          
          Array< OneD, NekDouble > lam(nedgetot), e_tmp;
          int nedgecoeffs = 0;
          for(i = 0; i < EdgeExp.num_elements(); ++i)
          {
              EdgeExp[i]->SetCoeffsToOrientation(v_GetCartersianEorient(i),ExpExp[i]->GetCoeffs(),e_tmp = lam + nedgecoeffs  );
              nedgecoeffs += EdgeExp[i]->GetNcoeffs();
          }

          DNekVec Lam(nedgetot,lam,eWrapper);
          DNekVec Out_d (ncoeffs,out_d,eWrapper);

          Out_d  = LamToQ*LamCoeffs;
#endif
      }

        enum BndToLocMatrixMapType
        {
            eBndToFullMatrixCG,
            eBndToBndMatrixCG,
            eBndToTraceMatrixDG
        };

        void Expansion2D::v_AddRobinMassMatrix(const int edge, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");
            ASSERTL1(inoutmat->GetRows() == inoutmat->GetColumns(),
                     "Assuming that input matrix was square");
            int i,j;
            int id1,id2;
            int order_e = m_edgeExp[edge]->GetNcoeffs();                    
         
            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;
            
            LocalRegions::MatrixKey mkey(StdRegions::eMass,StdRegions::eSegment, *m_edgeExp[edge], primCoeffs);
            DNekScalMat &edgemat = *m_edgeExp[edge]->GetLocMatrix(mkey);

            // Now need to identify a map which takes the local edge
            // mass matrix to the matrix stored in inoutmat;
            // This can currently be deduced from the size of the matrix
            
            // - if inoutmat.m_rows() == v_NCoeffs() it is a full
            //   matrix system
            
            // - if inoutmat.m_rows() == v_NumBndCoeffs() it is a
            //  boundary CG system

            // - if inoutmat.m_rows() == v_NumDGBndCoeffs() it is a
            //  trace DG system
            int rows = inoutmat->GetRows();

            if (rows == GetNcoeffs())
            {
                GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);
            }
            else if(rows == NumBndryCoeffs())
            {
                int nbndry = NumBndryCoeffs();
                Array<OneD,unsigned int> bmap(nbndry);

                GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);

                GetBoundaryMap(bmap);
                
                for(i = 0; i < order_e; ++i)
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
                map  = Array<OneD, unsigned int> (order_e);
                sign = Array<OneD, int> (order_e,1);
                
                for(i = 0; i < edge; ++i)
                {
                    cnt += GetEdgeNcoeffs(i);
                }
                
                for(i = 0; i < order_e; ++i)
                {
                    map[i] = cnt++;
                }
                // check for mapping reversal 
                if(GetEorient(edge) == StdRegions::eBackwards)
                {
                    switch(m_edgeExp[edge]->GetBasis(0)->GetBasisType())
                    {
                    case LibUtilities::eGLL_Lagrange:
                        reverse( map.get() , map.get()+order_e);
                        break;
                    case LibUtilities::eModified_A:
                        {
                            swap(map[0],map[1]);
                            for(i = 3; i < order_e; i+=2)
                            {
                                sign[i] = -1;
                            }  
                        }
                        break;
                    default:
                        ASSERTL0(false,"Edge boundary type not valid for this method");
                    }
                }
            }
            else
            {
                ASSERTL0(false,"Could not identify matrix type from dimension");
            }

            for(i = 0; i < order_e; ++i)
            {
                id1 = map[i];
                for(j = 0; j < order_e; ++j)
                {
                    id2 = map[j];
                    (*inoutmat)(id1,id2) +=  edgemat(i,j)*sign[i]*sign[j];
                }
            }
        }

        /**
         * Given an edge and vector of element coefficients:
         * - maps those elemental coefficients corresponding to the edge into
         *   an edge-vector.
         * - resets the element coefficients
         * - multiplies the edge vector by the edge mass matrix
         * - maps the edge coefficients back onto the elemental coefficients
         */
        void Expansion2D::v_AddRobinEdgeContribution(const int edgeid, const Array<OneD, const NekDouble> &primCoeffs, Array<OneD, NekDouble> &coeffs)
        {
            ASSERTL1(IsBoundaryInteriorExpansion(),
                     "Not set up for non boundary-interior expansions");
            int i,j;
            int order_e = m_edgeExp[edgeid]->GetNcoeffs();

            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;

            LocalRegions::MatrixKey mkey(StdRegions::eMass,StdRegions::eSegment, *m_edgeExp[edgeid], primCoeffs);
            DNekScalMat &edgemat = *m_edgeExp[edgeid]->GetLocMatrix(mkey);

            NekVector<NekDouble> vEdgeCoeffs (order_e);

            GetEdgeToElementMap(edgeid,v_GetEorient(edgeid),map,sign);

            for (i = 0; i < order_e; ++i)
            {
                vEdgeCoeffs[i] = coeffs[map[i]]*sign[i];
            }
            Vmath::Zero(GetNcoeffs(), coeffs, 1);

            vEdgeCoeffs = edgemat * vEdgeCoeffs;

            for (i = 0; i < order_e; ++i)
            {
                coeffs[map[i]] = vEdgeCoeffs[i]*sign[i];
            }

        }

    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion2D.cpp,v $
 *    Revision 1.24  2010/02/19 13:58:07  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.23  2010/01/11 20:52:16  cantwell
 *    Fixed HDGHelmholtz solver for embedded 2D expansion in 3D coordinate system.
 *    Fixed inconsistent MatrixTypeMap in StdRegions.
 *
 *    Revision 1.22  2009/12/17 17:48:22  bnelson
 *    Fixed visual studio compiler warning.
 *
 *    Revision 1.21  2009/12/15 18:09:02  cantwell
 *    Split GeomFactors into 1D, 2D and 3D
 *    Added generation of tangential basis into GeomFactors
 *    Updated ADR2DManifold solver to use GeomFactors for tangents
 *    Added <GEOMINFO> XML session section support in MeshGraph
 *    Fixed const-correctness in VmathArray
 *    Cleaned up LocalRegions code to generate GeomFactors
 *    Removed GenSegExp
 *    Temporary fix to SubStructuredGraph
 *    Documentation for GlobalLinSys and GlobalMatrix classes
 *
 *    Revision 1.20  2009/11/17 17:43:36  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.19  2009/11/16 16:27:48  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.18  2009/11/16 13:43:01  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.17  2009/11/13 16:18:34  sehunchun
 *    *** empty log message ***
 *
 *    Revision 1.16  2009/11/10 19:04:24  sehunchun
 *    Variable coefficients for HDG2D Solver
 *
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
