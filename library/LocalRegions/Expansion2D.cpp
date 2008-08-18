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

#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions 
    {
       

        void Expansion2D::AddEdgeNormBoundaryInt(const int edge, 
                                                 StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                                 Array<OneD, NekDouble> &Fx,  
                                                 Array<OneD, NekDouble> &Fy,  
                                                 Array<OneD, NekDouble> &outarray)
        {
            int i;
            int order_e = EdgeExp->GetNcoeffs();                    
            int nquad_e = EdgeExp->GetNumPoints(0);
            Array<OneD, const NekDouble> normals = EdgeExp->GetPhysNormals();
            Array<OneD,unsigned int>     map;
            Array<OneD,int>              sign;

            StdRegions::EdgeOrientation edgedir = v_GetEorient(edge);

            v_GetEdgeToElementMap(edge,edgedir,map,sign);

            ASSERTL1(v_GetCoordim() == 2,"Routine only set up for two-dimensions");
            
            Vmath::Vmul(nquad_e,&(normals[0]),1,&Fx[0],1,
                        &(EdgeExp->UpdatePhys())[0],1);
            Vmath::Vvtvp(nquad_e,&(normals[nquad_e]),1,
                         &Fy[0],1,&(EdgeExp->GetPhys())[0],1,
                         &(EdgeExp->UpdatePhys())[0],1);
            

            EdgeExp->IProductWRTBase(EdgeExp->GetPhys(),EdgeExp->UpdateCoeffs());
            // negate backward edge values due to inwards normal definition
            if(edgedir == StdRegions::eBackwards)
            {

                Vmath::Neg(order_e,EdgeExp->UpdateCoeffs(),1);
            }

            // add data to outarray if forward edge normal is outwards
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*EdgeExp->GetCoeff(i);
            }
        }

        // Given an array of trace expansion coefficients assuming all
        // edges are orientated in the eForwards direction re-order
        // entries to satisfy local edge orientation set by the
        // Geometry definition
        void Expansion2D::SetTraceToGeomOrientation(Array<OneD, NekDouble> &inout)
        {
            int i;
            int nedges  = v_GetNedges();
            Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    
            // Set up edge segment expansions from local geom info
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i] = v_GetEdgeExp(i,false);
            }

            SetTraceToGeomOrientation(EdgeExp,inout);
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
                                          Array<OneD,NekDouble> &outarray) 
        {
            int i,e,cnt;
            int order_e,nquad_e;
            int nedges = v_GetNedges();

            Array<OneD,NekDouble> normals;
            Array<OneD,StdRegions::StdExpansion1DSharedPtr> EdgeExp(nedges);
            
            // Set up edge segment expansions
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i] = v_GetEdgeExp(i);
            }

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


        void Expansion2D:: AddEdgeBoundaryInt(const int edge, 
                                              StdRegions::StdExpansion1DSharedPtr &EdgeExp,
                                              Array <OneD,NekDouble > &outarray)
        {
            
            int i;
            int order_e = EdgeExp->GetNcoeffs();                    
            Array<OneD,unsigned int> map;
            Array<OneD,int> sign;

            v_GetEdgeToElementMap(edge,v_GetEorient(edge),map,sign);

            EdgeExp->IProductWRTBase(EdgeExp->GetPhys(),EdgeExp->UpdateCoeffs());
            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[map[i]] += sign[i]*EdgeExp->GetCoeff(i);
            }
        }

        // Boundary terms associated with elemental Helmholtz matrix operations
        void Expansion2D::AddHDGHelmholtzMatrixBoundaryTerms(const NekDouble tau, 
                                                             const Array<OneD,
                                                             const NekDouble> &inarray,
                                                             Array<OneD,NekDouble> &outarray)
        {
            int i,e;
            int nbndry  = v_NumBndryCoeffs();
            int nquad0  = v_GetNumPoints(0);
            int nquad1  = v_GetNumPoints(1);
            int coordim = v_GetCoordim();
            int nedges  = v_GetNedges();

            Array<OneD, unsigned int>   emap;
            Array<OneD, int> sign;
            StdRegions::EdgeOrientation edgedir;

            Array<OneD,NekDouble>       in_phys(nquad0*nquad1);
            Array<OneD,StdRegions::StdExpansion1DSharedPtr >  EdgeExp(nedges);
            Array<OneD,Array<OneD,NekDouble> > deriv(3);
            
            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");

            // Set up edge segment expansions
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i] = v_GetEdgeExp(i);
            }

            //  Get physical solution. 
            v_BwdTrans(inarray,in_phys);

            // Calculate derivative for matrix terms.
            deriv[0] = Array<OneD,NekDouble>(nquad0*nquad1);
            deriv[1] = Array<OneD,NekDouble>(nquad0*nquad1);

            if(coordim == 2)
            {
                v_PhysDeriv(in_phys,deriv[0],deriv[1]);
            }
            else
            {
                deriv[2] = Array<OneD,NekDouble>(nquad0*nquad1);
                v_PhysDeriv(in_phys,deriv[0],deriv[1],deriv[2]);
            }

            // Loop over edges
            for(e = 0; e < nedges; ++e)
            {
                v_GetEdgePhysVals(e,EdgeExp[e],in_phys,EdgeExp[e]->UpdatePhys());
                AddHDGHelmholtzEdgeTerms(tau,e,EdgeExp,outarray);
                
                //=============================================================
                // Add -D^T M^{-1}G operation =-<n phi_i, n.d(in_phys)/dx]>
                //term which arise in matrix formulations but not rhs
                int nquad_e = EdgeExp[e]->GetNumPoints(0);                    
                int order_e = EdgeExp[e]->GetNcoeffs();                    
                Array<OneD,NekDouble> inval(nquad_e);
                Array<OneD,NekDouble> normals = EdgeExp[e]->GetPhysNormals();

                edgedir = v_GetEorient(e);
                v_GetEdgeToElementMap(e,edgedir,emap,sign);
                
                Vmath::Zero(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                
                for(i = 0; i < coordim; ++i)
                {
                    v_GetEdgePhysVals(e,EdgeExp[e],deriv[i],inval);
                    
                    Vmath::Vvtvp(nquad_e,&normals[i*nquad_e],1,&inval[0],1,
                                 &(EdgeExp[e]->UpdatePhys())[0],1,
                                 &(EdgeExp[e]->UpdatePhys())[0],1);
                }
                
                // negate backwards normal
                if(edgedir == StdRegions::eBackwards)
                {
                    Vmath::Neg(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                }

                // Fill edge and take inner product
                EdgeExp[e]->IProductWRTBase(EdgeExp[e]->GetPhys(),
                                            EdgeExp[e]->UpdateCoeffs());
                
                // Put data in out array
                for(i = 0; i < order_e; ++i)
                {
                    outarray[emap[i]] -= sign[i]*EdgeExp[e]->GetCoeff(i);
                }                    
            }
            //================================================================
        }

        //===============================================================
        // Boundary terms associated with elemental Helmholtz matrix
        // operations from the trace space
        void Expansion2D::AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                    const Array<OneD,const NekDouble> &inarray,
                                                    Array<OneD,NekDouble> &outarray)
        {
            int i;
            int nedges = v_GetNedges();
            Array<OneD, StdRegions::StdExpansion1DSharedPtr >  EdgeExp(nedges);
            
            // Set up edge segment expansions
            for(i = 0; i < nedges; ++i)
            {
                EdgeExp[i] = v_GetEdgeExp(i);
            }
            
            AddHDGHelmholtzTraceTerms(tau,inarray,EdgeExp,outarray);
        }
        
        
        // This method assumes that data in EdgeExp is orientated 
        // according to elemental counter clockwise format
        void Expansion2D::AddHDGHelmholtzTraceTerms(const NekDouble tau, 
                                                       const Array<OneD, const NekDouble> &inarray, Array<OneD,StdRegions::StdExpansion1DSharedPtr> &EdgeExp,  Array<OneD,NekDouble> &outarray)
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
                
                AddHDGHelmholtzEdgeTerms(tau,e,EdgeExp,outarray);
                
                cnt += order_e;
            }
        }


        void Expansion2D::AddHDGHelmholtzEdgeTerms(const NekDouble tau, 
                                                      const int edge,
                                                      Array <OneD, StdRegions::StdExpansion1DSharedPtr > &EdgeExp, 
                                                      Array <OneD,NekDouble > &outarray)
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
            for(n = 0; n < coordim; ++n)
            {
                //G;
                Vmath::Vmul(nquad_e,&normals[n*nquad_e],1,
                            &(EdgeExp[edge]->GetPhys())[0],1, &inval[0],1);

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
                
                switch(n)
                {
                case 0:
                    {
                        DNekScalMat &Dmat = *v_GetLocMatrix(StdRegions::eWeakDeriv0);
                        Coeffs = Coeffs  + Dmat*Tmpcoeff; 
                    }
                    break;
                case 1:
                    {
                        DNekScalMat &Dmat = *v_GetLocMatrix(StdRegions::eWeakDeriv1);
                        Coeffs = Coeffs  + Dmat*Tmpcoeff; 
                    }
                    break;
                case 2:
                    {
                        DNekScalMat &Dmat = *v_GetLocMatrix(StdRegions::eWeakDeriv2);
                        Coeffs = Coeffs  + Dmat*Tmpcoeff; 
                    }
                    break;
                }
            }
        }


        DNekMatSharedPtr Expansion2D::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmBndLam:
                {
                    int order_e, nquad_e;
                    int i,j,e,cnt;
                    int nbndry  = v_NumDGBndryCoeffs();
                    int coordim = v_GetCoordim();
                    int nedges  = v_GetNedges();
                    
                    Array<OneD,NekDouble>       work;
                    Array<OneD,const NekDouble> normals; 
                    Array<OneD,StdRegions::StdExpansion1DSharedPtr>  EdgeExp(nedges);
                    Array<OneD, NekDouble> lam(nbndry); 
                    
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    
                    Array<OneD,unsigned int>    emap;
                    Array<OneD, int>            sign;
                    StdRegions::EdgeOrientation edgedir;
                    
                    ASSERTL0(coordim < 3,"Needs to be  set up for expansion in 3 space");

                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    DNekScalMat &LamToU = *v_GetLocMatrix(StdRegions::eHybridDGLamToU, lambdaval,tau); 
                
                    // Matrix to map Lambda to Q0
                    DNekScalMat &LamToQ0 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ0, lambdaval,tau);

                    // Matrix to map Lambda to Q1
                    DNekScalMat &LamToQ1 = *v_GetLocMatrix(StdRegions::eHybridDGLamToQ1, lambdaval,tau);


                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < nedges; ++i)
                    {
                        EdgeExp[i] = v_GetEdgeExp(i);
                    }

                    // Set up matrix derived from <mu, Q_lam.n - \tau (
                    // U_lam - Lam) > 
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
          

                            Vmath::Vmul(nquad_e,normals,1,
                                        EdgeExp[e]->GetPhys(),1,
                                        work,1);

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
                                Vmath::Vvtvp(nquad_e,&normals[nquad_e],1,
                                             &EdgeExp[e]->GetPhys()[0],1,
                                             &work[0],1,&work[0],1);
                            }
                            else // subtrace values for negative normal
                            {
                                Vmath::Vvtvm(nquad_e,&normals[nquad_e],1,
                                             &EdgeExp[e]->GetPhys()[0],1,
                                             &work[0],1,&work[0],1);
                                Vmath::Neg(nquad_e,work,1);
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
                returnval = Expansion::GenMatrix(mkey);
                break;
            }
            
            return returnval;
        }
        
    } //end of namespace
} //end of namespace

/** 
 *    $Log: Expansion2D.cpp,v $
 *    Revision 1.1  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *
 **/
