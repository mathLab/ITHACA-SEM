///////////////////////////////////////////////////////////////////////////////
//
// File GenSegExp.cpp
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
// Description: GenSegExp routines
//
///////////////////////////////////////////////////////////////////////////////
#include <LocalRegions/GenSegExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        GenSegExp::GenSegExp(const LibUtilities::BasisKey &Ba, 
                             const SpatialDomains::Geometry1DSharedPtr &geom):
            SegExp(Ba,geom)
        {
        }

        GenSegExp::GenSegExp(const GenSegExp &S):
            SegExp(S),
            m_physNormal(S.m_physNormal),
            m_physBiNormal(S.m_physBiNormal)
        {
        }
        // by default the StdExpansion1D destructor will be called        
        GenSegExp::~GenSegExp()
        {
        }
        
        // Given a Geometry2D set normals to the same as the values
        // along edge \a edge. If NegateNormals is true then negate the
        // normals
        void GenSegExp::SetUpPhysNormals(const SpatialDomains::Geometry2DSharedPtr& Geom, const int edge, bool NegateNormals)
        {            
            int k;
            int coordim = Geom->GetCoordim();
            Array<TwoD, const NekDouble> normals = Geom->GetGeomFactors()->GetNormals();
            SpatialDomains::GeomType Gtype = Geom->GetGtype();

            // assume all coordinates have same basis
            LibUtilities::BasisSharedPtr CBasis = Geom->GetEdgeBasis(0,edge); 
            int nq_orig  = CBasis->GetNumPoints();
            int nq       = m_base[0]->GetNumPoints();
            
            m_physNormal = Array<OneD, NekDouble> (coordim*nq);
            
            if(Gtype == SpatialDomains::eDeformed)
            {
                Array<OneD,NekDouble> tmpnorm(nq_orig);
                
                for(k = 0;k < coordim; ++k)
                {
                    if(edge >= 2)
                    {
                        Vmath::Reverse(nq_orig,&(normals[edge])[k*nq_orig],
                                       1, &tmpnorm[0],1);
                    }
                    else
                    {
                        Vmath::Vcopy(nq_orig,&(normals[edge])[k*nq_orig],
                                     1,&tmpnorm[0],1);
                    }
                    
                    Interp1D(CBasis->GetBasisKey(), &tmpnorm[0],
                             GetBasis(0)->GetBasisKey(),
                             &m_physNormal[k*nq]);

                    if(edge >= 2)
                    {
                        Vmath::Reverse(nq,&m_physNormal[k*nq],1,
                                       &m_physNormal[k*nq],1);
                    }
                }
            }
            else
            {                                        
                for(k = 0;k < coordim; ++k)
                {
                    Vmath::Fill(nq,normals[edge][k],&m_physNormal[k*nq],1);
                }
            }
            
            if(NegateNormals)
            {
                Vmath::Neg(nq*coordim,m_physNormal,1);
            }
        }

    } // end of namespace    
}//end of namespace

// $Log: GenSegExp.cpp,v $
// Revision 1.1  2008/07/29 22:24:49  sherwin
// Generalised Segment expansion which include a normal and binormal at the physical quadrature points
//
