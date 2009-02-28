///////////////////////////////////////////////////////////////////////////////
//
// File GenExpList1D.cpp
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
// Description: Generalised Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GenExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
        GenExpList1D::GenExpList1D():
            ExpList1D()
        {
        }
        
        GenExpList1D::~GenExpList1D()
        {
        }
        
        GenExpList1D::GenExpList1D(const GenExpList1D &In):
            ExpList1D(In)
        {
        }


        GenExpList1D::GenExpList1D(const Array<OneD,const MultiRegions::ExpList1DSharedPtr> &bndConstraint,  
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond, 
                                   const StdRegions::StdExpansionVector &locexp, 
                                   SpatialDomains::MeshGraph2D &graph2D,
                                   const map<int,int> &periodicEdges)
        {
            int i,j,id, elmtid=0, edgeid;
            // Could replace these with maps
            map<int, int> EdgeDone;
            map<int, int> NormalSet;
            
            StdRegions::EdgeOrientation      orient;
            LocalRegions::GenSegExpSharedPtr Seg;
            SpatialDomains::ElementEdgeVectorSharedPtr con_elmt;
            SpatialDomains::Geometry2DSharedPtr AdjGeom;
            SpatialDomains::Geometry1DSharedPtr SegGeom;

            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            for(i = 0; i < bndCond.num_elements(); ++i)
            {
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey = bndConstraint[i]->GetExp(j)->GetBasis(0)->GetBasisKey();
                        SegGeom  = bndConstraint[i]->GetExp(j)->GetGeom1D();  
                        Seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(bkey, SegGeom);
                        con_elmt = graph2D.GetElementsFromEdge(SegGeom);
                        
                        EdgeDone[SegGeom->GetEid()] = elmtid;
                        
                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);   
                        
                    }
                }
            }
            
            // loop over all other edges and fill out other connectivities
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < locexp[i]->GetNedges(); ++j)
                {   
                    const SpatialDomains::Geometry1DSharedPtr& SegGeom = (locexp[i]->GetGeom2D())->GetEdge(j);
                    
                    id = SegGeom->GetEid();
                    
                    if(EdgeDone.count(id) == 0)
                    {
                        
                        LibUtilities::BasisKey EdgeBkey = locexp[i]->DetEdgeBasisKey(j);
                        
                        Seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(EdgeBkey, SegGeom);
                        
                        EdgeDone[id] = elmtid;

                        // set periodic edge 
                        if(periodicEdges.count(id) > 0)
                        {
                            EdgeDone[periodicEdges.find(id)->second] = elmtid;
                        }

                        Seg->SetElmtId(elmtid++);
                        
                        (*m_exp).push_back(Seg);
                    }

                    if(NormalSet.count(id) == 0)
                    {
                        Seg = boost::dynamic_pointer_cast<LocalRegions::GenSegExp>((*m_exp)[EdgeDone.find(id)->second]);
                        
                        // Set up normals at all Segment Quadrature points
                        Seg->SetUpPhysNormals(locexp[i],j);
                        
                        NormalSet[id] = 1;
                    }
                }
            }

            ExpList::SetCoeffPhys();
        }


        

        // Upwind the left and right states given by the Arrays Fwd
        // and Bwd using the vector quantity Vec and ouput the
        // upwinded value in the array upwind
        void GenExpList1D::Upwind(const Array<OneD, const Array<OneD, NekDouble> > &Vec,
                                  const Array<OneD, const NekDouble> &Fwd, 
                                  const Array<OneD, const NekDouble> &Bwd, 
                                  Array<OneD, NekDouble> &Upwind,
								  int direction)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals; 
            NekDouble Vn;

            // Assume whole array is of same coordimate dimention
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            ASSERTL1(Vec.num_elements() >= coordim,
              "Input vector does not have sufficient dimensions to match coordim");

            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                normals   = (*m_exp)[i]->GetPhysNormals();
                
                offset = m_phys_offset[i];

                for(j = 0; j < e_npoints; ++j)
                {
                    // Calculate normal velocity
                    Vn = 0.0;
                    for(k = 0; k < coordim; ++k)
                    {
                        Vn += Vec[k][offset+j]*normals[k*e_npoints + j];
                    }

                    // Upwind
                    if(direction*Vn > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }            
        }

        void GenExpList1D::Upwind(const Array<OneD, const NekDouble> &Vn, 
                                  const Array<OneD, const NekDouble> &Fwd, 
                                  const Array<OneD, const NekDouble> &Bwd, 
                                  Array<OneD, NekDouble> &Upwind,
								  int direction)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals; 

            // Assume whole array is of same coordimate dimention
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                offset = m_phys_offset[i];

                for(j = 0; j < e_npoints; ++j)
                {
                    // Upwind
                    if(direction*Vn[offset + j] > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }            
        }
		

        void GenExpList1D::GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals) 
        {
            
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> locnormals; 
            
            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to match coordim");
            
            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                locnormals   = (*m_exp)[i]->GetPhysNormals();
                
                offset = m_phys_offset[i];
                
                for(j = 0; j < e_npoints; ++j)
                {
                    for(k = 0; k < coordim; ++k)
                    {
                        normals[k][offset+j] = locnormals[k*e_npoints + j];
                    }
                }
            }
        }
        
    } //end of namespace
} //end of namespace

/**
* $Log: GenExpList1D.cpp,v $
* Revision 1.8  2008/10/29 22:46:35  sherwin
* Updates for const correctness and a few other bits
*
* Revision 1.7  2008/10/09 21:47:36  ehan
* Fixed error from the function Upwind().
*
* Revision 1.6  2008/10/04 19:54:28  sherwin
* Added upwind method using only the normal flux
*
* Revision 1.5  2008/09/09 15:06:03  sherwin
* Modifications related to curved elements.
*
* Revision 1.4  2008/08/22 09:42:32  pvos
* Updates for Claes' Shallow Water and Euler solver
*
* Revision 1.3  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.2  2008/07/31 11:17:13  sherwin
* Changed GetEdgeBasis with DetEdgeBasisKey
*
* Revision 1.1  2008/07/29 22:26:35  sherwin
* Generalised 1D Segment list which includes a normal direction at physical points
*
**/
