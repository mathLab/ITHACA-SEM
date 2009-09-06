///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.cpp
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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList1D.h>
#include <LibUtilities/Polylib/Polylib.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
        ExpList1D::ExpList1D():
            ExpList()
        {
        }
        
        ExpList1D::~ExpList1D()
        {
        }
        
        ExpList1D::ExpList1D(const ExpList1D &In):
            ExpList(In)
        {
        }

        ExpList1D::ExpList1D(const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraph1D &graph1D):
            ExpList()
        {
            int i,j, id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionVector &expansions = graph1D.GetExpansions();
            m_coeff_offset =  Array<OneD, int> (expansions.size());
            m_phys_offset  =  Array<OneD, int> (expansions.size());
            
            for(i = 0; i < expansions.size(); ++i)
            {                
                if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(expansions[i]->m_GeomShPtr))
                {
                    seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(Ba,SegmentGeom);
                    seg->SetElmtId(id++);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }  
            
                m_coeff_offset[i] = m_ncoeffs;
                m_phys_offset[i]  = m_npoints;
                m_ncoeffs += Ba.GetNumModes();
                m_npoints += Ba.GetNumPoints();
            } 
            
            ExpList::SetCoeffPhys();
        }

        ExpList1D::ExpList1D(SpatialDomains::MeshGraph1D &graph1D):
            ExpList()
        {
            int i,id=0;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionVector &expansions = graph1D.GetExpansions();
            m_coeff_offset = Array<OneD,int>(expansions.size());
            m_phys_offset  = Array<OneD,int>(expansions.size());
            
            for(i = 0; i < expansions.size(); ++i)
            {
                LibUtilities::BasisKey bkey = expansions[i]->m_BasisKeyVector[0];
                
                if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(expansions[i]->m_GeomShPtr))
                {
                    seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);
                    seg->SetElmtId(id++);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }  

                m_coeff_offset[i] = m_ncoeffs;
                m_phys_offset[i]  = m_npoints;
                m_ncoeffs += bkey.GetNumModes();
                m_npoints += bkey.GetNumPoints();
            }
      
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
            
        }
    
        ExpList1D::ExpList1D(const SpatialDomains::CompositeVector &domain, SpatialDomains::MeshGraph2D &graph2D):
            ExpList()
        {
            int i,j, nel,cnt,id=0;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;
            
            nel = 0;
            for(i = 0; i < domain.size(); ++i)
            {
                nel += (domain[i])->size();
            }

            m_coeff_offset = Array<OneD,int>(nel,0);
            m_phys_offset  = Array<OneD,int>(nel,0);

            cnt = 0;
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];
                
                for(j = 0; j < comp->size(); ++j)
                {                    
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        LibUtilities::BasisKey bkey = graph2D.GetEdgeBasisKey(SegmentGeom);
                        seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);
                        
                        seg->SetElmtId(id++);
                        (*m_exp).push_back(seg);  

                        m_coeff_offset[cnt] = m_ncoeffs; 
                        m_phys_offset[cnt]  = m_npoints; cnt++;
                        m_ncoeffs += bkey.GetNumModes();
                        m_npoints += bkey.GetNumPoints();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }  
                }
                
            } 
            
            ExpList::SetCoeffPhys();
        }

        ExpList1D::ExpList1D(const Array<OneD,const ExpList1DSharedPtr>  &bndConstraint, 
                             const Array<OneD, const SpatialDomains::BoundaryConditionType>  &bndTypes, 
                             const StdRegions::StdExpansionVector &locexp, SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j,k,cnt,id, elmtid=0;
            Array<OneD, int> EdgeDone(graph2D.GetNseggeoms(),0);
            
            // First loop over boundary conditions to renumber
            // Dirichlet boundaries
            cnt = 0;
            for(i = 0; i < bndTypes.num_elements(); ++i)
            {
                if(bndTypes[i] == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey = bndConstraint[i]->GetExp(j)->GetBasis(0)->GetBasisKey();
                        const SpatialDomains::Geometry1DSharedPtr& SegGeom = bndConstraint[i]->GetExp(j)->GetGeom1D();
                            
                        LocalRegions::SegExpSharedPtr Seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegGeom);
                            
                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);   

                        EdgeDone[SegGeom->GetEid()] = 1;
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
                        
                    if(!EdgeDone[id])
                    {
                        LibUtilities::BasisKey EdgeBkey = locexp[i]->DetEdgeBasisKey(j);
                        
                        LocalRegions::SegExpSharedPtr Seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(EdgeBkey, SegGeom);
                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);
                        
                        EdgeDone[id] = 1;
                    }
                }
            }

            // Set up offset information and array sizes
            ExpList::SetCoeffPhys();
            
        }


        void ExpList1D::SetBoundaryConditionExpansion(const SpatialDomains::MeshGraph1D &graph1D,
                                                      SpatialDomains::BoundaryConditions &bcs, 
                                                      const std::string variable,
                                                      Array<OneD, LocalRegions::PointExpSharedPtr> &bndCondExpansions,
                                                      Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions)
        {
            int i,j,k;
            int cnt  = 0;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();   
            
            LocalRegions::PointExpSharedPtr          locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond; 
            SpatialDomains::VertexComponentSharedPtr vert;

            int nbnd = bregions.size(); 
            
            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {  
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {       
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {
                            if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[j])[k]))
                            {
                                locPointExp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }
                        }
                    }
                } 
            } // end if Dirichlet
            
            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {    
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {     
                            if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[j])[k]))
                            {
                                locPointExp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }            
                        }
                    }
                }    
                else if((locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet) && 
                        (locBCond->GetBoundaryConditionType() != SpatialDomains::ePeriodic))
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }                  
            }
        }

        void ExpList1D::GetPeriodicVertices(const SpatialDomains::MeshGraph1D &graph1D,
                                              SpatialDomains::BoundaryConditions &bcs, 
                                              const std::string variable,
                                              map<int,int>& periodicVertices)
        {

            int i,j,k;
            
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int region1ID;
            int region2ID;

            SpatialDomains::Composite comp1;
            SpatialDomains::Composite comp2;

            SpatialDomains::VertexComponentSharedPtr vert1;
            SpatialDomains::VertexComponentSharedPtr vert2;
            
            SpatialDomains::BoundaryConditionShPtr locBCond; 

            // This std::map is a check so that the periodic pairs
            // are not treated twice
            map<int, int> doneBndRegions;

            int nbnd = bregions.size();
          
            for(i = 0; i < nbnd; ++i)
            {        
                locBCond = (*(bconditions[i]))[variable];  
                if(locBCond->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {    
                    region1ID = i;
                    region2ID = (boost::static_pointer_cast<SpatialDomains::PeriodicBoundaryCondition>(locBCond))->m_ConnectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {                    
                        ASSERTL0(bregions[region1ID]->size() == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions should be equal");
                    
                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];
                            
                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should be equal");
                            
                            for(k = 0; k < comp1->size(); k++)
                            {                                      
                                if(!(vert1 = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*comp1)[k]))||
                                   !(vert2 = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a VertexComponent failed");
                                } 

                                // Extract the periodic vertices
                                periodicVertices[vert1->GetVid()] = vert2->GetVid();
                            }
                        }
                    }
                    else
                    {
                        ASSERTL0(doneBndRegions[region1ID]==region2ID,
                                 "Boundary regions are not mutually periodic");
                    }
                    doneBndRegions[region2ID] = region1ID;
                }                  
            }
        }

        void ExpList1D::EvaluateBoundaryConditions(const NekDouble time,
                                                   Array<OneD, LocalRegions::PointExpSharedPtr> &bndCondExpansions,
                                                   Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndConditions)
        {
            int i;

            NekDouble x0;
            NekDouble x1;
            NekDouble x2;
            
            for(i = 0; i < bndCondExpansions.num_elements(); ++i)
            {
                bndCondExpansions[i]->GetCoords(x0,x1,x2);
                
                if(bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                { 
                    bndCondExpansions[i]->SetValue((boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(bndConditions[i])->
                                                    m_DirichletCondition).Evaluate(x0,x1,x2,time));
                }
                else if(bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                { 
                    bndCondExpansions[i]->SetValue((boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(bndConditions[i])->
                                                      m_NeumannCondition).Evaluate(x0,x1,x2,time));
                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }                
        }

		void ExpList1D::PostProcess(LibUtilities::KernelSharedPtr kernel, 
									Array<OneD,NekDouble> &inarray,
									Array<OneD,NekDouble> &outarray, 
									int elmId,
									NekDouble h)

		{
			int i,j,r;
			
			// get the local element expansion of the elmId element
			StdRegions::StdExpansionSharedPtr elmExp = GetExp(elmId);

			// Get the quadrature points and weights required for integration
			int quad_npoints = elmExp->GetTotPoints();
			LibUtilities::PointsKey quadPointsKey(quad_npoints, elmExp->GetPointsType(0));
			Array<OneD,NekDouble> quad_points = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			Array<OneD,NekDouble> quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();
			
			// Declare variable for the local kernel breaks
			int kernel_width = kernel->GetKernelWidth();
			Array<OneD,NekDouble> local_kernel_breaks(kernel_width+1);

			// Declare variable for the transformed quadrature points
			Array<OneD,NekDouble> mapped_quad_points(quad_npoints);
			
			// For each evaluation point
			for(i = 0; i < inarray.num_elements(); i++)
			{
				// Move the center of the kernel to the current point
				kernel->MoveKernelCenter(inarray[i],local_kernel_breaks);
							
				// Find the mesh breaks under the kernel support
				Array<OneD,NekDouble> mesh_breaks;
				kernel->FindMeshUnderKernel(local_kernel_breaks,h,mesh_breaks);
				
				// Sort the total breaks for integration purposes
				int total_nbreaks = local_kernel_breaks.num_elements() + 
									mesh_breaks.num_elements(); // number of the total breaks
				Array<OneD,NekDouble> total_breaks(total_nbreaks);
				kernel->Sort(local_kernel_breaks,mesh_breaks,total_breaks);
				
				// Integrate the product of kernel and function over the total breaks
				NekDouble integral_value = 0.0;
				for(j = 0; j < total_breaks.num_elements()-1; j++)
				{
					double a = total_breaks[j];
					double b = total_breaks[j+1];
					
					// Map the quadrature points to the appropriate interval
					for(r = 0; r < quad_points.num_elements(); r++)
					{
						mapped_quad_points[r] = (quad_points[r]+1.0)*0.5*(b-a)+a;
					}

					// Evaluate the function at the transformed quadrature points
					Array<OneD,NekDouble> u_value(quad_npoints);
					Array<OneD,NekDouble> coeffs = GetCoeffs();

					//PeriodicEval(mapped_quad_points,h,u->GetExpSize(),elmExp->GetBasisNumModes(0),coeffs,u_value);
					PeriodicEval(coeffs,mapped_quad_points,h,elmExp->GetBasisNumModes(0),u_value);
					
					// Evaluate the kernel at the transformed quadrature points
					Array<OneD,NekDouble> k_value(quad_npoints);
					kernel->EvaluateKernel(mapped_quad_points,h,k_value);

					// Integrate
					for(r = 0; r < quad_npoints; r++)
					{
						integral_value += (b-a)*0.5*k_value[r]*u_value[r]*quad_weights[r];
					}
				}
				outarray[i] = integral_value/h;
			}
		}

		void ExpList1D::PostProcess( LibUtilities::KernelSharedPtr kernel,
									 Array<OneD,NekDouble> &inarray,
									 Array<OneD,NekDouble> &outarray,
									 NekDouble h)
		{
			int i,j,r;
		
			// Get the quadrature points and weights required for integration
			int quad_npoints = (GetExp(0))->GetTotPoints();
			LibUtilities::PointsKey quadPointsKey(quad_npoints, GetExp(0)->GetPointsType(0));
			Array<OneD,NekDouble> quad_points = LibUtilities::PointsManager()[quadPointsKey]->GetZ();
			Array<OneD,NekDouble> quad_weights = LibUtilities::PointsManager()[quadPointsKey]->GetW();

			// Declare variable for the local kernel breaks
			int kernel_width = kernel->GetKernelWidth();
			Array<OneD,NekDouble> local_kernel_breaks(kernel_width+1);

			// Declare variable for the transformed quadrature points
			Array<OneD,NekDouble> mapped_quad_points(quad_npoints);
			
			// For each evaluation point
			for(i = 0; i < inarray.num_elements(); i++)
			{
				// Move the center of the kernel to the current point
				kernel->MoveKernelCenter(inarray[i],local_kernel_breaks);
							
				// Find the mesh breaks under the kernel support
				Array<OneD,NekDouble> mesh_breaks;
				kernel->FindMeshUnderKernel(local_kernel_breaks,h,mesh_breaks);
				
				// Sort the total breaks for integration purposes
				int total_nbreaks = local_kernel_breaks.num_elements() + 
									mesh_breaks.num_elements(); // number of the total breaks
				Array<OneD,NekDouble> total_breaks(total_nbreaks);
				kernel->Sort(local_kernel_breaks,mesh_breaks,total_breaks);
				
				// Integrate the product of kernel and function over the total breaks
				NekDouble integral_value = 0.0;
				for(j = 0; j < total_breaks.num_elements()-1; j++)
				{
					double a = total_breaks[j];
					double b = total_breaks[j+1];
					
					// Map the quadrature points to the appropriate interval
					for(r = 0; r < quad_points.num_elements(); r++)
					{
						mapped_quad_points[r] = (quad_points[r]+1.0)*0.5*(b-a)+a;
					}

					// Evaluate the function at the transformed quadrature points
					Array<OneD,NekDouble> u_value(quad_npoints);
					Array<OneD,NekDouble> coeffs = GetCoeffs();

					//PeriodicEval(mapped_quad_points,h,u->GetExpSize(),(u->GetExp(0))->GetBasisNumModes(0),coeffs,u_value);
					PeriodicEval(coeffs,mapped_quad_points,h,(GetExp(0))->GetBasisNumModes(0),u_value);
					
					// Evaluate the kernel at the transformed quadrature points
					Array<OneD,NekDouble> k_value(quad_npoints);
					kernel->EvaluateKernel(mapped_quad_points,h,k_value);

					// Integrate
					for(r = 0; r < quad_npoints; r++)
					{
						integral_value += (b-a)*0.5*k_value[r]*u_value[r]*quad_weights[r];
					}
					
				}
				
				outarray[i] = integral_value/h;
			}
			
		}

		void ExpList1D::PeriodicEval(Array<OneD,NekDouble> &inarray1, Array<OneD,NekDouble> &inarray2,
									 NekDouble h, int nmodes,
									 Array<OneD,NekDouble> &outarray)
		{
			int i,j,r;

			// Get the number of elements in the domain
			int num_elm = GetExpSize();

			// initializing the outarray
			for(i = 0; i < outarray.num_elements(); i++)
			{
				outarray[i] = 0.0;
			}
			
			// Make a copy for further modification
			int x_size = inarray2.num_elements();
			Array<OneD,NekDouble> x_values_cp(x_size);
			
			// Determining the element to which the x belongs
			Array<OneD,int> x_elm(x_size);
			for(i = 0; i < x_size; i++ )
			{
				x_elm[i] = floor(inarray2[i]/h);
			}

			// Clamp indices periodically
			for(i = 0; i < x_size; i++)
			{
				while(x_elm[i] < 0)
				{
					x_elm[i] += num_elm;
				}
				while(x_elm[i] >= num_elm)
				{
					x_elm[i] -= num_elm ;
				}
			}
			
			// Map the values of x to [-1 1] on its interval
			for(i = 0; i < x_size; i++)
			{
				x_values_cp[i] = (inarray2[i]/h - floor(inarray2[i]/h))*2 - 1.0;
			}

			//Evaluate the jocobi polynomials
			// (Evaluating the base at some points other than the quadrature points)
			// Should it be added to the base class????
			Array<TwoD,NekDouble> jacobi_poly(nmodes,x_size);
			for(i = 0; i < nmodes; i++)
			{	
				Polylib::jacobfd(x_size,x_values_cp.get(),jacobi_poly.get()+i*x_size,NULL,i,0.0,0.0);
			}

			// Evaluate the function values
			for(r = 0; r < nmodes; r++)
			{
				for(j = 0; j < x_size; j++)
				{
					int index = ((x_elm[j])*nmodes)+r;
					outarray[j] += inarray1[index]*jacobi_poly[r][j];
				}
			}

		}


    } //end of namespace
} //end of namespace

/**
* $Log: ExpList1D.cpp,v $
* Revision 1.36  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.35  2009/02/08 09:11:49  sherwin
* General updates to introduce multiple matrix definitions based on different boundary types
*
* Revision 1.34  2009/01/13 02:50:10  mirzaee
* Added definitions for the PostProcessing functions and PeriodicEval
*
* Revision 1.33  2008/09/09 15:06:03  sherwin
* Modifications related to curved elements.
*
* Revision 1.32  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.31  2008/07/31 11:17:13  sherwin
* Changed GetEdgeBasis with DetEdgeBasisKey
*
* Revision 1.30  2008/07/29 22:27:33  sherwin
* Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
*
* Revision 1.29  2008/07/12 17:31:39  sherwin
* Added m_phys_offset and rename m_exp_offset to m_coeff_offset
*
* Revision 1.28  2008/06/23 14:21:01  pvos
* updates for 1D ExpLists
*
* Revision 1.27  2008/05/14 18:06:50  sherwin
* mods to fix Seggeom to Geometry1D casting
*
* Revision 1.26  2008/05/13 22:06:58  sherwin
* Changed SegGeom to Geometry1D
*
* Revision 1.25  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.24  2008/03/18 14:14:13  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.23  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.22  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.21  2007/11/07 20:29:53  jfrazier
* Modified to use new expansion list contained in meshgraph.
*
* Revision 1.20  2007/09/25 14:25:29  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.19  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.18  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.17  2007/07/13 16:48:47  pvos
* Another HelmHoltz update (homogeneous dir BC multi-elemental solver does work)
*
* Revision 1.16  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.15  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.14  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
