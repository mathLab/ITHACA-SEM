///////////////////////////////////////////////////////////////////////////////
//
// File ContField3D.cpp
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
// Description: Field definition for 3D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

        ContField3D::ContField3D(void):
            ContExpList3D(),
            m_bndConstraint(),
            m_bndTypes()
        {
        }

        ContField3D::ContField3D(const ContField3D &In):
            ContExpList3D(In),
            m_bndConstraint(In.m_bndConstraint),
            m_bndTypes(In.m_bndTypes)
        {
        }

        ContField3D::ContField3D(SpatialDomains::MeshGraph3D &graph3D, SpatialDomains::BoundaryConditions &bcs, const int bc_loc):
            ContExpList3D(graph3D,false)  ,
            m_bndConstraint(),
            m_bndTypes()
        {
            GenerateField3D(graph3D, bcs, bcs.GetVariable(bc_loc));
            
             //TODO implement "LocalToGlobalMap3D"
             // m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,graph3D,m_bndConstraint, m_bndTypes);
        
        m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

         ContField3D::ContField3D(SpatialDomains::MeshGraph3D &graph3D, SpatialDomains::BoundaryConditions &bcs, const std::string variable):
            ContExpList3D(graph3D,false),
            m_bndConstraint(),
            m_bndTypes()
        {
            GenerateField3D(graph3D, bcs, variable);

               //TODO implement "LocalToGlobalMap3D"
            // m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp,graph3D, m_bndConstraint, m_bndTypes);
        
        m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc,
                                 const LibUtilities::PointsType TetNb):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,false)  ,
            m_bndConstraint(),
            m_bndTypes()
        {
            GenerateField3D(graph3D, bcs, bcs.GetVariable(bc_loc));

             //TODO implement "LocalToGlobalMap3D"
            // m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp, graph3D,m_bndConstraint, m_bndTypes);
        
        m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }

        ContField3D::ContField3D(const LibUtilities::BasisKey &Ba,
                                 const LibUtilities::BasisKey &Bb,
                                 const LibUtilities::BasisKey &Bc,
                                 SpatialDomains::MeshGraph3D &graph3D,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable,
                                 const LibUtilities::PointsType TetNb):
            ContExpList3D(Ba,Bb,Bc,graph3D,TetNb,false),
            m_bndConstraint(),
            m_bndTypes()
        {
            GenerateField3D(graph3D, bcs, variable);

             //TODO implement "LocalToGlobalMap3D"             
            // m_locToGloMap = MemoryManager<LocalToGlobalMap3D>::AllocateSharedPtr(m_ncoeffs,*m_exp, graph3D,m_bndConstraint,m_bndTypes);
        
        m_contNcoeffs = m_locToGloMap->GetTotGloDofs();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
        }


        ContField3D::~ContField3D()
        {
        }


//         void ContField3D::GenerateField3D(SpatialDomains::MeshGraph3D &graph3D,
//                                            SpatialDomains::BoundaryConditions &bcs,
//                                            const std::string variable)
//         {
//             int mycnt = 0;
// 
//             int i,j,m;
//             int cnt = 0;
//             int coeffcnt = 0;
//             int nbnd;
//             int nummodes;
//             int npoints;
//             SpatialDomains::BoundaryConditionShPtr locBCond;
//             SpatialDomains::BoundaryRegionShPtr locBregion;
//             SpatialDomains::SegGeomSharedPtr segGeom;         
//             LocalRegions::SegExpSharedPtr collSeg;
//             StdRegions::StdSegExpSharedPtr collStdSeg;
// 
//             SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
//             SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();
// 
//             MultiRegions::ExpList2DSharedPtr locExpList;
//             //LocalRegions::SegExpSharedPtr locSegExp;
//             LocalRegions::QuadExpSharedPtr locQuadExp;
//             LocalRegions::TriExpSharedPtr locTriExp;
// 
//             nbnd = bregions.size();
// 
//             m_bndConstraint = Array<OneD,MultiRegions::ExpList2DSharedPtr>(nbnd);
//             m_bndTypes = Array<OneD,SpatialDomains::BoundaryConditionType>(nbnd);
// 
//             // list Dirichlet boundaries first
//             for( i = 0; i < nbnd; ++i)
//             {
//                 locBCond = (*(bconditions[i]))[variable];
//                 if(locBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
//                 {
//                     locBregion = bregions[i];
// 
//                     //TODO: check this
//                     locExpList = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*locBregion, graph3D);
// 
//                     coeffcnt=0;
//                     for(j = 0; j < locExpList->GetExpSize(); j++)
//                     {
//                         locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>(locExpList->GetExp(j));
//                         locTriExp  = boost::dynamic_pointer_cast<LocalRegions::TriExp>(locExpList->GetExp(j));
//                         
//                         StdRegions::StdExpansion2DSharedPtr locExp =
//                             boost::static_pointer_cast<StdRegions::StdExpansion2D>(locExpList->GetExp(j));
// 
//                         nummodes = locExp->GetNcoeffs();
// 
//                         // Create a new BasisKey for projecting the dirichlet boundary conditions onto the boundary.
//                         // If the original BasisKey has N modes, the new BasisKey employs N GLL quadrature points such
//                         // that the FwdTrans using this new basis is in fact a collocation projection trough this
//                         // GLL-points. As a result, the expansion will be C0 continuous on the boundary.
// 
//                         // The PointsKey used for the (collocation) projection
//                         LibUtilities::PointsKey collPointsKey_0(nummodes,LibUtilities::eGaussLobattoLegendre);
//                         LibUtilities::PointsKey collPointsKey_1(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0);
//                         // The BasisKey used for the (collocation) projection
//                         LibUtilities::BasisKey collBasisKey_0(locExp->GetBasisType(0),nummodes,collPointsKey_0);
//                         LibUtilities::BasisKey collBasisKey_1(locExp->GetBasisType(1),nummodes,collPointsKey_1);
// 
//                         // Create a segment based on the new BasisKey in order to perfrom the projection
//                         //segGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(locSegExp->GetGeom());
//                         //collSeg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(collBasisKey, segGeom);
//                         //collStdSeg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(collBasisKey);
// 
//                         // Create a segment based on the new BasisKey in order to perfrom the projection
//                         SpatialDomains::Geometry2DSharedPtr geom = boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>(locExp->GetGeom());
//                         StdRegions::StdQuadExpSharedPtr collStd2D = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(collBasisKey_0, collBasisKey_1);
//                         LocalRegions::QuadExpSharedPtr coll2D = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(collBasisKey_0, collBasisKey_1, geom);
// 
//                         // Calculate the coordinates of the N GLL quadrature points
//                         Array<OneD,NekDouble> x0(nummodes, 0.0);
//                         Array<OneD,NekDouble> x1(nummodes, 0.0);
//                         Array<OneD,NekDouble> x2(nummodes, 0.0);
//                         coll2D->GetCoords(x0,x1,x2);
// 
//                         // Evaluate the Dirichlet boundary condition at the N GLL quadrature points
//                         for(m = 0; m < nummodes; ++m)
//                         {
//                             (coll2D->UpdatePhys())[m] = boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(locBCond)->
//                              m_DirichletCondition.Evaluate(x0[m],x1[m],x2[m]);
//                         }
//                         // Perform a FwdTrans() to calculate the expansion coefficients.
//                         // As both the original as the new expansion (used for the projection) are of the same order,
//                         // the coefficients will be identical and hence, can directly be stored at the right place
//                         Array<OneD,NekDouble> outarray;
//                         collStd2D->FwdTrans(coll2D->GetPhys(),outarray = (locExpList->UpdateCoeffs()) + coeffcnt);
// 
//                         coeffcnt += nummodes;
//                     }
// 
//                     m_bndConstraint[cnt] = locExpList;
//                     m_bndTypes[cnt++] = SpatialDomains::eDirichlet;
//                 } // end if Dirichlet
//             }
//             // list other boundaries
//             for(i = 0; i < nbnd; ++i)
//             {
//                 locBCond = (*(bconditions[i]))[variable];
// 
//                 if(locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
//                 {
//                     locBregion = bregions[i];
// 
//                     //TODO: check this
//                     locExpList = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(*locBregion,graph3D);
// 
//                     npoints = locExpList->GetPointsTot();
//                     Array<OneD,NekDouble> x0(npoints,0.0);
//                     Array<OneD,NekDouble> x1(npoints,0.0);
//                     Array<OneD,NekDouble> x2(npoints,0.0);
//                     locExpList->GetCoords(x0,x1,x2);
// 
//                     if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
//                     {
//                         for(m = 0; m < npoints; m++)
//                         {
//                             (locExpList->UpdatePhys())[m] = boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(locBCond)->
//                                 m_NeumannCondition.Evaluate(x0[m],x1[m],x2[m]);
//                         }
//                         locExpList->IProductWRTBase(*locExpList);
//                         m_bndConstraint[cnt] = locExpList;
//                         m_bndTypes[cnt++] = SpatialDomains::eNeumann;
//                     }
//                     else if(locBCond->GetBoundaryConditionType() == SpatialDomains::eRobin)
//                     {
//                         boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC =
//                             boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(locBCond);
//                         for(m = 0; m < npoints; m++)
//                         {
//                             (locExpList->UpdatePhys())[m] = -robinBC->m_a.Evaluate(x0[m],x1[m],x2[m])/
//                                 robinBC->m_b.Evaluate(x0[m],x1[m],x2[m]);
//                         }
//                         locExpList->IProductWRTBase(*locExpList);
//                         m_bndConstraint[cnt] = locExpList;
//                         m_bndTypes[cnt++] = SpatialDomains::eRobin;
//                     }
//                     else
//                     {
//                         ASSERTL0(false,"This type of BC not implemented yet");
//                     }
//                 }
//             }
//         }

        void ContField3D::FwdTrans(const ExpList &In)
        {
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,In);

            m_transState = eContinuous;
            m_physState = false;
        }

       // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
        void ContField3D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            // Note -1.0 term necessary to invert forcing function to
            // be consistent with matrix definition
            GlobalSolve(key,In,-1.0);
        }

       // TODO: check
       void ContField3D::GlobalSolve(const GlobalLinSysKey &key, const ExpList &Rhs, NekDouble ScaleForcing)
        {
            int i,j;
            int bndcnt=0;
            int NumDirBcs = m_locToGloMap->GetNumDirichletDofs();
            Array<OneD,NekDouble> sln;
            Array<OneD,NekDouble> init(m_contNcoeffs,0.0);
            Array<OneD,NekDouble> Dir_fce(m_contNcoeffs,0.0);

            //assume m_contCoeffs contains initial estimate
            // Set BCs in m_contCoeffs
            Blas::Dcopy(m_contNcoeffs, m_contCoeffs, 1, init, 1);
            for(i = 0; i < m_bndConstraint.num_elements(); ++i)
            {
                if(m_bndTypes[i] != SpatialDomains::eDirichlet)
                {
                    break;
                }
                for(j = 0; j < (m_bndConstraint[i])->GetNcoeffs(); j++)
                {
                    init[m_locToGloMap->GetBndCondMap(bndcnt++)] = (m_bndConstraint[i]->GetCoeffs())[j];
                }
            }
            GeneralMatrixOp(key, init, Dir_fce);

            // Set up forcing function
            IProductWRTBase(Rhs);
            
            // apply scaling of forcing term; (typically used to negate helmholtz forcing);
            Vmath::Smul(m_contNcoeffs, ScaleForcing, m_contCoeffs, 1, m_contCoeffs, 1);

            // Forcing function with Dirichlet conditions 
            Vmath::Vsub(m_contNcoeffs, m_contCoeffs, 1, Dir_fce, 1, m_contCoeffs, 1);

            // Forcing function with weak boundary conditions 
            for(i; i < m_bndConstraint.num_elements(); ++i)
            {
                for(j = 0; j < (m_bndConstraint[i])->GetNcoeffs(); j++)
                {
                    m_contCoeffs[m_locToGloMap->GetBndCondMap(bndcnt++)] +=  
                        (m_bndConstraint[i]->GetCoeffs())[j];
                }
            }

            if(m_contNcoeffs - NumDirBcs > 0)
            {

                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);
                       
                sln = m_contCoeffs+NumDirBcs;
                LinSys->Solve(sln,sln,*m_locToGloMap);
            }

            // Recover solution by addinig intial conditons
            Vmath::Zero(NumDirBcs, m_contCoeffs, 1);
            Vmath::Vadd(m_contNcoeffs, init, 1, m_contCoeffs, 1, m_contCoeffs, 1);

            m_transState = eContinuous;
            m_physState = false;
        }


        GlobalLinSysSharedPtr ContField3D::GetGlobalLinSys(const GlobalLinSysKey &mkey)
        {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(mkey);

            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalLinSys(mkey,m_locToGloMap);
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }


    

  } //end of namespace
} //end of namespace
