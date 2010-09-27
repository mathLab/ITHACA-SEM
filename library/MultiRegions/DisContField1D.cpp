///////////////////////////////////////////////////////////////////////////////
//
// File DisContField1D.cpp
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
// Description: Field definition for 1D domain with boundary
// conditions using LDG-H
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/DisContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class DisContField1D
         * This class augments the list of local expansions inherited from
         * ExpList1D with boundary conditions. Inter-element boundaries are
         * handled using an discontinuous Galerkin scheme.
         */

        /**
         * Constructs an empty expansion list with no boundary conditions.
         */
        DisContField1D::DisContField1D():
            ExpList1D(),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }


        /**
         * @todo Implement this constructor.
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   solnType    Type of global system to use.
         * @param   constructMap    ?
         */
        DisContField1D::DisContField1D(
                    SpatialDomains::MeshGraph1D &graph1D,
                    const GlobalSysSolnType solnType,
                    const bool constructMap):
            ExpList1D(graph1D),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }


        /**
         * Constructs a field as a copy of an existing field.
         * @param   In          Existing DisContField1D object to copy.
         */
        DisContField1D::DisContField1D(const DisContField1D &In):
            ExpList1D(In),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(In.m_bndCondExpansions),
            m_bndConditions(In.m_bndConditions),
            m_globalBndMat(In.m_globalBndMat),
            m_trace(In.m_trace),
            m_traceMap(In.m_traceMap)
        {
        }


        /**
         * An expansion list for the boundary expansions is generated first for
         * the field. These are subsequently evaluated for time zero. The trace
         * map is then constructed.
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   bc_loc      The index of the session variable associated
         *                      with the boundary conditions to enforce.
         * @param   solnType    Type of global system to use.
         */
        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const int bc_loc,
                    const GlobalSysSolnType solnType):
            ExpList1D(graph1D),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,
                                               bcs.GetVariable(bc_loc));
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,bcs.GetVariable(bc_loc),
                                periodicVertices);

            m_globalBndMat
                        = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            //GenerateFieldBnd1D(bcs,bcs.GetVariable(bc_loc));

            m_traceMap = MemoryManager<LocalToGlobalDGMap>
                                        ::AllocateSharedPtr(graph1D,*this,
                                                            solnType,
                                                            m_bndCondExpansions,
                                                            m_bndConditions);

            m_trace = Array<OneD,NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
        }


        /**
         * An expansion list for the boundary expansions is generated first for
         * the field. These are subsequently evaluated for time zero. The trace
         * map is then constructed.
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce.
         * @param   solnType    Type of global system to use.
         */
        DisContField1D::DisContField1D(SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable,
                    const GlobalSysSolnType solnType):
            ExpList1D(graph1D),
            m_numDirBndCondExpansions(0),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            GenerateBoundaryConditionExpansion(graph1D,bcs,variable);
            EvaluateBoundaryConditions();
            ApplyGeomInfo(graph1D);

            map<int,int> periodicVertices;
            GetPeriodicVertices(graph1D,bcs,variable,periodicVertices);

            m_globalBndMat
                        = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

            //GenerateFieldBnd1D(bcs,variable);

            m_traceMap = MemoryManager<LocalToGlobalDGMap>::
                AllocateSharedPtr(graph1D,*this,solnType,
                                  m_bndCondExpansions,m_bndConditions);

            m_trace = Array<OneD,NekDouble>(m_traceMap->GetNumLocalBndCoeffs());
        }


        /**
         * Generate the boundary condition expansion list
         * @param   graph1D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansions.
         * @param   bcs         Information about the enforced boundary
         *                      conditions.
         * @param   variable    The session variable associated with the
         *                      boundary conditions to enforce.
         */
        void DisContField1D::GenerateBoundaryConditionExpansion(
                    const SpatialDomains::MeshGraph1D &graph1D,
                    SpatialDomains::BoundaryConditions &bcs,
                    const std::string variable)
        {
            int i,j;
            int cnt  = 0;
            int cnt2 = 0;

            SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();
            // count the number of non-periodic boundary points
            for(i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType()
                                                != SpatialDomains::ePeriodic )
                {
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        cnt += (*bregions[i])[j]->size();
                    }

                    if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() == SpatialDomains::eDirichlet )
                    {
                        for(j = 0; j < bregions[i]->size(); j++)
                        {
                            cnt2 += (*bregions[i])[j]->size();
                        }
                    }
                }
            }

            m_numDirBndCondExpansions = cnt2;

            m_bndCondExpansions
                    = Array<OneD,LocalRegions::PointExpSharedPtr>(cnt);

            m_bndConditions
                    = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);

            SetBoundaryConditionExpansion(graph1D,bcs,variable,
                                           m_bndCondExpansions,
                                           m_bndConditions);
        }
#if 0
        void DisContField1D::GenerateFieldBnd1D(SpatialDomains::BoundaryConditions &bcs,   const std::string variable)
        {
            int i,nbnd;
            int cnt=0;
            LocalRegions::PointExpSharedPtr  p_exp;
            SpatialDomains::BoundaryConditionShPtr locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;
            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            nbnd = bregions.size();

            // count the number of non-periodic boundary regions
            for(int i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))[variable])->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    for(j = 0; j < bregions[i]->size(); ++j)
                    {
                        cnt += (*bregions[i])[j]->size();
                    }
                }
            }

            m_bndCondExpansions = Array<OneD,LocalRegions::PointExpSharedPtr>(cnt);
            m_bndConditions     = Array<OneD,SpatialDomains::BoundaryConditionShPtr>(cnt);

            // Set up matrix map
            m_globalBndMat   = MemoryManager<GlobalLinSysMap>::AllocateSharedPtr();

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
                                Array<OneD,NekDouble> coords(3,0.0);

                                p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                                vert->GetCoords(coords);
                                p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(locBCond)->m_DirichletCondition.Evaluate(coords[0],coords[1],coords[2]));

                                m_bndCondExpansions[cnt] = p_exp;
                                m_bndConditions[cnt++] = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,"dynamic cast to a vertex failed");
                            }
                }
            }

            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];

                if(locBCond->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        Array<OneD,NekDouble> coords(3,0.0);

                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        vert->GetCoords(coords);

                        if(locBCond->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                        {
                            p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(locBCond)
                                            ->m_NeumannCondition.Evaluate(coords[0],coords[1],coords[2]));
                            m_bndCondExpansions[cnt] = p_exp;
                            m_bndConditions[cnt++] = locBCond;
                        }
                        else if(locBCond->GetBoundaryConditionType() == SpatialDomains::eRobin)
                        {
                            boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC  =
                                boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(locBCond);
                            p_exp->SetValue(-robinBC->m_a.Evaluate(coords[0],coords[1],coords[2])/
                                            robinBC->m_b.Evaluate(coords[0],coords[1],coords[2]));
                            m_bndCondExpansions[cnt] = p_exp;
                            m_bndConditions[cnt++] = locBCond;
                        }
                        else
                        {
                            ASSERTL0(false,"This type of BC not implemented yet");
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex failed");
                    }
                }
            }
        }
#endif


        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   periodicVertices    Map into which the list of periodic
         *                      vertices is placed.
         */
        void DisContField1D::GetPeriodicVertices(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                      map<int,int>& periodicVertices)
        {

            int i,j,k;

            SpatialDomains::BoundaryRegionCollection &bregions
                    = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                    = bcs.GetBoundaryConditions();

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
                if(locBCond->GetBoundaryConditionType()
                        == SpatialDomains::ePeriodic)
                {
                    region1ID = i;
                    region2ID = (boost::static_pointer_cast<SpatialDomains::
                                    PeriodicBoundaryCondition>(locBCond))
                                        ->m_connectedBoundaryRegion;

                    if(doneBndRegions.count(region1ID)==0)
                    {
                        ASSERTL0(bregions[region1ID]->size()
                                    == bregions[region2ID]->size(),
                                 "Size of the 2 periodic boundary regions "
                                 "should be equal");

                        for(j = 0; j < bregions[region1ID]->size(); j++)
                        {
                            comp1 = (*(bregions[region1ID]))[j];
                            comp2 = (*(bregions[region2ID]))[j];

                            ASSERTL0(comp1->size() == comp2->size(),
                                     "Size of the 2 periodic composites should "
                                     "be equal");

                            for(k = 0; k < comp1->size(); k++)
                            {
                                if(!(vert1 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp1)[k]))||
                                   !(vert2 = boost::dynamic_pointer_cast
                                        <SpatialDomains::VertexComponent>(
                                            (*comp2)[k])))
                                {
                                    ASSERTL0(false,"dynamic cast to a "
                                                   "VertexComponent failed");
                                }

                                // Extract the periodic vertices
                                periodicVertices[vert1->GetVid()]
                                    = vert2->GetVid();
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

        /**
         * @param   graph1D     A mesh containing information about the domain
         *                      and the Spectral/hp element expansion.
         * @param   bcs         Information about the boundary conditions.
         * @param   variable    Specifies the field.
         * @param   bndCondExpansions   Array of ExpList1D objects each
         *                      containing a 1D spectral/hp element expansion
         *                      on a single boundary region.
         * @param   bncConditions   Array of BoundaryCondition objects which
         *                      contain information about the boundary
         *                      conditions on the different boundary regions.
         */
        void DisContField1D::SetBoundaryConditionExpansion(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                Array<OneD, LocalRegions::PointExpSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions)
        {
            int i,j,k;
            int cnt  = 0;

            SpatialDomains::BoundaryRegionCollection &bregions
                                                = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions
                                                = bcs.GetBoundaryConditions();

            LocalRegions::PointExpSharedPtr          locPointExp;
            SpatialDomains::BoundaryConditionShPtr   locBCond;
            SpatialDomains::VertexComponentSharedPtr vert;

            int nbnd = bregions.size();

            cnt=0;
            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];
                if(locBCond->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bregions[i]->size(); j++)
                    {
                        for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                        {
                            if(vert = boost::dynamic_pointer_cast
                                    <SpatialDomains::VertexComponent>(
                                        (*(*bregions[i])[j])[k]))
                            {
                                locPointExp
                                    = MemoryManager<LocalRegions::PointExp>
                                                ::AllocateSharedPtr(vert);
                                bndCondExpansions[cnt]  = locPointExp;
                                bndConditions[cnt++]    = locBCond;
                            }
                            else
                            {
                                ASSERTL0(false,
                                         "dynamic cast to a vertex failed");
                            }
                        }
                    }
                }
            } // end if Dirichlet

            // then, list the other (non-periodic) boundaries
            for(i = 0; i < nbnd; ++i)
            {
                locBCond = (*(bconditions[i]))[variable];

                switch(locBCond->GetBoundaryConditionType())
                {
                case SpatialDomains::eNeumann:
                case SpatialDomains::eRobin:
                    {
                        for(j = 0; j < bregions[i]->size(); j++)
                        {
                            for(k = 0; k < ((*bregions[i])[j])->size(); k++)
                            {
                                if(vert = boost::dynamic_pointer_cast
                                   <SpatialDomains::VertexComponent>(
                                        (*(*bregions[i])[j])[k]))
                                {
                                    locPointExp
                                        = MemoryManager<LocalRegions::PointExp>
                                        ::AllocateSharedPtr(vert);
                                    bndCondExpansions[cnt]  = locPointExp;
                                    bndConditions[cnt++]    = locBCond;
                                }
                                else
                                {
                                    ASSERTL0(false,
                                             "dynamic cast to a vertex failed");
                                }
                            }
                        }
                    }
                case SpatialDomains::eDirichlet: // do nothing for these types
                case SpatialDomains::ePeriodic:
                    break;
                default:
                    ASSERTL0(false,"This type of BC not implemented yet");
                    break;
                }
            }
        }


        /**
         *
         */
        DisContField1D::~DisContField1D()
        {
        }


        /**
         *
         */
        GlobalLinSysSharedPtr DisContField1D::GetGlobalBndLinSys(
                    const GlobalLinSysKey &mkey)
        {
            ASSERTL0(mkey.GetMatrixType() == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");

            ASSERTL1(mkey.GetGlobalSysSolnType() != eDirectFullMatrix,
                     "Full matrix global systems are not supported for HDG "
                     "expansions");

            ASSERTL1(mkey.GetGlobalSysSolnType()
                                        ==m_traceMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");

            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalBndMat->find(mkey);

            if(matrixIter == m_globalBndMat->end())
            {
                glo_matrix = GenGlobalBndLinSys(mkey,m_traceMap);
                (*m_globalBndMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }


        /**
         *
         */
        void DisContField1D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveDG(inarray, outarray, lambda, varLambda, varCoeff, 1);
        }


        /**
         *
         */
        void DisContField1D::v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                          NekDouble tau)
        {
            int i,n,cnt,nbndry;
            int nexp = GetExpSize();
            Array<OneD,NekDouble> f(m_ncoeffs);
            DNekVec F(m_ncoeffs,f,eWrapper);
            Array<OneD,NekDouble> e_f, e_l;

            //----------------------------------
            // Setup RHS Inner product
            //----------------------------------
            IProductWRTBase(inarray,f);
            Vmath::Neg(m_ncoeffs,f,1);

            //----------------------------------
            // Solve continuous Boundary System
            //----------------------------------
            int GloBndDofs = m_traceMap->GetNumGlobalBndCoeffs();
            int NumDirBCs  = m_traceMap->GetNumLocalDirBndCoeffs();
            int e_ncoeffs,id;

            GlobalMatrixKey HDGLamToUKey(StdRegions::eHybridDGLamToU,
                                         lambda,tau);

            const DNekScalBlkMatSharedPtr &HDGLamToU
                                                = GetBlockMatrix(HDGLamToUKey);

            Array<OneD,NekDouble> BndRhs(GloBndDofs,0.0);
            // Zero trace space
            Vmath::Zero(GloBndDofs,m_trace,1);

            int     LocBndCoeffs = m_traceMap->GetNumLocalBndCoeffs();
            Array<OneD, NekDouble> loc_lambda(LocBndCoeffs);
            DNekVec LocLambda(LocBndCoeffs,loc_lambda,eWrapper);

            //----------------------------------
            // Evaluate Trace Forcing
            //----------------------------------
            // Determing <u_lam,f> terms using HDGLamToU matrix
            for(cnt = n = 0; n < nexp; ++n)
            {
                nbndry = (*m_exp)[n]->NumDGBndryCoeffs();

                e_ncoeffs = (*m_exp)[n]->GetNcoeffs();
                e_f       = f+m_coeff_offset[n];
                e_l       = loc_lambda + cnt;

                // use outarray as tmp space
                DNekVec     Floc    (nbndry, e_l, eWrapper);
                DNekVec     ElmtFce (e_ncoeffs, e_f, eWrapper);
                Floc = Transpose(*(HDGLamToU->GetBlock(n,n)))*ElmtFce;

                cnt += nbndry;
            }

            // Assemble into global operator
            m_traceMap->AssembleBnd(loc_lambda,BndRhs);

            cnt = 0;
            // Copy Dirichlet boundary conditions into trace space
            for(i = 0; i < m_numDirBndCondExpansions; ++i)
            {
                id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                m_trace[id] = m_bndCondExpansions[i]->GetCoeff(0);
            }

            //Add weak boundary condition to trace forcing
            for(i = m_numDirBndCondExpansions;
                                    i < m_bndCondExpansions.num_elements(); ++i)
            {
                id = m_traceMap->GetBndCondCoeffsToGlobalCoeffsMap(i);
                BndRhs[id] += m_bndCondExpansions[i]->GetCoeff(0);
            }

            //----------------------------------
            // Solve trace problem
            //----------------------------------
            if(GloBndDofs - NumDirBCs > 0)
            {
                GlobalLinSysKey       key(StdRegions::eHybridDGHelmBndLam,
                                          m_traceMap,lambda,tau);
                GlobalLinSysSharedPtr LinSys = GetGlobalBndLinSys(key);
                LinSys->Solve(BndRhs,m_trace,m_traceMap);
            }

            //----------------------------------
            // Internal element solves
            //----------------------------------
            GlobalMatrixKey invHDGhelmkey(StdRegions::eInvHybridDGHelmholtz,
                                          lambda,tau);

            const DNekScalBlkMatSharedPtr& InvHDGHelm
                                                = GetBlockMatrix(invHDGhelmkey);
            DNekVec out(m_ncoeffs,outarray,eWrapper);
            Vmath::Zero(m_ncoeffs,outarray,1);

            // get local trace solution from BndSol
            m_traceMap->GlobalToLocalBnd(m_trace,loc_lambda);

            //  out =  u_f + u_lam = (*InvHDGHelm)*f + (LamtoU)*Lam
            out = (*InvHDGHelm)*F + (*HDGLamToU)*LocLambda;
        }


        /**
         *
         */
        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>&
                                            DisContField1D::v_GetBndConditions()
        {
            return m_bndConditions;
        }


        /**
         * Based on the expression \f$g(x,t)\f$ for the boundary conditions,
         * this function evaluates the boundary conditions for all boundaries
         * at time-level \a t.
         * @param   time        The time at which the boundary conditions
         *                      should be evaluated.
         * @param   bndCondExpansions   List of boundary expansions.
         * @param   bndConditions   Information about the boundary conditions.
         */
        void DisContField1D::v_EvaluateBoundaryConditions(
                                const NekDouble time, const NekDouble x2_in)
        {
            int i;

            NekDouble x0;
            NekDouble x1;
            NekDouble x2;

            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                m_bndCondExpansions[i]->GetCoords(x0,x1,x2);

                if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eDirichlet)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::DirichletBoundaryCondition>(m_bndConditions[i])
                             ->m_dirichletCondition).Evaluate(x0,x1,x2,time));
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eNeumann)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::NeumannBoundaryCondition>(m_bndConditions[i])
                             ->m_neumannCondition).Evaluate(x0,x1,x2,time));
                }
                else if(m_bndConditions[i]->GetBoundaryConditionType()
                        == SpatialDomains::eRobin)
                {
                    m_bndCondExpansions[i]->SetCoeff(
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinFunction).Evaluate(x0,x1,x2,time));

                    m_bndCondExpansions[i]->SetPhys(
                            (boost::static_pointer_cast<SpatialDomains
                             ::RobinBoundaryCondition>(m_bndConditions[i])
                             ->m_robinPrimitiveCoeff).Evaluate(x0,x1,x2,time));

                }
                else
                {
                    ASSERTL0(false,"This type of BC not implemented yet");
                }
            }
        }



        // Set up a list of element ids and edge ids that link to the
        // boundary conditions
        void DisContField1D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &VertID)
        {
            map<int, int> VertGID;
            int i,n,id;
            int bid,cnt,Vid;
            int nbcs = m_bndConditions.num_elements();

            // make sure arrays are of sufficient length
            if(ElmtID.num_elements() != nbcs)
            {
                ElmtID = Array<OneD, int>(nbcs,-1);
            }
            else
            {
                fill(ElmtID.get(), ElmtID.get()+nbcs, -1);
            }

            if(VertID.num_elements() != nbcs)
            {
                VertID = Array<OneD, int>(nbcs);
            }

            // setup map of all global ids along boundary
            for(cnt = n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                Vid =  m_bndCondExpansions[n]->GetGeom()->GetVid();
                VertGID[Vid] = cnt++;
            }

            // Loop over elements and find verts that match;
            for(cnt = n = 0; n < GetExpSize(); ++n)
            {
                for(i = 0; i < (*m_exp)[n]->GetNverts(); ++i)
                {
                    id = (*m_exp)[n]->GetGeom()->GetVid(i);

                    if(VertGID.count(id) > 0)
                    {
                        bid = VertGID.find(id)->second;
                        ASSERTL1(ElmtID[bid] == -1,"Edge already set");
                        ElmtID[bid] = n;
                        VertID[bid] = i;
                        cnt ++;
                    }
                }
            }

            ASSERTL1(cnt == nbcs,"Failed to visit all boundary condtiions");
        }

        /**
         * Search through the edge expansions and identify which ones
         * have Robin/Mixed type boundary conditions. If find a Robin
         * boundary then store the edge id of the boundary condition
         * and the array of points of the physical space boundary
         * condition which are hold the boundary condition primitive
         * variable coefficient at the quatrature points
         *
         * \return std map containing the robin boundary condition
         * info using a key of the element id
         *
         * There is a next member to allow for more than one Robin
         * boundary condition per element
         */
        map<int, RobinBCInfoSharedPtr> DisContField1D::GetRobinBCInfo(void)
        {
            int i;
            map<int, RobinBCInfoSharedPtr> returnval;
            Array<OneD, int> ElmtID,VertID;
            GetBoundaryToElmtMap(ElmtID,VertID);

            for(i = 0; i < m_bndCondExpansions.num_elements(); ++i)
            {
                LocalRegions::PointExpSharedPtr locExpList;

                if(m_bndConditions[i]->GetBoundaryConditionType() == SpatialDomains::eRobin)
                {
                    int elmtid;
                    Array<OneD, NekDouble> Array_tmp;

                    locExpList = m_bndCondExpansions[i];

                    RobinBCInfoSharedPtr rInfo = MemoryManager<RobinBCInfo>::AllocateSharedPtr(VertID[i],Array_tmp = locExpList->GetPhys());

                    elmtid = ElmtID[i];
                    // make link list if necessary (not likely in
                    // 1D but needed in 2D & 3D)
                    if(returnval.count(elmtid) != 0)
                    {
                        rInfo->next = returnval.find(elmtid)->second;
                    }
                    returnval[elmtid] = rInfo;
                }
            }

            return returnval;
        }

    } // end of namespace
} //end of namespace
