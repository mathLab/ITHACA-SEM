///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapDG.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Local to Global Base Class mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <MultiRegions/AssemblyMap/AssemblyCommDG.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/PointExp.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

#include <boost/core/ignore_unused.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        AssemblyMapDG::AssemblyMapDG():
            m_numDirichletBndPhys(0)
        {
        }

        AssemblyMapDG::~AssemblyMapDG()
        {
        }

        AssemblyMapDG::AssemblyMapDG(
            const LibUtilities::SessionReaderSharedPtr                &pSession,
            const SpatialDomains::MeshGraphSharedPtr                  &graph,
            const ExpListSharedPtr                                    &trace,
            const ExpList                                             &locExp,
            const Array<OneD, const MultiRegions::ExpListSharedPtr>         &bndCondExp,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond,
            const PeriodicMap                                         &periodicTrace,
            const std::string variable):
            AssemblyMap(pSession,variable)
        {
            boost::ignore_unused(graph);

            int i, j, k, cnt, id, id1, gid;
            int order_e   = 0;
            int nTraceExp = trace->GetExpSize();
            int nbnd      = bndCondExp.size();

            LocalRegions::ExpansionSharedPtr  exp;
            LocalRegions::ExpansionSharedPtr  bndExp;
            SpatialDomains::GeometrySharedPtr traceGeom;

            const LocalRegions::ExpansionVector expList = *(locExp.GetExp());
            int nel = expList.size();

            map<int, int> meshTraceId;

            m_signChange = true;

            // determine mapping from geometry edges to trace
            for(i = 0; i < nTraceExp; ++i)
            {
                meshTraceId[trace->GetExp(i)->GetGeom()->GetGlobalID()] = i;
            }

            // Count total number of trace elements
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += expList[i]->GetNtrace();
            }

            Array<OneD, LocalRegions::ExpansionSharedPtr> tracemap(cnt);
            m_elmtToTrace = Array<
                OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >(nel);

            // set up trace expansions links;
            for(cnt = i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = tracemap + cnt;
                exp = expList[i];
                
                for(j = 0; j < exp->GetNtrace(); ++j)
                {
                    id = exp->GetGeom()->GetTid(j);

                    if(meshTraceId.count(id) > 0)
                    {
                        m_elmtToTrace[i][j] =
                            trace->GetExp(meshTraceId.find(id)->second);
                    }
                    else
                    {
                        ASSERTL0(false, "Failed to find trace map");
                    }
                }

                cnt += exp->GetNtrace();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                        SpatialDomains::ePeriodic)
                {
                    continue;
                }
                cnt += bndCondExp[i]->GetExpSize();
            }

            set<int> dirTrace;

            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys  = 0;

            for(i = 0; i < bndCondExp.size(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndExp    = bndCondExp[i]->GetExp(j);
                    traceGeom = bndExp->GetGeom();
                    id        = traceGeom->GetGlobalID();

                    if(bndCond[i]->GetBoundaryConditionType() ==
                           SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs += bndExp->GetNcoeffs();
                        m_numDirichletBndPhys  += bndExp->GetTotPoints();
                        dirTrace.insert(id);
                    }
                }
            }

            // Set up integer mapping array and sign change for each degree of
            // freedom + initialise some more data members.
            m_staticCondLevel           = 0;
            m_lowestStaticCondLevel     = 0;
            m_numPatches                = nel;
            m_numLocalBndCoeffsPerPatch = Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch = Array<OneD, unsigned int>(nel);

            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                int BndCoeffs = expList[i]->NumDGBndryCoeffs();
                nbndry += BndCoeffs;
                m_numLocalIntCoeffsPerPatch[i] = 0;
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) BndCoeffs;
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            m_numLocalBndCoeffs     = nbndry;
            m_numLocalCoeffs        = nbndry;
            m_localToGlobalBndMap   = Array<OneD, int>       (nbndry);
            m_localToGlobalBndSign  = Array<OneD, NekDouble> (nbndry,1);

            // Set up array for potential mesh optimsation
            Array<OneD,int> traceElmtGid(nTraceExp, -1);
            int nDir = 0;
            cnt = 0;

            // We are now going to construct a graph of the mesh which can be
            // reordered depending on the type of solver we would like to use.
            typedef boost::adjacency_list<
                boost::setS, boost::vecS, boost::undirectedS> BoostGraph;

            BoostGraph boostGraphObj;
            int trace_id, trace_id1;
            int dirOffset = 0;

            // make trace renumbering map where first solved trace starts
            // at 0 so we can set up graph.
            for(i = 0; i < nTraceExp; ++i)
            {
                id = trace->GetExp(i)->GetGeom()->GetGlobalID();

                if (dirTrace.count(id) == 0)
                {
                    // Initial put in element ordering (starting from zero) into
                    // traceElmtGid
                    boost::add_vertex(boostGraphObj);
                    traceElmtGid[i] = cnt++;
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    traceElmtGid[i] = dirOffset;
                    dirOffset      += trace->GetExp(i)->GetNcoeffs();
                    nDir++;
                }
            }

            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                exp = expList[i];
                
                for(j = 0; j < exp->GetNtrace(); ++j)
                {
                    // Add trace to boost graph for non-Dirichlet Boundary
                    traceGeom = m_elmtToTrace[i][j]->GetGeom();
                    id        = traceGeom->GetGlobalID();
                    trace_id  = meshTraceId.find(id)->second;

                    if(dirTrace.count(id) == 0)
                    {
                        for(k = j+1; k < exp->GetNtrace(); ++k)
                        {
                            traceGeom = m_elmtToTrace[i][k]->GetGeom();
                            id1       = traceGeom->GetGlobalID();
                            trace_id1 = meshTraceId.find(id1)->second;

                            if(dirTrace.count(id1) == 0)
                            {
                                boost::add_edge((size_t)traceElmtGid[trace_id],
                                                (size_t)traceElmtGid[trace_id1],
                                                boostGraphObj);
                            }
                        }
                    }
                }
            }

            int                                 nGraphVerts = nTraceExp - nDir;
            Array<OneD, int>                    perm (nGraphVerts);
            Array<OneD, int>                    iperm(nGraphVerts);
            Array<OneD, int>                    vwgts(nGraphVerts);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;

            for(i = 0; i < nGraphVerts; ++i)
            {
                vwgts[i] = trace->GetExp(i+nDir)->GetNcoeffs();
            }

            if(nGraphVerts)
            {
                switch(m_solnType)
                {
                    case eDirectFullMatrix:
                    case eIterativeFull:
                    case eIterativeStaticCond:
                    case eXxtFullMatrix:
                    case eXxtStaticCond:
                    case ePETScFullMatrix:
                    case ePETScStaticCond:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                        break;
                    }
                    case eDirectMultiLevelStaticCond:
                    case eIterativeMultiLevelStaticCond:
                    case eXxtMultiLevelStaticCond:
                    case ePETScMultiLevelStaticCond:
                    {
                        MultiLevelBisectionReordering(boostGraphObj,perm,iperm,
                                                      bottomUpGraph);
                        break;
                    }
                    default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }
            }

            // Recast the permutation so that it can be used as a map from old
            // trace ID to new trace ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < nTraceExp - nDir; ++i)
            {
                traceElmtGid[perm[i]+nDir] = cnt;
                cnt += trace->GetExp(perm[i]+nDir)->GetNcoeffs();
            }

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                exp = expList[i];

                for(j = 0; j < exp->GetNtrace(); ++j)
                {
                    traceGeom = m_elmtToTrace[i][j]->GetGeom();
                    id        = traceGeom->GetGlobalID();
                    gid       = traceElmtGid[meshTraceId.find(id)->second];

                    const int nDim = exp->GetNumBases();

                    if (nDim == 1)
                    {
                        order_e = 1;
                        m_localToGlobalBndMap[cnt] = gid;
                    }
                    else if (nDim == 2)
                    {
                        order_e = exp->GetEdgeNcoeffs(j);
                    
                        if(exp->GetEorient(j) == StdRegions::eForwards)
                        {
                            for(k = 0; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[k+cnt] = gid + k;
                            }
                        }
                        else
                        {
                            switch(m_elmtToTrace[i][j]->GetBasisType(0))
                            {
                                case LibUtilities::eModified_A:
                                {
                                    // reverse vertex order
                                    m_localToGlobalBndMap[cnt]   = gid + 1;
                                    m_localToGlobalBndMap[cnt+1] = gid;
                                    for (k = 2; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[k+cnt] = gid + k;
                                    }

                                    // negate odd modes
                                    for(k = 3; k < order_e; k+=2)
                                    {
                                        m_localToGlobalBndSign[cnt+k] = -1.0;
                                    }
                                    break;
                                }
                                case LibUtilities::eGLL_Lagrange:
                                {
                                    // reverse  order
                                    for(k = 0; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                    }
                                    break;
                                }
                                case LibUtilities::eGauss_Lagrange:
                                {
                                    // reverse  order
                                    for(k = 0; k < order_e; ++k)
                                    {
                                        m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                                    }
                                    break;
                                }
                                default:
                                {
                                    ASSERTL0(false,"Boundary type not permitted");
                                }
                            }
                        }
                    }
                    else if (nDim == 3)
                    {
                        order_e = exp->GetFaceNcoeffs(j);

                        std::map<int, int> orientMap;

                        Array<OneD, unsigned int> elmMap1 (order_e);
                        Array<OneD,          int> elmSign1(order_e);
                        Array<OneD, unsigned int> elmMap2 (order_e);
                        Array<OneD,          int> elmSign2(order_e);

                        StdRegions::Orientation fo = exp->GetForient(j);

                        // Construct mapping which will permute global IDs
                        // according to face orientations.
                        exp->GetFaceToElementMap(j,fo,elmMap1,elmSign1);
                        exp->GetFaceToElementMap(
                            j,StdRegions::eDir1FwdDir1_Dir2FwdDir2,elmMap2,elmSign2);

                        for (k = 0; k < elmMap1.size(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < elmMap2.size(); ++l)
                            {
                                if (elmMap1[k] == elmMap2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            m_localToGlobalBndMap [k+cnt] = gid + orientMap[k];
                            m_localToGlobalBndSign[k+cnt] = elmSign2[orientMap[k]];
                        }
                    }

                    cnt += order_e;
                }
            }

            // set up identify map for lcoal to lcoal 
            m_localToLocalBndMap    = Array<OneD,int>(m_numLocalBndCoeffs);

            //local to bnd map is just a copy 
            for(i = 0; i < m_numLocalBndCoeffs; ++i)
            {
                m_localToLocalBndMap[i] = i;
            }


            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs    = m_numGlobalBndCoeffs;

            CalculateBndSystemBandWidth();

            // set up m_bndCondCoeffsToLocalTraceMap 
            // Number of boundary expansions
            int nbndexp = 0;
            int bndTotal = 0;
            int bndOffset;
            int traceOffset;

            cnt = 0; 
            for(i = 0; i < nbnd; ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }
                cnt += bndCondExp[i]->GetNcoeffs();
                nbndexp  += bndCondExp[i]->GetExpSize();
            }
            
            m_bndCondCoeffsToLocalTraceMap = Array<OneD, int>(cnt); 
            m_bndCondIDToGlobalTraceID = Array<OneD, int>(nbndexp);

            cnt = 0;
            for(i = 0; i < bndCondExp.size(); ++i)
            {
                if (bndCond[i]->GetBoundaryConditionType() ==
                    SpatialDomains::ePeriodic)
                {
                    continue;
                }

                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    bndExp    = bndCondExp[i]->GetExp(j);
                    id        = bndExp->GetGeom()->GetGlobalID();

                    int meshId = meshTraceId.find(id)->second;
                    m_bndCondIDToGlobalTraceID[cnt++] = meshId;

                    // initialy set up map with global bnd location
                    // and then use the localToGlobalBndMap to invert
                    // since this is a one to one mapping on boundaries
                    traceOffset = traceElmtGid[meshId];
                    bndOffset   = bndCondExp[i]->GetCoeff_Offset(j) + bndTotal;


                    for(k = 0; k < bndExp->GetNcoeffs(); ++k)
                    {
                        m_bndCondCoeffsToLocalTraceMap[bndOffset+k] =
                            traceOffset + k;
                    }
                }
                bndTotal += bndCondExp[i]->GetNcoeffs();
            }

            // generate an inverse local to global bnd map;
            map<int,int> invLocToGloMap;
            for(i = 0; i < nbndry; ++i)
            {
                invLocToGloMap[m_localToGlobalBndMap[i]] = i;
            }

            // reset bndCondCoeffToLocalTraceMap to hold local rather
            // than global reference
            for(i = 0; i < m_bndCondCoeffsToLocalTraceMap.size(); ++i)
            {
                m_bndCondCoeffsToLocalTraceMap[i] =
                    invLocToGloMap[m_bndCondCoeffsToLocalTraceMap[i]];
            }
            
            // Now set up mapping from global coefficients to universal.
            ExpListSharedPtr tr = std::dynamic_pointer_cast<ExpList>(trace);
            SetUpUniversalDGMap   (locExp);

            m_assemblyComm = AssemblyCommDGSharedPtr(
                MemoryManager<AssemblyCommDG>::AllocateSharedPtr(
                    locExp, tr, m_elmtToTrace, bndCondExp, bndCond, periodicTrace));

            if ((m_solnType == eDirectMultiLevelStaticCond ||
                 m_solnType == eIterativeMultiLevelStaticCond ||
                 m_solnType == eXxtMultiLevelStaticCond ||
                 m_solnType == ePETScMultiLevelStaticCond)
                && nGraphVerts)
            {
                if (m_staticCondLevel < (bottomUpGraph->GetNlevels() - 1))
                {
                    Array<OneD, int> vwgts_perm(nGraphVerts);

                    for (int i = 0; i < nGraphVerts; i++)
                    {
                        vwgts_perm[i] = vwgts[perm[i]];
                    }

                    bottomUpGraph->ExpandGraphWithVertexWeights(vwgts_perm);
                    m_nextLevelLocalToGlobalMap = MemoryManager<AssemblyMap>::
                        AllocateSharedPtr(this, bottomUpGraph);
                }
            }

            m_hash = hash_range(m_localToGlobalBndMap.begin(),
                                m_localToGlobalBndMap.end());
        }

        /**
         * Constructs a mapping between the process-local global numbering and
         * a universal numbering of the trace space expansion. The universal
         * numbering is defined by the mesh edge IDs to enforce consistency
         * across processes.
         *
         * @param       locExp  List of local elemental expansions.
         */
        void AssemblyMapDG::SetUpUniversalDGMap(const ExpList &locExp)
        {
            LocalRegions::ExpansionSharedPtr locExpansion;
            int cnt       = 0;
            int id        = 0;
            int order_e   = 0;
            int vGlobalId = 0;
            int maxDof    = 0;
            int dof       = 0;
            int nDim      = 0;
            int i,j,k;

            const LocalRegions::ExpansionVector &locExpVector = *(locExp.GetExp());

            // Initialise the global to universal maps.
            m_globalToUniversalBndMap = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);
            m_globalToUniversalBndMapUnique = Nektar::Array<OneD, int>(m_numGlobalBndCoeffs, -1);

            // Loop over all the elements in the domain and compute max
            // DOF. Reduce across all processes to get universal maximum.
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = locExpVector[i];
                nDim = locExpansion->GetShapeDimension();

                // Loop over all edges of element i
                if (nDim == 1)
                {
                    maxDof = (1 > maxDof ? 1 : maxDof);
                }
                else if (nDim == 2)
                {
                    for (j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        dof    = locExpansion->GetEdgeNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
                else if (nDim == 3)
                {
                    for (j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        dof    = locExpansion->GetFaceNcoeffs(j);
                        maxDof = (dof > maxDof ? dof : maxDof);
                    }
                }
            }
            m_comm->AllReduce(maxDof, LibUtilities::ReduceMax);

            // Now have trace edges Gid position
            cnt = 0;
            for(i = 0; i < locExpVector.size(); ++i)
            {
                locExpansion = locExpVector[i];
                nDim = locExpansion->GetShapeDimension();

                // Populate mapping for each edge of the element.
                if (nDim == 1)
                {
                    int nverts = locExpansion->GetNverts();
                    for(j = 0; j < nverts; ++j)
                    {
                        LocalRegions::PointExpSharedPtr locPointExp =
                            m_elmtToTrace[i][j]->as<LocalRegions::PointExp>();
                        id = locPointExp->GetGeom()->GetGlobalID();
                        vGlobalId = m_localToGlobalBndMap[cnt+j];
                        m_globalToUniversalBndMap[vGlobalId]
                            = id * maxDof + j + 1;
                    }
                    cnt += nverts;
                }
                else if (nDim == 2)
                {
                    for(j = 0; j < locExpansion->GetNedges(); ++j)
                    {
                        LocalRegions::SegExpSharedPtr locSegExp =
                            m_elmtToTrace[i][j]->as<LocalRegions::SegExp>();

                        id  = locSegExp->GetGeom()->GetGlobalID();
                        order_e = locExpansion->GetEdgeNcoeffs(j);

                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);

                        locExpansion->GetEdgeToElementMap(j, StdRegions::eForwards, map1, sign1);
                        locExpansion->GetEdgeToElementMap(j, locExpansion->GetEorient(j), map2, sign2);

                        for (k = 0; k < map1.size(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.size(); ++l)
                            {
                                if (map1[k] == map2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
                else if (nDim == 3)
                {
                    for(j = 0; j < locExpansion->GetNfaces(); ++j)
                    {
                        LocalRegions::Expansion2DSharedPtr locFaceExp =
                                m_elmtToTrace[i][j]
                                           ->as<LocalRegions::Expansion2D>();

                        id  = locFaceExp->GetGeom()->GetGlobalID();
                        order_e = locExpansion->GetFaceNcoeffs(j);

                        map<int,int> orientMap;
                        Array<OneD, unsigned int> map1(order_e), map2(order_e);
                        Array<OneD, int> sign1(order_e), sign2(order_e);

                        locExpansion->GetFaceToElementMap(j, StdRegions::eDir1FwdDir1_Dir2FwdDir2, map1, sign1);
                        locExpansion->GetFaceToElementMap(j, locExpansion->GetForient(j), map2, sign2);

                        for (k = 0; k < map1.size(); ++k)
                        {
                            // Find the elemental co-efficient in the original
                            // mapping.
                            int idx = -1;
                            for (int l = 0; l < map2.size(); ++l)
                            {
                                if (map1[k] == map2[l])
                                {
                                    idx = l;
                                    break;
                                }
                            }

                            ASSERTL2(idx != -1, "Problem with face to element map!");
                            orientMap[k] = idx;
                        }

                        for(k = 0; k < order_e; ++k)
                        {
                            vGlobalId = m_localToGlobalBndMap[k+cnt];
                            m_globalToUniversalBndMap[vGlobalId]
                                = id * maxDof + orientMap[k] + 1;
                        }
                        cnt += order_e;
                    }
                }
            }

            // Initialise GSlib and populate the unique map.
            Array<OneD, long> tmp(m_globalToUniversalBndMap.size());
            for (i = 0; i < m_globalToUniversalBndMap.size(); ++i)
            {
                tmp[i] = m_globalToUniversalBndMap[i];
            }
            m_bndGsh = m_gsh = Gs::Init(tmp, m_comm);
            Gs::Unique(tmp, m_comm);
            for (i = 0; i < m_globalToUniversalBndMap.size(); ++i)
            {
                m_globalToUniversalBndMapUnique[i] = (tmp[i] >= 0 ? 1 : 0);
            }
        }

        void AssemblyMapDG::RealignTraceElement(
            Array<OneD, int>        &toAlign,
            StdRegions::Orientation  orient,
            int                      nquad1,
            int                      nquad2)
        {
            if (orient == StdRegions::eBackwards)
            {
                ASSERTL1(nquad2 == 0, "nquad2 != 0 for reorienation");
                for (int i = 0; i < nquad1/2; ++i)
                {
                    swap(toAlign[i], toAlign[nquad1-1-i]);
                }
            }
            else if (orient != StdRegions::eForwards)
            {
                ASSERTL1(nquad2 != 0, "nquad2 == 0 for reorienation");

                Array<OneD, int> tmp(nquad1*nquad2);

                // Copy transpose.
                if (orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1; ++j)
                        {
                            tmp[i*nquad1 + j] = toAlign[j*nquad2 + i];
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1; ++j)
                        {
                            tmp[i*nquad1 + j] = toAlign[i*nquad1 + j];
                        }
                    }
                }

                if (orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
                    orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    // Reverse x direction
                    for (int i = 0; i < nquad2; ++i)
                    {
                        for (int j = 0; j < nquad1/2; ++j)
                        {
                            swap(tmp[i*nquad1 + j],
                                 tmp[i*nquad1 + nquad1-j-1]);
                        }
                    }
                }

                if (orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
                    orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
                    orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
                {
                    // Reverse y direction
                    for (int j = 0; j < nquad1; ++j)
                    {
                        for (int i = 0; i < nquad2/2; ++i)
                        {
                            swap(tmp[i*nquad1 + j],
                                 tmp[(nquad2-i-1)*nquad1 + j]);
                        }
                    }
                }
                Vmath::Vcopy(nquad1*nquad2, tmp, 1, toAlign, 1);
            }
        }

        int AssemblyMapDG::v_GetLocalToGlobalMap(const int i) const
        {
            return m_localToGlobalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMap(const int i) const
        {
            return m_globalToUniversalBndMap[i];
        }

        int AssemblyMapDG::v_GetGlobalToUniversalMapUnique(const int i) const
        {
            return m_globalToUniversalBndMapUnique[i];
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetLocalToGlobalMap()
        {
            return m_localToGlobalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMap()
        {
            return m_globalToUniversalBndMap;
        }

        const Array<OneD,const int>& AssemblyMapDG::v_GetGlobalToUniversalMapUnique()
        {
            return m_globalToUniversalBndMapUnique;
        }

        NekDouble AssemblyMapDG::v_GetLocalToGlobalSign(
                    const int i) const
        {
            return GetLocalToGlobalBndSign(i);
        }

        void AssemblyMapDG::v_LocalToGlobal(
                    const Array<OneD, const NekDouble>& loc,
                    Array<OneD,       NekDouble>& global,
                    bool useComm ) const
        {
            LocalBndToGlobal(loc,global,useComm);
        }

        void AssemblyMapDG::v_GlobalToLocal(
                    const Array<OneD, const NekDouble>& global,
                          Array<OneD,       NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        void AssemblyMapDG::v_GlobalToLocal(
                    const NekVector<NekDouble>& global,
                          NekVector<      NekDouble>& loc) const
        {
            GlobalToLocalBnd(global,loc);
        }

        void AssemblyMapDG::v_Assemble(
                    const Array<OneD, const NekDouble> &loc,
                          Array<OneD,       NekDouble> &global) const
        {
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_Assemble(
                    const NekVector<NekDouble>& loc,
                          NekVector<      NekDouble>& global) const
        {
            AssembleBnd(loc,global);
        }

        void AssemblyMapDG::v_UniversalAssemble(
                      Array<OneD,     NekDouble>& pGlobal) const
        {
            Gs::Gather(pGlobal, Gs::gs_add, m_gsh);
        }

        void AssemblyMapDG::v_UniversalAssemble(
                      NekVector<      NekDouble>& pGlobal) const
        {
            UniversalAssemble(pGlobal.GetPtr());
        }

        int AssemblyMapDG::v_GetFullSystemBandWidth() const
        {
            return GetBndSystemBandWidth();
        }

        int AssemblyMapDG::GetNumDirichletBndPhys()
        {
            return m_numDirichletBndPhys;
        }

        Array<OneD, LocalRegions::ExpansionSharedPtr>&
                    AssemblyMapDG::GetElmtToTrace(const int i)
        {
            ASSERTL1(i >= 0 && i < m_elmtToTrace.size(),
                     "i is out of range");
            return m_elmtToTrace[i];
        }

        Array<OneD, Array< OneD, LocalRegions::ExpansionSharedPtr> >&
                    AssemblyMapDG::GetElmtToTrace()
        {
            return m_elmtToTrace;
        }

        AssemblyCommDGSharedPtr AssemblyMapDG::GetAssemblyCommDG()
        {
            return m_assemblyComm;
        }

    } //namespace
} // namespace
