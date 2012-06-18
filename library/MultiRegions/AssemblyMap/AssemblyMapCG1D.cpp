///////////////////////////////////////////////////////////////////////////////
//
// File AssemblyMapCG1D.cpp
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
// Description: C0-continuous assembly mappings specific to 1D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMapCG1D.h>
#include <LocalRegions/PointExp.h>
#include <LocalRegions/SegExp.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class AssemblyMapCG1D
         * Mappings are created for three possible global solution types:
         *  - Direct full matrix
         *  - Direct static condensation
         *  - Direct multi-level static condensation
         * In the latter case, mappings are created recursively for the
         * different levels of static condensation.
         *
         * These mappings are used by GlobalLinSys to generate the global
         * system.
         */

        /**
         *
         */
        AssemblyMapCG1D::AssemblyMapCG1D(
                const LibUtilities::SessionReaderSharedPtr &pSession):
            AssemblyMapCG(pSession)
        {
        }


        /**
         *
         */
        AssemblyMapCG1D::AssemblyMapCG1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const ExpListSharedPtr>
                                                            &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                            &bndConditions,
                const map<int,int>& periodicVerticesId):
            AssemblyMapCG(pSession)
        {
            SetUp1DExpansionC0ContMap(numLocalCoeffs,
                                      locExp,
                                      bndCondExp,
                                      bndConditions,
                                      periodicVerticesId);

            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG1D::AssemblyMapCG1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const int numLocalCoeffs,
                const ExpList &locExp):
            AssemblyMapCG(pSession)
        {
            SetUp1DExpansionC0ContMap(numLocalCoeffs, locExp);
            CalculateBndSystemBandWidth();
            CalculateFullSystemBandWidth();
        }


        /**
         *
         */
        AssemblyMapCG1D::~AssemblyMapCG1D()
        {
        }


        /**
         * Construction of the local->global map is achieved in several stages.
         * A mesh vertex renumbering is constructed as follows
         *  - Domain Dirichlet boundaries are numbered first.
         *  - Domain periodic boundaries are numbered next and given
         *    consistent global indices.
         *  - Remaining vertices are ordered last.
         * This mesh vertex renumbering is then used to generate the first
         * part of the local to global mapping of the degrees of freedom. The
         * remainder of the map consists of the element-interior degrees of
         * freedom. This leads to the block-diagonal submatrix as each
         * element-interior mode is globally orthogonal to modes in all other
         * elements.
         *
         * The boundary condition mapping is generated from the same vertex
         * renumbering.
         */
        void AssemblyMapCG1D::SetUp1DExpansionC0ContMap(
                const int numLocalCoeffs,
                const ExpList &locExp,
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                                                                &bndCondExp,
                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                &bndConditions,
                const map<int,int>& periodicVerticesId)
        {
            int i,j;
            int cnt = 0;
            int meshVertId;
            int meshVertId2;
            int graphVertId = 0;
            int globalId;

            const StdRegions::StdExpansionVector &locExpVector = *(locExp.GetExp());

            LocalRegions::SegExpSharedPtr locSegExp;

            int nbnd = bndCondExp.num_elements();
            m_staticCondLevel                = 0;
            m_signChange                     = false;
            m_numPatches                     = locExpVector.size();
            m_numLocalCoeffs                 = numLocalCoeffs;
            m_numLocalBndCoeffs              = 2*locExpVector.size();
            m_localToGlobalMap               = Array<OneD, int>(m_numLocalCoeffs,-1);
            m_localToGlobalBndMap            = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int>(nbnd);
            m_numLocalBndCoeffsPerPatch      = Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch      = Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i)
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) locExpVector[locExp.GetOffset_Elmt_Id(i)]->GetNcoeffs() - locExpVector[locExp.GetOffset_Elmt_Id(i)]->NumBndryCoeffs();
            }

            // The only unique identifiers of the vertices of the mesh
            // are the vertex id (stored in their corresponding
            // Geometry object).  However, setting up a global
            // numbering based on these id's will not lead to a
            // suitable or optimal numbering. Mainly because: - we
            // want the Dirichlet DOF's to be listed first - we want
            // an optimal global numbering of the remaining DOF's
            // (strategy still need to be defined but can for example
            // be: minimum bandwith or minimum fill-in of the
            // resulting global system matrix)
            //
            // That's why the vertices be rearranged. Currently, this
            // is done in the following way: The vertices of the mesh
            // are considered as vertices of a graph. (We then will
            // use algorithms to reorder these graph-vertices -
            // although not implemented yet in 1D).
            //
            // A container is used to store the graph vertex id's of
            // the different mesh vertices. This is implemented as a
            // STL map such that the graph vertex id can later be
            // retrieved by the unique mesh vertex id's which serve as
            // the key of the map.
            map<int, int> vertReorderedGraphVertId;
            map<int,int>::const_iterator mapConstIt;

            // STEP 1: Order the Dirichlet vertices first
            m_numGlobalDirBndCoeffs = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndConditions[i]->GetBoundaryConditionType()==SpatialDomains::eDirichlet)
                {
                    meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                    vertReorderedGraphVertId[meshVertId] = graphVertId++;
                    m_numGlobalDirBndCoeffs++;
                    m_numLocalDirBndCoeffs++;
                }
            }

            // STEP 2: Order the periodic vertices next
            // This allows to give corresponding DOF's the same
            // global ID
            for(mapConstIt = periodicVerticesId.begin(); mapConstIt != periodicVerticesId.end(); mapConstIt++)
            {
                meshVertId  = mapConstIt->first;
                meshVertId2 = mapConstIt->second;

                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");
                ASSERTL0(vertReorderedGraphVertId.count(meshVertId) == 0,
                         "This periodic boundary vertex has been specified before");

                vertReorderedGraphVertId[meshVertId]  = graphVertId;
                vertReorderedGraphVertId[meshVertId2] = graphVertId++;
            }


            // STEP 3: List the other vertices
            // STEP 4: Set up simple map based on vertex and edge id's
            for(i = 0; i < locExpVector.size(); ++i)
            {
                cnt = locExp.GetCoeff_Offset(i);
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(locExpVector[i]))
                {
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {
                        meshVertId = (locSegExp->GetGeom1D())->GetVid(j);
                        if(vertReorderedGraphVertId.count(meshVertId) == 0)
                        {
                            vertReorderedGraphVertId[meshVertId] = graphVertId++;
                        }

                        m_localToGlobalMap[cnt+locSegExp->GetVertexMap(j)] =
                            vertReorderedGraphVertId[meshVertId];
                    }
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }

            // Set up the mapping for the boundary conditions
            for(i = 0; i < nbnd; i++)
            {
                meshVertId = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = vertReorderedGraphVertId[meshVertId];
            }

            // Setup interior mapping and the boundary map
            globalId = graphVertId;
            m_numGlobalBndCoeffs = globalId;

            cnt=0;
            for(i = 0; i < m_numLocalCoeffs; ++i)
            {
                if(m_localToGlobalMap[i] == -1)
                {
                    m_localToGlobalMap[i] = globalId++;
                }
                else
                {
                    m_localToGlobalBndMap[cnt++]=m_localToGlobalMap[i];
                }
            }
            m_numGlobalCoeffs = globalId;

            SetUpUniversalC0ContMap(locExp);

            m_hash = boost::hash_range(m_localToGlobalMap.begin(), m_localToGlobalMap.end());
        }




    } // namespace
} // namespace
