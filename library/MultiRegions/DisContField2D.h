///////////////////////////////////////////////////////////////////////////////
//
// File DisContField2D.h
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
// Description: Field definition in two-dimensions for a discontinuous
// LDG-H expansion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H
#define NEKTAR_LIBS_MULTIREGIONS_DISCONTFIELD2D_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/GenExpList1D.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/BoundaryConditions.h>
#include <SpatialDomains/SegGeom.h>
#include <MultiRegions/GlobalLinSys.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class DisContField2D: public ExpList2D
        {
        public:
            DisContField2D();

            DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                SpatialDomains::BoundaryConditions &bcs, 
                           const int bc_loc = 0);

            DisContField2D(SpatialDomains::MeshGraph2D &graph2D,
                           SpatialDomains::BoundaryConditions &bcs, 
                           const std::string variable);
            
            DisContField2D(const DisContField2D &In);

            ~DisContField2D();
            
            void GenerateFieldBnd2D(SpatialDomains::MeshGraph2D &graph2D,
                                    SpatialDomains::BoundaryConditions &bcs, 
                                    const std::string variable);

            void SetUpTraceMappings(SpatialDomains::MeshGraph2D &graph2D);

            inline GenExpList1DSharedPtr &GetTrace(void)
            {
                return m_trace;
            }

            
            void HelmSolve(DisContField2D &Fce, NekDouble lambda);

            GlobalLinSysSharedPtr GenGlobalBndLinSys(const GlobalLinSysKey &mkey);
            GlobalLinSysSharedPtr GetGlobalBndLinSys(const GlobalLinSysKey &mkey);
        
            void GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                    Array<OneD,NekDouble> &Bwd);

            void ExtractTracePhys()
            {
                ExtractTracePhys(m_trace->UpdatePhys());
            }

            void ExtractTracePhys(Array<OneD,NekDouble> &outarray);

            void AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                                  Array<OneD, const NekDouble> &Fy, 
                                  Array<OneD, NekDouble> &outarray);

        protected:

        private:
            int m_numTraceDirichletBCs;       ///< Number of Coeff space BCs
            int m_numTraceDirichletPhysBCs;   ///< Number of physical space BCs
            Array<OneD,MultiRegions::ExpList1DSharedPtr>       m_bndConstraint;
            Array<OneD,SpatialDomains::BoundaryConditionType>  m_bndTypes;
            GlobalLinSysMapShPtr                               m_globalBndMat;
            GenExpList1DSharedPtr                              m_trace;

            Array<OneD,Array<OneD,LocalRegions::GenSegExpSharedPtr> > m_elmtToTrace;
            // NOTE This should all go into a class structure
            Array<OneD,Array<OneD, int> > m_bndEidToTraceEid;  ///< Boundary list Expansion ID to Trace list Expansion ID
            Array<OneD,Array<OneD, AdjacentEdgeOrientation> > m_bndExpAdjacentOrient;  ///< Boundary Expansion adjacent edge orientation 
            Array<OneD, int > m_elmtTraceMap;
            Array<OneD, int > m_elmtTraceSign;

        };

        typedef boost::shared_ptr<DisContField2D>   DisContField2DSharedPtr;

    } //end of namespace
} //end of namespace

#endif // MULTIERGIONS_DISCONTFIELD2D_H
