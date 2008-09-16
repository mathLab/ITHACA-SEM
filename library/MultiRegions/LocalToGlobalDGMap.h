///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalDGMap.h
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
// Description: Local to Global DG mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////
#ifndef MULTIREGIONS_LOCALTOGLOBALDGMAP_H
#define MULTIREGIONS_LOCALTOGLOBALDGMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/GenExpList1D.h>

#include <LocalRegions/PointExp.h>

namespace Nektar
{
    namespace MultiRegions 
    {

        class LocalToGlobalDGMap: public LocalToGlobalBaseMap
        {
        public:
            LocalToGlobalDGMap();

            ~LocalToGlobalDGMap();

            LocalToGlobalDGMap( const SpatialDomains::MeshGraph1D &graph1D,
                                const boost::shared_ptr<StdRegions::StdExpansionVector> &exp1D, 
                                const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndConstraint,
                                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond);
            
            LocalToGlobalDGMap(SpatialDomains::MeshGraph2D &graph2D, 
                               const GenExpList1DSharedPtr &trace, 
                               const boost::shared_ptr<StdRegions::StdExpansionVector> &exp2D, 
                               const Array<OneD, MultiRegions::ExpList1DSharedPtr> &bndContraint, 
                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond,
                               const map<int,int> &periodicEdges);
            
            int GetNumDirichletBndPhys()
            {
                return m_numDirichletBndPhys;
            }

            Array<OneD, StdRegions::StdExpansion1DSharedPtr> GetElmtToTrace(const int i)
            {
                ASSERTL1(i >= 0 && i < m_elmtToTrace.num_elements(),
                         "i is out of range");
                return m_elmtToTrace[i];
            }

            Array<OneD, Array< OneD, StdRegions::StdExpansion1DSharedPtr> > GetElmtToTrace()
            {
                return m_elmtToTrace;
            }
            
            AdjacentTraceOrientation GetBndExpAdjacentOrient(const int i)
            {
                return m_bndExpAdjacentOrient[i];
            }
            
        protected:

        private:
            int m_numDirichletBndPhys;  // Number of physical dirichlet boundary values in trace

            Array<OneD, Array<OneD, StdRegions::StdExpansion1DSharedPtr> > m_elmtToTrace;  //< list of edge expansions for a given element 
            
            Array<OneD, AdjacentTraceOrientation > m_bndExpAdjacentOrient;
        };
        
        typedef boost::shared_ptr<LocalToGlobalDGMap>  LocalToGlobalDGMapSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //LOCALTOGLOBALDGMAP_H


/** $Log: LocalToGlobalDGMap.h,v $
/** Revision 1.1  2008/08/18 08:16:23  sherwin
/** First version of this new class container for mappings
/**
 */

