///////////////////////////////////////////////////////////////////////////////
//
// File LocTraceToTraceMap.h
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
// Description: Local Trace to general trace map information
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCTRACETOTRACEMAP_H
#define MULTIREGIONS_LOCTRACETOTRACEMAP_H

#include <LocalRegions/Expansion.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        enum InterpLocTraceToTrace
        {
            eNoInterp,
            eInterpDir0,
            eInterpEndPtDir0,
            eInterpDir1,
            eInterpEndPtDir1,
            eInterpBothDirs,
            eInterpEndPtDir0InterpDir1
        };

        typedef boost::tuple<
            LibUtilities::PointsKey,
            LibUtilities::PointsKey,
            LibUtilities::PointsKey,
            LibUtilities::PointsKey> TraceInterpPoints;
        
        struct cmpop
        {
            bool operator()(
                TraceInterpPoints const &a,
                TraceInterpPoints const &b) const
            {
                if (a.get<0>() < b.get<0>())
                {
                    return true;
                }

                if (b.get<0>() < a.get<0>())
                {
                    return false;
                }

                if (a.get<1>() < b.get<1>())
                {
                    return true;
                }
                if (b.get<1>() < a.get<1>())
                {
                    return false;
                }

                if (a.get<2>() < b.get<2>())
                {
                    return true;
                }

                if (b.get<2>() < a.get<2>())
                {
                    return false;
                }

                if (a.get<3>() < b.get<3>())
                {
                    return true;
                }

                
                return false;
            }
        };
        
        class LocTraceToTraceMap
        {
        public:
            // Constructor
            MULTI_REGIONS_EXPORT LocTraceToTraceMap(
                const ExpList &locExp,
                const ExpListSharedPtr &trace,
                const Array<OneD, Array<OneD, LocalRegions::
                    ExpansionSharedPtr> >&elmtToTrace,
                const vector<bool> &LeftAdjacents);
            
            // Destructor
            MULTI_REGIONS_EXPORT virtual ~LocTraceToTraceMap();
            
            MULTI_REGIONS_EXPORT void Setup2D(
                const ExpList &locExp,
                const ExpListSharedPtr &trace,
                const Array<OneD, Array<OneD, LocalRegions::
                    ExpansionSharedPtr> >&elmtToTrace,
                const vector<bool> &LeftAdjacents);
            
            MULTI_REGIONS_EXPORT void Setup3D(
                const ExpList &locExp,
                const ExpListSharedPtr &trace,
                const Array<OneD, Array<OneD, LocalRegions::
                    ExpansionSharedPtr> >&elmtToTrace,
                const vector<bool> &LeftAdjacents);
            
            MULTI_REGIONS_EXPORT void LocTracesFromField(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> faces);

            
            MULTI_REGIONS_EXPORT void FwdLocTracesFromField(
                const Array<OneD, const NekDouble> &field,
                      Array<OneD,       NekDouble> faces);
            
            MULTI_REGIONS_EXPORT void InterpLocEdgesToTrace(
                const int dir,
                const Array<OneD, const NekDouble> &locfaces,
                      Array<OneD,       NekDouble> edges);
            
            MULTI_REGIONS_EXPORT void InterpLocFacesToTrace(
                const int dir,
                const Array<OneD, const NekDouble> &locfaces,
                      Array<OneD,       NekDouble> faces);

            MULTI_REGIONS_EXPORT void AddTraceCoeffsToFieldCoeffs(
                const Array<OneD, const NekDouble> &trace,
                      Array<OneD,       NekDouble> &field);

            MULTI_REGIONS_EXPORT void AddTraceCoeffsToFieldCoeffs(
                const int dir,
                const Array<OneD, const NekDouble> &race,
                      Array<OneD,       NekDouble> &field);

            MULTI_REGIONS_EXPORT inline int GetNFwdLocTracePts()
            {
                return m_nFwdLocTracePts; 
            }

            MULTI_REGIONS_EXPORT inline int GetNLocTracePts()
            {
                return m_nLocTracePts; 
            }
        private:
            // the number of data points in m_fieldtoLocTraceMap
            // that are associated with Fwd Trace
            int m_nFwdLocTracePts;  
            int m_nLocTracePts;
            int m_nTracePts;
            int m_expdim;
            
            // mapping from field to local trace
            Array<OneD, int> m_fieldToLocTraceMap; 
            
            Array<OneD, Array<OneD, int> > m_LocTraceToTraceMap;
            Array<OneD, Array<OneD, InterpLocTraceToTrace> >   m_interpTrace;
            Array<OneD, Array<OneD, DNekMatSharedPtr> >        m_interpTraceI0;
            Array<OneD, Array<OneD, DNekMatSharedPtr> >        m_interpTraceI1;
            Array<OneD, Array<OneD, TraceInterpPoints> >       m_interpPoints;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_interpEndPtI0;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_interpEndPtI1;
            Array<OneD, Array<OneD, int> > m_interpNfaces; 
            
            int m_nTraceCoeffs[2];
            Array<OneD, Array<OneD, int> > m_traceCoeffsToElmtMap;
            Array<OneD, Array<OneD, int> > m_traceCoeffsToElmtTrace;
            Array<OneD, Array<OneD, int> > m_traceCoeffsToElmtSign;
        };

        typedef boost::shared_ptr<LocTraceToTraceMap>
            LocTraceToTraceMapSharedPtr;

    } // end of namespace
} // end of namespace

#endif 
