///////////////////////////////////////////////////////////////////////////////
//
// File LocTraceToTraceMap.cpp
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
// Description: Local trace to general trace mapping information
//
///////////////////////////////////////////////////////////////////////////////
#include <MultiRegions/AssemblyMap/LocTraceToTraceMap.h>
#include <MultiRegions/ExpList.h>

#include <LibUtilities/Foundations/ManagerAccess.h> 

namespace Nektar
{
    namespace MultiRegions
    {
        
        LocTraceToTraceMap::LocTraceToTraceMap(const ExpList &locExp,
                                 const ExpListSharedPtr &trace,
                                               const Array<OneD, Array<OneD, StdRegions::StdExpansionSharedPtr> >&elmtToTrace,
                                               const vector<bool> &LeftAdjacentFaces)
        {
            m_LocTraceToTraceMap = Array<OneD, Array<OneD, int> >(2);
            m_interpTrace = Array<OneD, Array<OneD, InterpLocTraceToTrace> >(2);
            m_interpTraceI0 = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
            m_interpTraceI1 = Array<OneD, Array<OneD, DNekMatSharedPtr> >(2);
            m_interpPoints  = Array<OneD, Array<OneD, TraceInterpPoints> >(2);
            m_interpNfaces  = Array<OneD, Array<OneD, int> >(2); 

            LocalRegions::Expansion3DSharedPtr exp3d;
            int cnt, n, e, phys_offset;
            
            const boost::shared_ptr<LocalRegions::ExpansionVector> exp =
                locExp.GetExp();
            
            int nexp = exp->size();
            // count number of faces and points required for maps 
            int nFwdPts = 0;
            int nBwdPts = 0; 
            m_nFwdLocTracePts = 0;
            m_nLocTracePts = 0; 
            m_nTracePts = trace->GetTotPoints();
            for(cnt = n = 0; n < nexp; ++n)
            {
                exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();
                
                for(int i = 0; i < exp3d->GetNfaces(); ++i, ++cnt)
                {
                    int nLocPts = exp3d->GetFaceNumPoints(i);
                    m_nLocTracePts += nLocPts; 
                    
                    if(LeftAdjacentFaces[cnt])
                    {
                        nFwdPts += elmtToTrace[n][i]->GetTotPoints();
                        m_nFwdLocTracePts += nLocPts; 
                    }
                    else
                    {
                        nBwdPts += elmtToTrace[n][i]->GetTotPoints();
                    }
                    
                }
            }
            m_fieldToLocTraceMap = Array<OneD, int>(m_nLocTracePts);
            
            m_LocTraceToTraceMap[0] = Array<OneD, int>(nFwdPts);
            m_LocTraceToTraceMap[1] = Array<OneD, int>(nBwdPts);
            
            // gather information about trace interpolations. 
            map<TraceInterpPoints, vector<pair<int,int> >, cmpop > TraceInterpMap;
            map<TraceInterpPoints, vector<pair<int,int> >, cmpop >::iterator it; 
            
            vector<vector<int> > TraceOrder;
            TraceOrder.resize(nexp);
            int nface;
            // generate a map of similar traces with the same
            // interpolation requirements.
            for(cnt = n = 0; n < nexp; ++n)
            {
                exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();
                nface = exp3d->GetNfaces();
                TraceOrder[n].resize(nface);
                for(e = 0; e < nface; ++e,++cnt)
                {
                    StdRegions::StdExpansionSharedPtr face = elmtToTrace[n][e]; 
                    StdRegions::Orientation orient = exp3d->GetFaceOrient(e);
                    
                    LibUtilities::PointsKey fromPointsKey0;
                    LibUtilities::PointsKey fromPointsKey1;
                    LibUtilities::PointsKey toPointsKey0;
                    LibUtilities::PointsKey toPointsKey1;
                    
                    int dir0 = exp3d->GetGeom3D()->GetDir(e,0);
                    int dir1 = exp3d->GetGeom3D()->GetDir(e,1);
                    
                    fromPointsKey0 = exp3d->GetBasis(dir0)->GetPointsKey();
                    fromPointsKey1 = exp3d->GetBasis(dir1)->GetPointsKey();
                    
                    if(orient < StdRegions::eDir1FwdDir2_Dir2FwdDir1)
                    {
                        toPointsKey0 = face->GetBasis(0)->GetPointsKey();
                        toPointsKey1 = face->GetBasis(1)->GetPointsKey();
                    }
                    else // transpose points key evaluation 
                    {
                        toPointsKey0 = face->GetBasis(1)->GetPointsKey();
                        toPointsKey1 = face->GetBasis(0)->GetPointsKey();
                    }
                    
                    TraceInterpPoints fpoint(fromPointsKey0,fromPointsKey1,
                                             toPointsKey0,toPointsKey1);
                    
                    pair<int,int> epf(n,e);
                    TraceInterpMap[fpoint].push_back(epf);
                    
                    TraceOrder[n][e] = cnt; 
                }
            }
            
            int nInterpType = TraceInterpMap.size();
            for(int i = 0; i < 2; ++i)
            {
                m_interpTrace[i]   = Array<OneD, InterpLocTraceToTrace>(nInterpType);
                m_interpTraceI0[i] = Array<OneD, DNekMatSharedPtr>     (nInterpType);
                m_interpTraceI1[i] = Array<OneD, DNekMatSharedPtr>     (nInterpType);
                m_interpPoints[i]  = Array<OneD, TraceInterpPoints>    (nInterpType);
                m_interpNfaces[i]  = Array<OneD, int>                  (nInterpType,0);
            }
            
            int nfacepts,nfacepts1; 
            int cnt1 = 0; 
            int cnt2 = 0; 
            int cntFwd = 0;
            int cntBwd = 0;
            int cntFwd1 = 0;
            int cntBwd1 = 0;
            int set;
            Array<OneD, int> faceids;
            Array<OneD, int> locTraceToTraceMap;
            cnt = 0 ; 
            for(it = TraceInterpMap.begin(); it != TraceInterpMap.end(); ++it,++cnt1)
            {
                LibUtilities::PointsKey fromPointsKey0 = it->first.get<0>();
                LibUtilities::PointsKey fromPointsKey1 = it->first.get<1>();
                LibUtilities::PointsKey toPointsKey0   = it->first.get<2>();
                LibUtilities::PointsKey toPointsKey1   = it->first.get<3>();
                
                bool fwdSet = false;
                bool bwdSet = false;
                
                for(int f = 0; f < it->second.size(); ++f,++cnt2)
                {
                    n = it->second[f].first;
                    e = it->second[f].second; 
                    
                    StdRegions::StdExpansionSharedPtr face = elmtToTrace[n][e]; 
                    
                    exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();
                    phys_offset = locExp.GetPhys_Offset(n);
                    
                    // mapping of new face order to one that loops
                    // over elmts then faces set up mapping of
                    // faces in standard cartesian order
                    exp3d->GetFacePhysMap(e, faceids);
                    nfacepts  = exp3d->GetFaceNumPoints(e);
                    nfacepts1 = face->GetTotPoints();
                    
                    StdRegions::Orientation orient = exp3d->GetFaceOrient(e);
                    
                    exp3d->ReOrientFacePhysMap(elmtToTrace[n][e]->GetNverts(), 
                                               orient,
                                               toPointsKey0.GetNumPoints(),
                                               toPointsKey1.GetNumPoints(),
                                               locTraceToTraceMap);
                    
                    int offset = trace->GetPhys_Offset(elmtToTrace[n][e]->GetElmtId());
                    
                    if(LeftAdjacentFaces[TraceOrder[n][e]])
                    {
                        
                        for(int i = 0; i < nfacepts; ++i)
                        {
                            m_fieldToLocTraceMap[cntFwd + i] = phys_offset + faceids[i];
                        }
                        
                        for(int i = 0; i < nfacepts1; ++i)
                        {
                            m_LocTraceToTraceMap[0][cntFwd1 + i] = offset + locTraceToTraceMap[i];
                        }
                        
                        cntFwd  += nfacepts; 
                        cntFwd1 += nfacepts1; 
                        set = 0; 
                    }
                    else
                    {
                        for(int i = 0; i < nfacepts; ++i)
                        {
                            m_fieldToLocTraceMap[m_nFwdLocTracePts + cntBwd +i] = phys_offset + faceids[i];
                        }
                        
                        for(int i = 0; i < nfacepts1; ++i)
                        {
                            m_LocTraceToTraceMap[1][cntBwd1 + i] = offset + locTraceToTraceMap[i];
                        }
                        
                        cntBwd  += nfacepts;                             
                        cntBwd1 += nfacepts1;                             
                        set = 1;
                    }
                    
                    m_interpNfaces[set][cnt1] += 1;  
                    

                    if(((fwdSet == false)&&(set == 0))||
                       ((bwdSet == false)&&(set == 1)))
                    {
                        m_interpPoints[set][cnt1] = it->first; 
                        
                        if(fromPointsKey0 == toPointsKey0)
                        {
                            if(fromPointsKey1 == toPointsKey1)
                            {
                                m_interpTrace[set][cnt1] = eNoInterp;
                            }
                            else
                            {
                                m_interpTrace[set][cnt1] = eInterpDir1;
                                m_interpTraceI1[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey1]->GetI(toPointsKey1);
                            }
                        }
                        else
                        {
                            if(fromPointsKey1 == toPointsKey1)
                            {
                                m_interpTrace[set][cnt1] = eInterpDir0;
                                m_interpTraceI0[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey0]-> GetI(toPointsKey0);          
                            }
                            else
                            {
                                m_interpTrace[set][cnt1] = eInterpBothDirs;
                                m_interpTraceI0[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey0]->GetI(toPointsKey0);                            
                                m_interpTraceI1[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey1]->GetI(toPointsKey1);           
                            }
                        }
                        
                        if(set == 0)
                        {
                            fwdSet = true;
                        }
                        else
                        {
                            bwdSet = true;
                        }
                    }
                }
                
            }
        }
        
        //desctructor
        LocTraceToTraceMap::~LocTraceToTraceMap()
        {};

        void LocTraceToTraceMap::LocTracesFromField(const Array<OneD, const NekDouble> &field, Array<OneD, NekDouble> faces)
        {
            Vmath::Gathr(m_fieldToLocTraceMap.num_elements(),
                         field,m_fieldToLocTraceMap,faces);
        }


        void LocTraceToTraceMap::FwdLocTracesFromField(const Array<OneD, const NekDouble> &field, Array<OneD, NekDouble> faces)
        {
            Vmath::Gathr(m_nFwdLocTracePts, field,m_fieldToLocTraceMap,faces);
        }


        /// interpolate local faces to trace face point distribution
        void LocTraceToTraceMap::InterpLocFacesToTrace(const int dir, const Array<OneD, const NekDouble> &locfaces, Array<OneD, NekDouble> faces)
        {
            ASSERTL1(dir < 2,"option dir out of range, dir=0 is fwd, dir=1 is bwd");

            int cnt1 = 0;
            int cnt  = 0;

            // tmp space assuming forward map is of size of trace
            Array<OneD, NekDouble> tmp(m_nTracePts);

            for(int i = 0; i < m_interpTrace[dir].num_elements(); ++i)
            {
                if(m_interpNfaces[dir][i]) // check there are faces to interpolate
                {
                    // get to/from points 
                    LibUtilities::PointsKey fromPointsKey0 = m_interpPoints[dir][i].get<0>();
                    LibUtilities::PointsKey fromPointsKey1 = m_interpPoints[dir][i].get<1>();
                    LibUtilities::PointsKey toPointsKey0 = m_interpPoints[dir][i].get<2>();
                    LibUtilities::PointsKey toPointsKey1 = m_interpPoints[dir][i].get<3>();
                    int fnp0 = fromPointsKey0.GetNumPoints();
                    int fnp1 = fromPointsKey1.GetNumPoints();
                    int tnp0 = toPointsKey0.GetNumPoints();
                    int tnp1 = toPointsKey1.GetNumPoints();
                
                    int nfromfacepts = m_interpNfaces[dir][i]*fnp0*fnp1;
                    
                    // Do interpolation here if required
                    switch(m_interpTrace[dir][i])
                    {
                    case eNoInterp: // Just copy
                        {
                            Vmath::Vcopy(nfromfacepts,locfaces.get()+cnt,1,
                                         tmp.get()+cnt1,1);
                        }
                        break;
                    case eInterpDir0:
                        {
                            DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                            Blas::Dgemm('N', 'N', tnp0, tnp1, fnp0, 1.0, 
                                        I0->GetPtr().get(),
                                        tnp0, locfaces.get()+cnt, fnp0, 0.0, 
                                        tmp.get()+cnt1, tnp0);  
                        }
                        break;
                    case eInterpDir1:
                        {
                            DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];
                            for(int j = 0; j <  m_interpNfaces[dir][i]; ++j)
                            {
                            Blas::Dgemm('N', 'T', tnp0, tnp1, fnp1, 1.0, 
                                        locfaces.get()+cnt +j*fnp0*fnp1, tnp0, 
                                        I1->GetPtr().get(), 
                                        tnp1, 0.0,  tmp.get()+cnt1 + j*tnp0*tnp1, tnp0); 
                        }
                        }
                        break;
                    case eInterpBothDirs:
                        {
                            DNekMatSharedPtr I0 = m_interpTraceI0[dir][i];
                            DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];
                            Array<OneD, NekDouble > wsp(m_interpNfaces[dir][i]*fnp0*tnp1*fnp0);
                            
                            for(int j = 0; j <  m_interpNfaces[dir][i]; ++j)
                            {
                                Blas::Dgemm('N', 'T', fnp0, tnp1, fnp1, 1.0,
                                            locfaces.get()+cnt+j*fnp0*fnp1, fnp0,
                                            I1->GetPtr().get(), 
                                            tnp1, 0.0,  wsp.get()+j*fnp0*tnp1, fnp0);     
                        }
                            
                            Blas::Dgemm('N', 'N', tnp0, tnp1*m_interpNfaces[dir][i], fnp0, 1.0, 
                                        I0->GetPtr().get(),
                                        tnp0, wsp.get(), fnp0, 0.0, 
                                        tmp.get()+cnt1, tnp0);     
                        }
                        break;
                    }
                    cnt  += nfromfacepts; 
                    cnt1 += m_interpNfaces[dir][i]*tnp0*tnp1;
                }
            }

            Vmath::Scatr(m_LocTraceToTraceMap[dir].num_elements(),
                         tmp.get(),m_LocTraceToTraceMap[dir].get(),
                         faces.get());
        }
    } //namespace
} // namespace
