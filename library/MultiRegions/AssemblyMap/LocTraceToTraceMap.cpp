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
#include <LocalRegions/Expansion2D.h> 

namespace Nektar
{
    namespace MultiRegions
    {
        
        LocTraceToTraceMap::LocTraceToTraceMap(const ExpList &locExp,
                                 const ExpListSharedPtr &trace,
                                               const Array<OneD, Array<OneD, LocalRegions::ExpansionSharedPtr> >&elmtToTrace,
                                               const vector<bool> &LeftAdjacentFaces)
        {
            m_LocTraceToTraceMap = Array<OneD, Array<OneD, int> >(2);
            m_interpTrace = Array<OneD, Array<OneD, InterpLocTraceToTrace> >(2);
            m_interpTraceI0 = Array<OneD, Array<OneD, DNekMatSharedPtr> > (2);
            m_interpTraceI1 = Array<OneD, Array<OneD, DNekMatSharedPtr> > (2);
            m_interpEndPtI0 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
            m_interpEndPtI1 = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(2);
            m_interpPoints  = Array<OneD, Array<OneD, TraceInterpPoints> >(2);
            m_interpNfaces  = Array<OneD, Array<OneD, int> >(2); 

            m_traceCoeffsToElmtMap   = Array<OneD, Array<OneD, int> >(2);
            m_traceCoeffsToElmtTrace = Array<OneD, Array<OneD, int> >(2);
            m_traceCoeffsToElmtSign  = Array<OneD, Array<OneD, int> >(2);

            LocalRegions::Expansion3DSharedPtr exp3d;
            int cnt, n, e, phys_offset;
            
            const boost::shared_ptr<LocalRegions::ExpansionVector> exp =
                locExp.GetExp();
            
            int nexp = exp->size();
            // count number of faces and points required for maps 
            int nFwdPts = 0;
            int nBwdPts = 0; 
            int nFwdCoeffs = 0;
            int nBwdCoeffs = 0; 
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
                        nFwdPts    += elmtToTrace[n][i]->GetTotPoints();
                        nFwdCoeffs += elmtToTrace[n][i]->GetNcoeffs();
                        m_nFwdLocTracePts += nLocPts; 
                    }
                    else
                    {
                        nBwdPts += elmtToTrace[n][i]->GetTotPoints();
                        nBwdCoeffs += elmtToTrace[n][i]->GetNcoeffs();
                    }
                    
                }
            }

            m_fieldToLocTraceMap = Array<OneD, int>(m_nLocTracePts);
            
            m_LocTraceToTraceMap[0] = Array<OneD, int>(nFwdPts);
            m_LocTraceToTraceMap[1] = Array<OneD, int>(nBwdPts);

            m_nTraceCoeffs[0] = nFwdCoeffs;
            m_nTraceCoeffs[1] = nBwdCoeffs;
            m_traceCoeffsToElmtMap[0]   = Array<OneD, int>(nFwdCoeffs+nBwdCoeffs);
            m_traceCoeffsToElmtMap[1]   = m_traceCoeffsToElmtMap[0] + nFwdCoeffs;
            m_traceCoeffsToElmtTrace[0] = Array<OneD, int>(nFwdCoeffs+nBwdCoeffs);
            m_traceCoeffsToElmtTrace[1] = m_traceCoeffsToElmtTrace[0] + nFwdCoeffs;
            m_traceCoeffsToElmtSign[0]  = Array<OneD, int>(nFwdCoeffs+nBwdCoeffs);
            m_traceCoeffsToElmtSign[1]  = m_traceCoeffsToElmtSign[0] + nFwdCoeffs;
            
            // gather information about trace interpolations. 
            map<TraceInterpPoints, vector<pair<int,int> >, cmpop > TraceInterpMap;
            map<TraceInterpPoints, vector<pair<int,int> >, cmpop >::iterator it; 
            
            vector<vector<int> > TraceOrder;
            TraceOrder.resize(nexp);
            int nface;
            int fwdcnt = 0;
            int bwdcnt = 0;
            // generate a map of similar traces with the same
            // interpolation requirements.
            for(cnt = n = 0; n < nexp; ++n)
            {
                exp3d = (*exp)[n]->as<LocalRegions::Expansion3D>();
                nface = exp3d->GetNfaces();
                TraceOrder[n].resize(nface);

                int coeffoffset = locExp.GetCoeff_Offset(n);
                for(e = 0; e < nface; ++e,++cnt)
                {
                    StdRegions::StdExpansionSharedPtr face = elmtToTrace[n][e]; 
                    StdRegions::Orientation orient = exp3d->GetForient(e);
                    
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

                    // setup for coefficient mapping from trace normal flux to elements
                    Array<OneD, unsigned int> map;
                    Array<OneD, int> sign;
                    exp3d->GetFaceToElementMap(e, orient, map, sign);

                    int order_f = face->GetNcoeffs();
                    int foffset = trace->GetCoeff_Offset(face->GetElmtId());
                    
                    int fac = (*exp)[n]->FaceNormalNegated(e) ? -1.0 : 1.0;
                    
                    if (exp3d->GetFaceExp(e)->GetRightAdjacentElementExp())
                    {
                        if (exp3d->GetFaceExp(e)->GetRightAdjacentElementExp()->GetGeom3D()
                            ->GetGlobalID() == exp3d->GetGeom3D()->GetGlobalID())
                        {
                            fac = -1.0; 
                        }
                    }

                    if(LeftAdjacentFaces[cnt])
                    {
                        for(int i = 0; i < order_f; ++i)
                        {
                            m_traceCoeffsToElmtMap  [0][fwdcnt]   = coeffoffset + map[i];
                            m_traceCoeffsToElmtTrace[0][fwdcnt]   = foffset + i;
                            m_traceCoeffsToElmtSign [0][fwdcnt++] = fac*sign[i];
                        }
                    }
                    else
                    {
                        for(int i = 0; i < order_f; ++i)
                        {
                            m_traceCoeffsToElmtMap  [1][bwdcnt]   = coeffoffset + map[i];
                            m_traceCoeffsToElmtTrace[1][bwdcnt]   = foffset + i;
                            m_traceCoeffsToElmtSign [1][bwdcnt++] = fac*sign[i];
                        }
                    }
                }
            }
            
            int nInterpType = TraceInterpMap.size();
            for(int i = 0; i < 2; ++i)
            {
                m_interpTrace[i]   = Array<OneD, InterpLocTraceToTrace>(nInterpType);
                m_interpTraceI0[i] = Array<OneD, DNekMatSharedPtr>     (nInterpType);
                m_interpTraceI1[i] = Array<OneD, DNekMatSharedPtr>     (nInterpType);
                m_interpEndPtI0[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
                m_interpEndPtI1[i] = Array<OneD, Array<OneD, NekDouble> >(nInterpType);
                m_interpPoints[i]  = Array<OneD, TraceInterpPoints>    (nInterpType);
                m_interpNfaces[i]  = Array<OneD, int>                (nInterpType,0);
            }
            
            int nfacepts,nfacepts1; 
            int cnt1    = 0; 
            int cnt2    = 0; 
            int cntFwd  = 0;
            int cntBwd  = 0;
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
                    
                    StdRegions::Orientation orient = exp3d->GetForient(e);
                    
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

                                // check to see if we can just
                                // interpolate endpoint
                                if((fromPointsKey1.GetPointsType() == 
                                    LibUtilities::eGaussRadauMAlpha1Beta0)
                                   &&(toPointsKey1.GetPointsType() == 
                                      LibUtilities::eGaussLobattoLegendre))
                                {
                                    if(fromPointsKey1.GetNumPoints() +1 ==
                                       toPointsKey1.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] = eInterpEndPtDir1;
                                        int fnp1 = fromPointsKey1.GetNumPoints();
                                        int tnp1 = toPointsKey1.GetNumPoints();
                                        m_interpEndPtI1[set][cnt1] = Array<OneD, NekDouble>(fnp1);
                                        Vmath::Vcopy(fnp1,
                                                     m_interpTraceI1[set][cnt1]->GetPtr().get()
                                                     +tnp1-1,tnp1,
                                                     &m_interpEndPtI1[set][cnt1][0],1);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(fromPointsKey1 == toPointsKey1)
                            {
                                m_interpTrace[set][cnt1] = eInterpDir0;
                                m_interpTraceI0[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey0]-> GetI(toPointsKey0); 

                                // check to see if we can just
                                // interpolate endpoint
                                if((fromPointsKey0.GetPointsType() == 
                                    LibUtilities::eGaussRadauMAlpha1Beta0)
                                   &&(toPointsKey0.GetPointsType() == 
                                      LibUtilities::eGaussLobattoLegendre))
                                {
                                    if(fromPointsKey0.GetNumPoints() +1 ==
                                       toPointsKey0.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] = eInterpEndPtDir0;

                                        int fnp0 = fromPointsKey0.GetNumPoints();
                                        int tnp0 = toPointsKey0.GetNumPoints();
                                        m_interpEndPtI0[set][cnt1] = Array<OneD, NekDouble>(fnp0);
                                        Vmath::Vcopy(fnp0,
                                                     m_interpTraceI0[set][cnt1]->GetPtr().get()
                                                     +tnp0-1,tnp0,
                                                     &m_interpEndPtI0[set][cnt1][0],1);
                                    }
                                }


                            }
                            else
                            {
                                m_interpTrace[set][cnt1] = eInterpBothDirs;
                                m_interpTraceI0[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey0]->GetI(toPointsKey0);                            
                                m_interpTraceI1[set][cnt1] = LibUtilities::PointsManager()[fromPointsKey1]->GetI(toPointsKey1);           

                                // check to see if we can just
                                // interpolate endpoint
                                if((fromPointsKey0.GetPointsType() == 
                                    LibUtilities::eGaussRadauMAlpha1Beta0)
                                   &&(toPointsKey0.GetPointsType() == 
                                      LibUtilities::eGaussLobattoLegendre))
                                {
                                    if(fromPointsKey0.GetNumPoints() +1 ==
                                       toPointsKey0.GetNumPoints())
                                    {
                                        m_interpTrace[set][cnt1] = eInterpEndPtDir0InterpDir1;
                                        int fnp0 = fromPointsKey0.GetNumPoints();
                                        int tnp0 = toPointsKey0.GetNumPoints();
                                        m_interpEndPtI0[set][cnt1] = Array<OneD, NekDouble>(fnp0);
                                        Vmath::Vcopy(fnp0,
                                                     m_interpTraceI0[set][cnt1]->GetPtr().get()
                                                     +tnp0-1,tnp0,
                                                     &m_interpEndPtI0[set][cnt1][0],1);
                                    }
                                }
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
                    case eInterpEndPtDir0:
                        {
                            int nfaces = m_interpNfaces[dir][i];
                            for(int k = 0; k < fnp0; ++k)
                            {
                                Vmath::Vcopy(nfaces*fnp1,locfaces.get()+cnt+k,
                                             fnp0,tmp.get()+cnt1+k,tnp0);
                            }
                            Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
                            Blas::Dgemv('T',  fnp0, tnp1*m_interpNfaces[dir][i],
                                        1.0, tmp.get()+cnt1, tnp0,I0.get(),1,  0.0, 
                                        tmp.get()+cnt1+ tnp0-1, tnp0);     
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
                    case eInterpEndPtDir1:
                        {
                            Array<OneD, NekDouble> I1 = m_interpEndPtI1[dir][i];
                            for(int j = 0; j < m_interpNfaces[dir][i]; ++j)
                            {
                                // copy all points
                                Vmath::Vcopy(fnp0*fnp1, locfaces.get() + cnt + 
                                             j*fnp0*fnp1,1,
                                             tmp.get()+cnt1+j*tnp0*tnp1,1);
                                
                                // interpolate end points 
                                for(int k = 0; k < tnp0; ++k)
                                {
#if 0 
                                    tmp[cnt1+k+(j+1)*tnp0*tnp1-tnp0]
                                        = Blas::Ddot(fnp1,locfaces.get()+cnt+
                                                     j*fnp0*fnp1+k,fnp0,
                                                     I1->GetPtr().get()+tnp1-1,
                                                     tnp1);
#else
                                    tmp[cnt1+k+(j+1)*tnp0*tnp1-tnp0]
                                        = Blas::Ddot(fnp1,locfaces.get()+cnt+
                                                     j*fnp0*fnp1+k,fnp0,&I1[0],1);
#endif
                                }
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
                            
                            Blas::Dgemm('N', 'N', tnp0, 
                                        tnp1*m_interpNfaces[dir][i], fnp0, 1.0, 
                                        I0->GetPtr().get(),
                                        tnp0, wsp.get(), fnp0, 0.0, 
                                        tmp.get()+cnt1, tnp0);     
                        }
                        break;
                    case eInterpEndPtDir0InterpDir1:
                        {
                            DNekMatSharedPtr I1 = m_interpTraceI1[dir][i];
                            
                            for(int j = 0; j <  m_interpNfaces[dir][i]; ++j)
                            {
                                Blas::Dgemm('N', 'T', fnp0, tnp1, fnp1, 1.0,
                                            locfaces.get()+cnt+j*fnp0*fnp1, fnp0,
                                            I1->GetPtr().get(), tnp1, 0.0,  
                                            tmp.get()+cnt1+j*tnp0*tnp1, tnp0);     
                            }

                            Array<OneD, NekDouble> I0 = m_interpEndPtI0[dir][i];
#if 0
                            for(int j = 0; j < tnp1*m_interpNfaces[dir][i]; ++j)
                            {
                                tmp[cnt1+(j+1)*tnp0-1] = Blas::Ddot(fnp0,tmp.get()+
                                                                    cnt1 + j*tnp0,1,
                                                                    I0.get(),1);
                            }
#else
                            Blas::Dgemv('T',  fnp0, tnp1*m_interpNfaces[dir][i],
                                        1.0, tmp.get()+cnt1, tnp0,I0.get(),1,  0.0, 
                                        tmp.get()+cnt1+ tnp0-1, tnp0);     
#endif
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

        void LocTraceToTraceMap::AddTraceCoeffsToFieldCoeffs(const Array<OneD, const NekDouble> &trace, Array<OneD, NekDouble> &field)
        {
            int nvals = m_nTraceCoeffs[0] + m_nTraceCoeffs[1];
            for(int i = 0; i < nvals; ++i)
            {
                field[m_traceCoeffsToElmtMap[0][i]] +=
                    m_traceCoeffsToElmtSign[0][i] * 
                    trace[m_traceCoeffsToElmtTrace[0][i]];
            }
        }


        void LocTraceToTraceMap::AddTraceCoeffsToFieldCoeffs(const int dir, const Array<OneD, const NekDouble> &trace, Array<OneD, NekDouble> &field)
        {
            int nvals = m_nTraceCoeffs[dir];
            for(int i = 0; i < nvals; ++i)
            {
                field[m_traceCoeffsToElmtMap[dir][i]] +=
                    m_traceCoeffsToElmtSign[dir][i] * 
                    trace[m_traceCoeffsToElmtTrace[dir][i]];
            }
        }
    } //namespace
} // namespace
