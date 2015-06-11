////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>

#include <LibUtilities/MeshUtils/CurveMesh.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{
namespace MeshUtils {
    
    void CurveMesh::Mesh()
    {
        m_cadcurve->Bounds(m_bounds);
        m_curvelength = m_cadcurve->Length(m_bounds[0],m_bounds[1]);
        
        m_numSamplePoints = int(m_curvelength/m_octree->GetMinDelta())+1;
        
        GetSampleFunction();
        
        Ae = 0.0;
        
        for(int i = 0; i < m_numSamplePoints-1; i++)
        {
            Ae+=m_dst[1][1]*(1.0/m_dst[i][0]+1.0/m_dst[i+1][0])/2.0;
        }
        
        Ne=round(Ae);
        
        GetPhiFunction();
        
        meshsvalue = Array<OneD, NekDouble>(Ne+1);
        meshsvalue[0]=0.0;
        meshsvalue[Ne]=m_curvelength;
        
        for(int i = 1; i < Ne; i++)
        {
            bool iterate = true;
            int k = i;
            NekDouble ski = meshsvalue[i-1];
            NekDouble lastSki;
            while(iterate)
            {
                NekDouble rhs =EvaluateDS(ski)/Ae*(EvaluatePS(ski)-k);
                lastSki=ski;
                ski=ski-rhs;
                if(abs(lastSki-ski)<1E-10)
                {
                    iterate = false;
                }
            }
            
            meshsvalue[i]=ski;
        }
        
        m_meshpoints = Array<OneD, Array<OneD, NekDouble> >(Ne+1);
        
        for(int i = 0; i < Ne; i++)
        {
            NekDouble t = m_cadcurve->tAtArcLength(meshsvalue[i]);
            Array<OneD, NekDouble> loc;
            m_cadcurve->P(t,loc);
            m_meshpoints[i]=loc;
        }
        NekDouble t = m_bounds[1];
        Array<OneD, NekDouble> loc;
        m_cadcurve->P(t,loc);
        m_meshpoints[Ne]=loc;
    }
    
    void CurveMesh::GetPhiFunction()
    {
        m_ps = Array<OneD, Array<OneD, NekDouble> >(m_numSamplePoints);
        Array<OneD, NekDouble> newPhi(2);
        
        newPhi[0]=0.0;
        newPhi[1]=0.0;
        
        m_ps[0]=newPhi;
        
        NekDouble runningInt=0.0;
        
        for(int i = 1; i < m_numSamplePoints; i++)
        {
            runningInt+=(1.0/m_dst[i-1][0]+1.0/m_dst[i][0])/2.0*m_dst[1][1];
            newPhi[0]=Ne/Ae*runningInt;
            newPhi[1]=m_dst[i][1];
            m_ps[i]=newPhi;
        }
        
    }
    
    NekDouble CurveMesh::EvaluateDS(NekDouble s)
    {
        int a=0;
        int b=0;
        
        if(s==0)
        {
            return m_dst[0][0];
        }
        else if(s==m_curvelength)
        {
            return m_dst[m_numSamplePoints-1][0];
        }
        
        for(int i = 0; i < m_numSamplePoints-2; i++)
        {
            if(m_dst[i][1]<s && m_dst[i+1][1]>=s)
            {
                a=i;
                b=i+1;
                break;
            }
        }
        
        NekDouble s1 = m_dst[a][1];
        NekDouble s2 = m_dst[b][1];
        NekDouble d1 = m_dst[a][0];
        NekDouble d2 = m_dst[b][0];
        
        NekDouble m  = (d2-d1)/(s2-s1);
        NekDouble c = d2- m*s2;
        
        return m*s+c;
    }
    
    NekDouble CurveMesh::EvaluatePS(NekDouble s)
    {
        int a=0;
        int b=0;
        
        if(s==0)
        {
            return m_ps[0][0];
        }
        else if(s==m_curvelength)
        {
            return m_ps[m_numSamplePoints-1][0];
        }
        
        for(int i = 0; i < m_numSamplePoints-2; i++)
        {
            if(m_ps[i][1]<s && m_ps[i+1][1]>=s)
            {
                a=i;
                b=i+1;
                break;
            }
        }
        
        NekDouble s1 = m_ps[a][1];
        NekDouble s2 = m_ps[b][1];
        NekDouble d1 = m_ps[a][0];
        NekDouble d2 = m_ps[b][0];
        
        NekDouble m  = (d2-d1)/(s2-s1);
        NekDouble c = d2- m*s2;
        
        return m*s+c;
    }
    
    void CurveMesh::GetSampleFunction()
    {
        m_dst = Array<OneD, Array<OneD, NekDouble> >(m_numSamplePoints);
        Array<OneD, NekDouble> loc;
        m_cadcurve->P(m_bounds[0],loc);
        
        Array<OneD, NekDouble> dsti(3);
        
        dsti[0]=m_octree->Query(loc);
        dsti[1]=0.0;
        dsti[2]=m_bounds[0];
        
        m_dst[0] = dsti;
        
        NekDouble targetLength = m_curvelength/(m_numSamplePoints-1);
        
        for(int i = 1; i < m_numSamplePoints; i++)
        {
            dsti[1]=i*targetLength;
            NekDouble t = m_cadcurve->tAtArcLength(dsti[1]);
            
            m_cadcurve->P(t,loc);
            
            dsti[0]=m_octree->Query(loc);
            dsti[2]=t;
            
            m_dst[i]=dsti;
        }
    }
    
    
            
}
}
}

