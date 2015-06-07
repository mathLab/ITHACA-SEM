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


#include <LibUtilities/MeshUtils/Octree.h>

using namespace std;
namespace Nektar{
    namespace LibUtilities{
        namespace MeshUtils {
            
            void Octree::Build(const NekDouble &m_minDelta,
                               const NekDouble &m_maxDelta,
                               const NekDouble &m_eps)
            {
                BoundingBox = Array<OneD, NekDouble> (6);
                m_cad->GetBoundingBox(BoundingBox);
                
                CompileCuravturePointList(m_minDelta,m_maxDelta,m_eps);
                
                /*octant newOctant((BoundingBox[1]+BoundingBox[0])/2,
                                 (BoundingBox[3]+BoundingBox[2])/2,
                                 (BoundingBox[5]+BoundingBox[4])/2,
                                 (BoundingBox[1]-BoundingBox[0])/2,
                                 (BoundingBox[3]-BoundingBox[2])/2,
                                 (BoundingBox[5]-BoundingBox[4])/2,
                                 NULL, 0, P, curvaturePointList, octantList);
                
                octantList.push_back(newOctant);
                //parent created.
                
                cout << endl << "Parent created. Dividing based on geometry" << endl;
                
                if(octantList[0].Divide())
                {
                    octantList[0].LeafFalse();
                    subdivide(0, P, curvaturePointList, octantList);
                }
                
                int ct=0;
                
                for(int i = 0; i < octantList.size(); i++)
                {
                    if(octantList[i].isLeaf()){ct++;}
                    if(octantList[i].GetLevel()>P.maxLevel)
                    {
                        P.maxLevel=octantList[i].GetLevel();
                    }
                }
                
                cout << endl << "No. octant leaves:: Max level" << endl;
                cout << ct << " " << P.maxLevel << endl;*/
            }
            
            
            void Octree::CompileCuravturePointList(const NekDouble &m_minDelta,
                                                   const NekDouble &m_maxDelta,
                                                   const NekDouble &m_eps)
            {
                NekDouble MaxDim = 0.0;
                if(BoundingBox[1]-BoundingBox[0]>MaxDim)
                    MaxDim = BoundingBox[1]-BoundingBox[0];
                if(BoundingBox[3]-BoundingBox[2]>MaxDim)
                    MaxDim = BoundingBox[3]-BoundingBox[2];
                if(BoundingBox[5]-BoundingBox[4]>MaxDim)
                    MaxDim = BoundingBox[5]-BoundingBox[4];
                
                
                
                int ns = MaxDim/m_minDelta;
                
                for(int i = 1; i <= m_cad->GetNumSurf(); i++)
                {
                    Array<OneD, NekDouble> ParameterPlaneBounds;
                    m_cad->GetParameterPlaneBounds(i,ParameterPlaneBounds);
                    
                    NekDouble du = (ParameterPlaneBounds[1]-ParameterPlaneBounds[0])
                                    /(ns-1);
                    NekDouble dv = (ParameterPlaneBounds[3]-ParameterPlaneBounds[2])
                                        /(ns-1);
                    
                    for(int j = 0; j < ns; j++)
                    {
                        for(int k = 0; k < ns; k++)
                        {
                            NekDouble u = du*j + ParameterPlaneBounds[0];
                            NekDouble v = dv*k + ParameterPlaneBounds[2];
                            
                            Array<OneD, NekDouble> N;
                            Array<OneD, NekDouble> r;
                            
                            m_cad->N(i,u,v,N);
                            m_cad->D2(i,u,v,r);
                            
                            if(N[0]==0 && N[1]==0 && N[2]==0)
                            {
                                continue;
                            }
                            
                            NekDouble E = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
                            NekDouble F = r[3]*r[6] + r[4]*r[7] + r[5]*r[8];
                            NekDouble G = r[6]*r[6] + r[7]*r[7] + r[8]*r[8];
                            NekDouble e = N[3]*r[9] + N[4]*r[10] + N[5]*r[11];
                            NekDouble f = N[3]*r[15] + N[4]*r[16] + N[5]*r[17];
                            NekDouble g = N[3]*r[12] + N[4]*r[13] + N[5]*r[14];
                            
                            if(E*G-F*F<1E-30)
                            {
                                continue;
                            }
                            
                            NekDouble K, H;
                            
                            K = (e*g-f*f)/(E*G-F*F);
                            H = 0.5*(e*G-2*f*F+g*E)/(E*G-F*F);
                            
                            NekDouble kv[2];
                            kv[0] = abs(H + sqrt(H*H-K));
                            kv[1] = abs(H - sqrt(H*H-K));
                            
                            if(kv[0] != 0 || kv[1] != 0)
                            {
                                CurvaturePointSharedPtr newCPoint =
                                MemoryManager<CurvaturePoint>::AllocateSharedPtr
                                (r[0],r[1],r[2],
                                 1.0/(kv[0] > kv[1] ? kv[0] : kv[1]),
                                 N[0],N[1],N[2]);
                                
                                m_cpList.push_back(newCPoint);
                            }else
                            {
                                CurvaturePointSharedPtr newCPoint =
                                MemoryManager<CurvaturePoint>::AllocateSharedPtr
                                (r[0],r[1],r[2],
                                 N[0],N[1],N[2]);
                                m_cpList.push_back(newCPoint);
                            }
                        }
                    }
                }
                
                for(int i = 0; i < m_cpList.size(); i++)
                {
                    m_cpList[i]->Process(m_minDelta,m_maxDelta,m_eps);
                }
                
            }
            
        }
    }
}
