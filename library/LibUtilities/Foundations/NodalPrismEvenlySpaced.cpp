///////////////////////////////////////////////////////////////////////////////
//
// File NodalPrismEvenlySpaced.cpp
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
// Description: 3D Nodal Prism Evenly Spaced Point Definitions
//
///////////////////////////////////////////////////////////////////////////////


#include <LibUtilities/Foundations/NodalPrismEvenlySpaced.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <vector>

namespace Nektar
{
    namespace LibUtilities
    {
        namespace
        {
            bool isVertex(int x, int y, int z, int npts)
            {
                return (x == 0        && y == 0        && z == 0       ) || 
                       (x == (npts-1) && y == 0        && z == 0       ) || 
                       (x == (npts-1) && y == (npts-1) && z == 0       ) || 
                       (x == 0        && y == (npts-1) && z == 0       ) || 
                       (x == 0        && y == 0        && z == (npts-1)) || 
                       (x == 0        && y == (npts-1) && z == (npts-1));
            }
            
            bool isEdge_01(int x, int y, int z, int npts)
            {
                return y == 0 && z == 0;
            }

            bool isEdge_12(int x, int y, int z, int npts)
            {
                return x == (npts-1) && z == 0;
            }

            bool isEdge_23(int x, int y, int z, int npts)
            {
                return y == (npts-1) && z == 0;
            }
            
            bool isEdge_30(int x, int y, int z, int npts)
            {
                return x == 0 && z == 0;
            }

            bool isEdge_04(int x, int y, int z, int npts)
            {
                return x == 0 && y == 0;
            }
            
            bool isEdge_14(int x, int y, int z, int npts)
            {
                return x + z == (npts-1) && y == 0;
            }
            
            bool isEdge_25(int x, int y, int z, int npts)
            {
                return x + z == (npts-1) && y == (npts-1);
            }
            
            bool isEdge_35(int x, int y, int z, int npts)
            {
                return x == 0 && y == (npts-1);
            }
            
            bool isEdge_45(int x, int y, int z, int npts)
            {
                return x == 0 && z == (npts-1);
            }

            bool isEdge(int x, int y, int z, int npts){
                return isEdge_01(x,y,z,npts) || isEdge_12(x,y,z,npts) || 
                       isEdge_23(x,y,z,npts) || isEdge_30(x,y,z,npts) ||
                       isEdge_04(x,y,z,npts) || isEdge_14(x,y,z,npts) ||
                       isEdge_25(x,y,z,npts) || isEdge_35(x,y,z,npts) ||
                       isEdge_45(x,y,z,npts);
            }

            bool isFace_0123(int x, int y, int z, int npts)
            {
                return z == 0;
            }

            bool isFace_014(int x, int y, int z, int npts)
            {
                return y == 0;
            }

            bool isFace_1254(int x, int y, int z, int npts)
            {
                return x + z == npts-1;
            }

            bool isFace_325(int x, int y, int z, int npts)
            {
                return y == (npts-1);
            }
            
            bool isFace_0354(int x, int y, int z, int npts)
            {
                return y == 0;
            }

            bool isFace(int x, int y, int z, int npts){
                return isFace_0123(x,y,z,npts) || isFace_014(x,y,z,npts) ||
                       isFace_1254(x,y,z,npts) || isFace_325(x,y,z,npts) ||
                       isFace_0354(x,y,z,npts);
            }
        }

        // Calculate evenly spaced number of points
        void NodalPrismEvenlySpaced::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();
            
            // Populate m_points
            unsigned int npts = GetNumPoints();
            NekDouble delta = 2.0/(npts - 1.0);
            for(unsigned int z=0, index=0; z<npts; ++z){
                for(int y=0; y<npts; ++y){
                    for(int x=0; x<npts-z; ++x, ++index){
                        NekDouble xi = -1.0 + x*delta;
                        NekDouble yi = -1.0 + y*delta;
                        NekDouble zi = -1.0 + z*delta;

                        m_points[0][index] = xi;
                        m_points[1][index] = yi;
                        m_points[2][index] = zi;
                    }
                }                
            }
            
            NodalPointReorder3d();            
        }

        void NodalPrismEvenlySpaced::NodalPointReorder3d()
        {
            unsigned int npts = GetNumPoints();
            using std::vector;
            vector<int> vertex;
            vector<int> iEdge_01;   // interior edge 0
            vector<int> iEdge_12;   // interior edge 1
            vector<int> iEdge_23;   // interior edge 2
            vector<int> iEdge_30;   // interior edge 3
            vector<int> iEdge_04;   // interior edge 4
            vector<int> iEdge_14;   // interior edge 5
            vector<int> iEdge_25;   // interior edge 6
            vector<int> iEdge_35;   // interior edge 7
            vector<int> iEdge_45;   // interior edge 8
            vector<int> iFace_0123; // interior face 0
            vector<int> iFace_014;  // interior face 1 
            vector<int> iFace_1254; // interior face 2
            vector<int> iFace_325;  // interior face 3
            vector<int> iFace_0354; // interior face 4
            vector<int> interiorVolumePoints; // interior volume points
            vector<int> map;

            // Build the lattice prism left to right - bottom to top
            for(int z=0, index=0; z<npts; ++z){
                for(int y=0; y<npts; ++y){
                    for(int x=0; x<npts-z; ++x, ++index){
                        if (isVertex(x,y,z,npts))
                        {
                            vertex.push_back(index);
                        } 
                        else if (isEdge(x,y,z,npts))
                        {
                            if (isEdge_01(x,y,z,npts))
                            {
                                iEdge_01.push_back(index);
                            }
                            else if (isEdge_12(x,y,z,npts))
                            {
                                iEdge_12.push_back(index);
                            } 
                            else if (isEdge_23(x,y,z,npts))
                            {
                                iEdge_23.push_back(index);
                            }
                            else if (isEdge_30(x,y,z,npts))
                            {
                                iEdge_30.push_back(index);
                            }
                            else if (isEdge_04(x,y,z,npts))
                            {
                                iEdge_04.push_back(index);
                            }
                            else if (isEdge_14(x,y,z,npts))
                            {
                                iEdge_14.push_back(index);
                            }
                            else if (isEdge_25(x,y,z,npts))
                            {
                                iEdge_25.push_back(index);
                            }
                            else if (isEdge_35(x,y,z,npts))
                            {
                                iEdge_35.push_back(index);
                            }
                            else if (isEdge_45(x,y,z,npts))
                            {
                                iEdge_45.push_back(index);
                            }
                        } 
                        else if (isFace(x,y,z,npts)) 
                        {
                            if (isFace_0123(x,y,z,npts))
                            {
                                iFace_0123.push_back(index);
                            } 
                            else if (isFace_014(x,y,z,npts))
                            {
                                iFace_014.push_back(index);
                            } 
                            else if (isFace_1254(x,y,z,npts))
                            {
                                iFace_1254.push_back(index);
                            } 
                            else if (isFace_325(x,y,z,npts))
                            {
                                iFace_325.push_back(index);
                            }
                            else if (isFace_0354(x,y,z,npts))
                            {
                                iFace_0354.push_back(index);
                            }
                        } 
                        else 
                        {
                            interiorVolumePoints.push_back(index);
                        }
                    }
                }
            }
           
            for (unsigned int n=0; n<vertex.size(); ++n)
            {
                map.push_back(vertex[n]);
            }
            
            for (unsigned int n=0; n<iEdge_01.size(); ++n)
            {
                map.push_back(iEdge_01[n]);
            }

            for (unsigned int n=0; n<iEdge_12.size(); ++n)
            {
                map.push_back(iEdge_12[n]);
            }

            for (unsigned int n=0; n<iEdge_23.size(); ++n)
            {
                map.push_back(iEdge_23[n]);
            }

            for (unsigned int n=0; n<iEdge_30.size(); ++n)
            {
                map.push_back(iEdge_30[n]);
            }

            for (unsigned int n=0; n<iEdge_04.size(); ++n)
            {
                map.push_back(iEdge_04[n]);
            }

            for (unsigned int n=0; n<iEdge_14.size(); ++n)
            {
                map.push_back(iEdge_14[n]);
            }

            for (unsigned int n=0; n<iEdge_25.size(); ++n)
            {
                map.push_back(iEdge_25[n]);
            }

            for (unsigned int n=0; n<iEdge_35.size(); ++n)
            {
                map.push_back(iEdge_35[n]);
            }

            for (unsigned int n=0; n<iEdge_45.size(); ++n)
            {
                map.push_back(iEdge_45[n]);
            }

            for (unsigned int n=0; n<iFace_0123.size(); ++n)
            {
                map.push_back(iFace_0123[n]);
            }

            for (unsigned int n=0; n<iFace_014.size(); ++n)
            {
                map.push_back(iFace_014[n]);
            }

            for(unsigned int n=0; n<iFace_1254.size(); ++n)
            {
                map.push_back(iFace_1254[n]);
            }

            for(unsigned int n=0; n<iFace_325.size(); ++n)
            {
                map.push_back(iFace_325[n]);
            }

            for(unsigned int n=0; n<iFace_0354.size(); ++n)
            {
                map.push_back(iFace_0354[n]);
            }

            for(unsigned int n=0; n<interiorVolumePoints.size(); ++n)
            {
                map.push_back(interiorVolumePoints[n]);
            }

            Array<OneD, NekDouble> points[3];
            points[0] = Array<OneD, NekDouble>(GetTotNumPoints());
            points[1] = Array<OneD, NekDouble>(GetTotNumPoints());
            points[2] = Array<OneD, NekDouble>(GetTotNumPoints());
            
            for(unsigned int index=0; index<map.size(); ++index)
            {
                points[0][index] = m_points[0][index];
                points[1][index] = m_points[1][index];
                points[2][index] = m_points[2][index];
            }

            for(unsigned int index=0; index<map.size(); ++index)
            {
                m_points[0][index] = points[0][map[index]];
                m_points[1][index] = points[1][map[index]];
                m_points[2][index] = points[2][map[index]];
            }
        }

        void NodalPrismEvenlySpaced::CalculateWeights()
        {            

        }


        // ////////////////////////////////////////
        //        CalculateInterpMatrix()
        void NodalPrismEvenlySpaced::CalculateInterpMatrix(const Array<OneD, const NekDouble>& xia,
                                                         const Array<OneD, const NekDouble>& yia,
                                                         const Array<OneD, const NekDouble>& zia,
                                                         Array<OneD, NekDouble>& interp)
        {
            ASSERTL0(false, "Not yet implemented");
        }

        // ////////////////////////////////////////
        //        CalculateDerivMatrix()
        void NodalPrismEvenlySpaced::CalculateDerivMatrix()
        {

        } 

        boost::shared_ptr<PointsBaseType> NodalPrismEvenlySpaced::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalPrismEvenlySpaced>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }
    } // end of namespace 
} // end of namespace 

