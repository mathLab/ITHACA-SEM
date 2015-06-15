////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.h
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


#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_TRIANGLEINTERFACE_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_TRIANGLEINTERFACE_H

#include <boost/shared_ptr.hpp>

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#define ANSI_DECLARATORS
#define TRILIBRARY
#define VOID int

#include <LibUtilities/MeshUtils/ExternalLib/triangle.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar {
    namespace LibUtilities {
        namespace MeshUtils {
            
            class TriangleInterface
            {
            public:
                friend class MemoryManager<TriangleInterface>;
                
                TriangleInterface(const std::vector<std::vector<
                                  std::vector<NekDouble> > > &boundingloops,
                                  const std::vector<std::vector<NekDouble> >
                                  &centers,
                                  const std::vector<std::vector<NekDouble> >
                                  &stienerp) : m_boundingloops(boundingloops),
                                  m_centers(centers), m_stienerpoints(stienerp)
                {
                    meshloaded = false;
                };
                TriangleInterface()
                {
                    meshloaded = false;
                };
                
                void Assign(const std::vector<std::vector<
                            std::vector<NekDouble> > > &boundingloops,
                            const std::vector<std::vector<NekDouble> >
                            &centers,
                            const std::vector<std::vector<NekDouble> >
                            &stienerp)
                {
                    if(meshloaded)
                    {
                        freetri();
                    }
                    m_boundingloops = boundingloops;
                    m_centers = centers;
                    m_stienerpoints = stienerp;
                }
                
                ~TriangleInterface()
                {
                    if(meshloaded)
                        freetri();
                }
                
                void Mesh(bool Quiet = true, bool Quality = false);
                
                void Extract(int &np,
                             int &nt,
                             Array<OneD, Array<OneD, NekDouble> > &Points,
                             Array<OneD, Array<OneD, int> > &Connec);
                
            private:
                
                void freetri();
                
                std::vector<std::vector<std::vector<NekDouble> > >
                        m_boundingloops;
                std::vector<std::vector<NekDouble> > m_centers;
                std::vector<std::vector<NekDouble> > m_stienerpoints;
                
                bool meshloaded;
                
                struct triangulateio in,out;
                
            };
            
            typedef boost::shared_ptr<TriangleInterface>
                            TriangleInterfaceSharedPtr;
        }
    }
}

#endif
