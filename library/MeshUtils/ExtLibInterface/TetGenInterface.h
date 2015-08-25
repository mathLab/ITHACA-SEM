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

#ifndef NEKTAR_LIB_UTILITIES_MESHUTILS_TETGENINTERFACE_H
#define NEKTAR_LIB_UTILITIES_MESHUTILS_TETGENINTERFACE_H

#include <boost/shared_ptr.hpp>

#define TETLIBRARY
#include <tetgen.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <MeshUtils/MeshElem.hpp>

namespace Nektar{
namespace MeshUtils{

    class TetGenInterface
    {
    public:
        friend class MemoryManager<TetGenInterface>;

        TetGenInterface()
        {
        };

        void Assign(const std::vector<int> &nis,
                    const std::map<int, MeshTriSharedPtr> &t,
                    const std::map<int, MeshNodeSharedPtr> &n,
                    const std::vector<int> &stiner)
        {
            m_nodesinsurface = nis;
            Tris = t;
            Nodes = n;
            m_stienerpoints = stiner;
            meshloaded = false;
        }

        void Mesh(bool Quiet = true, bool Quality = false);

        void Extract(int &numtet, Array<OneD, Array<OneD, int> > &tetconnect);

    private:

        void freetet();

        tetgenio surface, additional, output;

        std::map<int, MeshTriSharedPtr> Tris;
        std::map<int, MeshNodeSharedPtr> Nodes;
        std::vector<int> m_stienerpoints;
        std::vector<int> m_nodesinsurface;

        std::map<int, int> nodemap;
        std::map<int, int> nodemapr;

        bool meshloaded;

    };

    typedef boost::shared_ptr<TetGenInterface> TetGenInterfaceSharedPtr;
}
}
#endif
