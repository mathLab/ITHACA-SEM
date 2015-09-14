////////////////////////////////////////////////////////////////////////////////
//
//  File: TetGenInterface.h
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
//  Description: class for interacting with tetgen
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_EXTLIBINTERFACE_TETGENINTERFACE_H
#define NEKTAR_MESHUTILS_EXTLIBINTERFACE_TETGENINTERFACE_H

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

        /**
         * @brief default constructor
         */
        TetGenInterface()
        {
        };

        /**
         * @brief assign parameters for meshing
         */
        void InitialMesh(const std::vector<int> &nis,
                         const std::vector<NekDouble> &ndel,
                         std::map<int, MeshTriSharedPtr> &Tris,
                         std::map<int, MeshNodeSharedPtr> &Nodes);

        void GetNewPoints(const std::vector<int> &nis, std::vector<Array<OneD, NekDouble> > &newp);

        void RefineMesh(const std::vector<int> &nis,
                        const std::vector<NekDouble> &ndel,
                        std::map<int, MeshTriSharedPtr> &Tris,
                        std::map<int, MeshNodeSharedPtr> &Nodes,
                        const std::vector<NekDouble> &newndel);
        /**
         * @brief extract mesh information
         */
        void Extract(int &numtet, Array<OneD, Array<OneD, int> > &tetconnect);

        void AddNodes(const std::vector<int> &nis, std::map<int, MeshNodeSharedPtr> &Nodes);

    private:

        /**
         * @brief clear previous mesh
         */
        void freetet();

        ///tetgen objects
        tetgenio surface, output, input;
        /// map from meshconvert id to tetgen id
        std::map<int, int> nodemap;
        /// map in reverse
        std::map<int, int> nodemapr;
        /// has a mesh been peviously loaded
        bool meshloaded;
};

typedef boost::shared_ptr<TetGenInterface> TetGenInterfaceSharedPtr;

}
}
#endif
