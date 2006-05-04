////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshGraph3D.h,v $
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
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H

#include <list>

namespace Nektar
{
    namespace SpatialDomains
    {
        class EdgeComponent;
        class TriFaceComponent;
        class QuadFaceComponent;
        class TetGeom;
        class PyrGeom;
        class PrimsGeom;
        class HexGeom;

        class MeshGraph3D: public MeshGraph
        {
        public:
            MeshGraph3D();
            virtual ~MeshGraph3D();

        protected:

        private:
            std::list<EdgeComponent*>     m_ecomps;
            std::list<TriFaceComponent*>  m_tfcomps;
            std::list<QuadFaceComponent*> m_qfcomps;
            std::list<TetGeom*>           m_tetgeoms;
            std::list<PyrGeom*>           m_pyrgeoms;
            std::list<PrismGeom*>         m_prismgeoms;
            std::list<HexGeom*>           m_hexgeoms;

            virtual void ReadPrivate(std::string &infilename);
            virtual void WritePrivate(std::string &outfilename);
        };
    }; // end of namespace
}; // end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH3D_H

//
// $Log: MeshGraph3D.h,v $
// Revision 1.6  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.5  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.4  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.3  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//
