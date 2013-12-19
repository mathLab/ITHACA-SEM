///////////////////////////////////////////////////////////////////////////////
//
// File PointExp.h
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
// Description: Definition of a Point expansion 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POINTEXP_H
#define POINTEXP_H

#include <StdRegions/StdPointExp.h>
#include <SpatialDomains/PointGeom.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/Expansion0D.h>

namespace Nektar
{
    namespace LocalRegions
    {
        class PointExp: virtual public StdRegions::StdPointExp, virtual public Expansion0D
        {
        public:
            LOCAL_REGIONS_EXPORT PointExp(const SpatialDomains::PointGeomSharedPtr &m_geom);
            LOCAL_REGIONS_EXPORT ~PointExp(void);

        protected:
            virtual void v_GetCoords(
                Array<OneD, NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
            {
                Array<OneD, NekDouble> coords(3);
                SpatialDomains::PointGeomSharedPtr v =
                    boost::dynamic_pointer_cast<SpatialDomains::PointGeom>(m_geom);
                v->GetCoords(coords);
                coords_0[0] = coords[0];
                coords_1[0] = coords[1];
                coords_2[0] = coords[2];
            }
        };
        
        // type defines for use of PointExp in a boost vector
        typedef boost::shared_ptr<PointExp> PointExpSharedPtr;
        typedef std::vector<PointExpSharedPtr> PointExpVector;
        typedef std::vector<PointExpSharedPtr>::iterator PointExpVectorIter;
        
        const static Array<OneD, PointExpSharedPtr> NullPointExpSharedPtrArray;
    } //end of namespace
} //end of namespace

#endif // POINTEXP_H
