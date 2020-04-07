///////////////////////////////////////////////////////////////////////////////
//
// File PointExp.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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

#include <LocalRegions/PointExp.h>

namespace Nektar
{
    namespace LocalRegions
    {

        PointExp::PointExp(const SpatialDomains::PointGeomSharedPtr &geom):
            StdExpansion  (1,0),
            StdExpansion0D(),
            StdRegions::StdPointExp(),
            Expansion     (geom),
            Expansion0D   (geom)
        {
            m_ncoeffs = 1;
        }

        PointExp::~PointExp(void)
        {
        }

        void PointExp::v_GetCoords(Array<OneD,NekDouble> &coords_0,
                Array<OneD, NekDouble> &coords_1,
                Array<OneD, NekDouble> &coords_2)
        {
            ASSERTL1(coords_0.size() > 0,
                     "Coords_0 is of insufficient size.");
            ASSERTL1(GetCoordim() < 2 || coords_1.size() > 0,
                     "Coords_1 is of insufficient size.");
            ASSERTL1(GetCoordim() < 3 || coords_2.size() > 0,
                     "Coords_2 is of insufficient size.");

            SpatialDomains::PointGeomSharedPtr v = std::dynamic_pointer_cast<SpatialDomains::PointGeom>(m_geom);
            NekDouble tmp;
            switch(GetCoordim())
            {
                case 1:
                    v->GetCoords(coords_0[0], tmp, tmp);
                    break;
                case 2:
                    v->GetCoords(coords_0[0], coords_1[0], tmp);
                    break;
                case 3:
                    v->GetCoords(coords_0[0], coords_1[0], coords_2[0]);
                    break;
            }
        }

        void PointExp::v_NormVectorIProductWRTBase(
            const Array<OneD, const NekDouble> &Fx,
                  Array<OneD,       NekDouble> &outarray)
        {
            const Array<OneD, const Array<OneD, NekDouble> >
                 &normals =
                    GetLeftAdjacentElementExp()->
                        GetVertexNormal(GetLeftAdjacentElementVertex());
            outarray[0] = Fx[0]*normals[0][0];
        }

    }
}
