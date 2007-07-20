///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/PointExp.h,v $
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

#include <LocalRegions/LocalRegions.hpp>
#include <SpatialDomains/MeshComponents.h>

namespace Nektar
{
    namespace LocalRegions
    {

    class PointExp
    {
        
    public:
        PointExp(const SpatialDomains::VertexComponentSharedPtr &m_geom);
        ~PointExp(void);

            inline NekDouble GetValue(void)
            {
                return m_value;
            }

            inline void  SetValue(const NekDouble &value)
            {
                m_value = value;
            }

    protected:
            NekDouble      m_value; //!< Array containing expansion coefficients
        SpatialDomains::VertexComponentSharedPtr m_geom;
    private:

    };

    // type defines for use of PointExp in a boost vector
    typedef ptr<PointExp> PointExpSharedPtr;
    typedef std::vector<PointExpSharedPtr> PointExpVector;
    typedef std::vector<PointExpSharedPtr>::iterator PointExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // POINTEXP_H
