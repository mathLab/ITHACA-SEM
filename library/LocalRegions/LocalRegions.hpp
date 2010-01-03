///////////////////////////////////////////////////////////////////////////////
//
// File LocalRegions.hpp
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
// Description: General header for localregions containing enum and constants 
//
///////////////////////////////////////////////////////////////////////////////
#ifndef LOCALREGIONS_H
#define LOCALREGIONS_H

#include <string>

namespace Nektar
{    
    namespace LocalRegions
    {
        enum GeomState
        {
            eNotFilled,
            ePtsFilled,
        };

        // Defines a "fast find"
        // Assumes that first/last define the beginning/ending of
        // a continuous range of classes, and that start is
        // an iterator between first and last

        template<class InputIterator, class EqualityComparable>
        InputIterator find(InputIterator first, InputIterator last,
            InputIterator startingpoint,
            const EqualityComparable& value)
        {

            InputIterator val;

            if(startingpoint == first)
            {
                val = find(first,last,value);
            }
            else
            {
                val = find(startingpoint,last,value);
                if(val == last)
                {
                    val = find(first,startingpoint,value);
                    if(val == startingpoint)
                    {
                        val = last;
                    }
                }
            }
            return val;
        }

    } // end of namespace
} // end of namespace

#endif //LOCALREGIONS_H

/** 
*    $Log: LocalRegions.hpp,v $
*    Revision 1.6  2009/11/11 14:56:12  cantwell
*    Moved detailed connectivity documentation out of code.
*    Updated MultiRegions documentation.
*
*    Revision 1.5  2008/09/09 15:05:09  sherwin
*    Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
*
*    Revision 1.4  2008/06/05 15:06:42  pvos
*    Added documentation
*
*    Revision 1.3  2007/04/08 03:33:30  jfrazier
*    Minor reformatting and fixing SharedArray usage.
*
*    Revision 1.2  2006/05/06 20:36:16  sherwin
*    Modifications to get LocalRegions/Project1D working
*
*    Revision 1.1  2006/05/04 18:58:45  kirby
*    *** empty log message ***
*
*    Revision 1.2  2006/03/13 11:17:03  sherwin
*
*    First compiing version of Demos in SpatialDomains and LocalRegions. However they do not currently seem to execute properly
*
*    Revision 1.1  2006/03/12 07:43:31  sherwin
*
*    First revision to meet coding standard. Needs to be compiled
*
**/
