///////////////////////////////////////////////////////////////////////////////
//
// File StdRegions.hpp
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
// Description: Definition of enum lists and constants
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDREGIONS_HPP
#define STDREGIONS_HPP


// Headers from LibUtilities needed in StdRegions
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicUtils/Lapack.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>

#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>

#include <string>
#include <loki/Factory.h>

namespace Nektar
{
    namespace StdRegions
    {

        /** Number of workspace arrays for expansions and size of arrays;
        *  first value should be set to the maximum number of threads of
        *  execution used per process; second should be set to the number
        *  of workspaecs required for a single thread and the third value
        *  should be set to hold 3D point set per element
        *
        *  This data is used in StdWKSpace.(h/cpp)
        */
        namespace NekConstants
        {
            /** Tolerance to within which a point is considered to be located */
            const double kEvaluateTol  = 1e-12;

            // constants from NodalBasisManager.h & .cpp
            const int kMaxSym  = 5;
            const int kMaxBary = 4;
            const int kMaxDim  = 3;
            const int kBufSize = 1000;
        };

        enum GeomType
        {
            eRegular,
            eDeformed,
            eMovingRegular,
        };

        enum MatrixType
        {
            eMassMatrix,
            eLapMatrix,
            eNBasisTrans
        };

        const char* const MatrixTypeMap[] = 
	{
            "MassMatrix",
            "LapMatrix",
            "NBasisTrans"
        };

        /** enum list of StdExpansion regions */
        enum ShapeType
        {
            eNoShapeType,
            eSegment,
            eTriangle,
            eQuadrilateral,
            eTetrahedron,
            ePyramid,
            ePrism,
            eHexahedron,
            SIZE_ShapeType
        };
	

        const char* const ShapeTypeMap[] = 
        {
            "NoShapeType",
            "Segment",
            "Triangle",
            "Quadrilateral",
            "Tetrahedron",
            "Pyramid",
            "Prism",
            "Hexahedron"
        };
	



        // Hold the dimension of each of the types of shapes.
        const unsigned int ShapeTypeDimMap[SIZE_ShapeType] =
        {
            0,  // Unknown
            1,  // eSegment
            2,  // eTriangle
            2,  // eQuadrilateral
            3,  // eTetrahedron
            3,  // ePyramid
            3,  // ePrism
            3,  // eHexahedron
        };

        enum EdgeOrientation
        {
            eForwards,
            eBackwards
        };

        enum FaceOrientation
        {
            eDir1FwdDir1_Dir2FwdDir2,
            eDir1FwdDir1_Dir2BwdDir2,
            eDir1BwdDir1_Dir2FwdDir2,
            eDir1BwdDir1_Dir2BwdDir2,
            eDir1FwdDir2_Dir2FwdDir1,
            eDir1FwdDir2_Dir2BwdDir1,
            eDir1BwdDir2_Dir2FwdDir1,
            eDir1BwdDir2_Dir2BwdDir1
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

#endif //STDREGIONS_H

/**
* $Log: StdRegions.hpp,v $
* Revision 1.9  2007/02/23 19:26:08  jfrazier
* General bug fix and formatting.
*
* Revision 1.8  2007/02/22 22:02:28  sherwin
* Update with executing StdMatManager
*
* Revision 1.7  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.6  2007/02/17 03:40:20  jfrazier
* Couple changes to reflect additions and corrections to reflect linear algebra calls.
*
* Revision 1.5  2007/01/28 18:34:24  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.4  2007/01/20 22:35:21  sherwin
* Version with StdExpansion compiling
*
* Revision 1.3  2007/01/17 16:05:41  pvos
* updated doxygen documentation
*
* Revision 1.2  2006/06/01 13:43:20  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:33  kirby
* *** empty log message ***
*
* Revision 1.8  2006/03/23 21:55:02  jfrazier
* Declared NekConstants namespace to replace NekConstants class.
*
* Revision 1.7  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.6  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.5  2006/03/05 22:11:03  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.4  2006/03/04 20:26:55  bnelson
* Added comments after #endif.
*
* Revision 1.3  2006/02/26 23:37:30  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.2  2006/02/19 13:26:13  sherwin
*
* Coding standard revisions so that libraries compile
*
* Revision 1.1  2006/02/15 08:06:36  sherwin
*
* Put files into coding standard (although they do not compile)
*
*
**/
