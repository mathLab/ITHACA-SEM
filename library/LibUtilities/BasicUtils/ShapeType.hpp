////////////////////////////////////////////////////////////////////////////////
//
//  File:  ShapeType.h
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
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHAPE_TYPE_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_SHAPE_TYPE_H

#include <algorithm>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#ifdef min
#undef min
#endif

namespace Nektar
{
    namespace LibUtilities
    {

        // Types of geometry types.
        enum ShapeType
        {
            eNoShapeType,
            ePoint,
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
            "NoGeomShapeType",
            "Point",
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
            0,  // ePoint
            1,  // eSegment
            2,  // eTriangle
            2,  // eQuadrilateral
            3,  // eTetrahedron
            3,  // ePyramid
            3,  // ePrism
            3,  // eHexahedron
        };


        // Dimensions of coefficients for each space
        namespace StdTriData
        {
            inline int getNumberOfCoefficients(int Na, int Nb)
            {
                ASSERTL0(Na <= Nb, "order in 'a' direction is higher "
                         "than order in 'b' direction");
                return Na*(Na+1)/2 + Na*(Nb-Na);
            }
        }


        namespace StdTetData
        {
            /**
             * Adds up the number of cells in a truncated Nc by Nc by Nc
             * pyramid, where the longest Na rows and longest Nb columns are
             * kept. Example: (Na, Nb, Nc) = (3, 4, 5); The number of
             * coefficients is the sum of the elements of the following
             * matrix:
             *
             * |5  4  3  2  0|
             * |4  3  2  0   |
             * |3  2  0      |
             * |0  0         |
             * |0            |
             * 
             * Sum = 28 = number of tet coefficients.
             */
            inline int getNumberOfCoefficients(int Na, int Nb, int Nc)
            {
                int nCoef = 0;
                for (int a = 0; a < Na; ++a)
                {
                    for (int b = 0; b < Nb - a; ++b)
                    {
                        for (int c = 0; c < Nc - a - b; ++c)
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }


        namespace StdPyrData
        {
            inline int getNumberOfCoefficients(int Na, int Nb, int Nc)
            {
                int nCoef = 0;
                for (int c = 0; c < Nc; ++c)
                {
                    for (int b = 0; b < std::min(Nc-c,Nb); ++b)
                    {
                        for (int a = 0 ; a < std::min(Nc-c,Na); ++a)
                        {
                            ++nCoef;
                        }
                    }
                }
                /*
                for (int a = 0; a < Na; ++a)
                {
                    for (int b = 0; b < Nb; ++b)
                    {
                        for (int c = 0; c < Nc - a - b; ++c)
                        {
                            ++nCoef;
                        }
                    }
                }
                */
                //std::cout << "Na = " << Na << " Nb = " << Nb << " Nc = " << Nc << " nCoef = " << nCoef << std::endl;
                return nCoef;
            }
        }

        namespace StdPrismData
        {
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for (int a = 0; a < Na; ++a)
                {
                    for (int b = 0; b < Nb; ++b)
                    {
                        for (int c = 0; c < Nc - a; ++c)
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

    }
}

#endif

