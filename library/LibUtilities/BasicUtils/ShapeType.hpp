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
#include <vector>

#include <boost/core/ignore_unused.hpp>

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


        namespace StdSegData
        {
            inline int getNumberOfCoefficients(int Na)
            {
                return Na;
            }

            inline int getNumberOfBndCoefficients(int Na)
            {
                boost::ignore_unused(Na);
                return 2;
            }
        }

        // Dimensions of coefficients for each space
        namespace StdTriData
        {
            inline int getNumberOfCoefficients(int Na, int Nb)
            {
                // Note these assertions have been set to > 0 because
                // it can also be used to evaluate face expansion
                // order
                ASSERTL2(Na > 0, "Order in 'a' direction must be > 0.");
                ASSERTL2(Nb > 0, "Order in 'b' direction must be > 0.");
                ASSERTL1(Na <= Nb, "order in 'a' direction is higher "
                         "than order in 'b' direction");
                return Na*(Na+1)/2 + Na*(Nb-Na);
            }

            inline int getNumberOfBndCoefficients(int Na, int Nb)
            {
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL1(Na <= Nb, "order in 'a' direction is higher "
                         "than order in 'b' direction");
                return (Na-1) + 2*(Nb-1);
            }
        }

        namespace StdQuadData
        {
            inline int getNumberOfCoefficients(int Na, int Nb)
            {
                // Note these assertions have been set to > 0 because
                // it can also be used to evaluate face expansion
                // order
                ASSERTL2(Na > 0, "Order in 'a' direction must be > 0.");
                ASSERTL2(Nb > 0, "Order in 'b' direction must be > 0.");
                return Na*Nb;
            }

            inline int getNumberOfBndCoefficients(int Na, int Nb)
            {
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                return 2*(Na-1) + 2*(Nb-1);
            }
        }



        namespace StdHexData
        {
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
                return Na*Nb*Nc;
            }

            inline int getNumberOfBndCoefficients(int Na, int Nb, int Nc)
            {
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
                return 2*Na*Nb + 2*Na*Nc + 2*Nb*Nc
                        - 4*(Na + Nb + Nc) + 8;
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
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "order in 'a' direction is higher "
                         "than order in 'c' direction");
                ASSERTL1(Nb <= Nc, "order in 'b' direction is higher "
                         "than order in 'c' direction");
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

            inline int getNumberOfBndCoefficients(int Na, int Nb, int Nc)
            {
                ASSERTL2(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL2(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL2(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "order in 'a' direction is higher "
                         "than order in 'c' direction");
                ASSERTL1(Nb <= Nc, "order in 'b' direction is higher "
                         "than order in 'c' direction");

                int nCoef =    Na*(Na+1)/2 + (Nb-Na)*Na // base
                          +    Na*(Na+1)/2 + (Nc-Na)*Na // front
                          + 2*(Nb*(Nb+1)/2 + (Nc-Nb)*Nb)// 2 other sides
                          - Na - 2*Nb - 3*Nc            // less edges
                          + 4;                          // plus vertices

                return nCoef;
            }
        }


        namespace StdPyrData
        {
            inline int getNumberOfCoefficients(int Na, int Nb, int Nc)
            {
                ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                        "than order in 'c' direction.");
                ASSERTL1(Nb <= Nc, "Order in 'b' direction is higher "
                        "than order in 'c' direction.");

                // Count number of coefficients explicitly.
                int nCoeff = 0;

                // Count number of interior tet modes
                for (int a = 0; a < Na; ++a)
                {
                    for (int b = 0; b < Nb; ++b)
                    {
                        for (int c = 0; c < Nc - std::max(a,b); ++c)
                        {
                            ++nCoeff;
                        }
                    }
                }
                return nCoeff;
            }

            inline int getNumberOfBndCoefficients(int Na, int Nb, int Nc)
            {
                ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                        "than order in 'c' direction.");
                ASSERTL1(Nb <= Nc, "Order in 'b' direction is higher "
                        "than order in 'c' direction.");

                return Na*Nb                        // base
                     + 2*(Na*(Na+1)/2 + (Nc-Na)*Na) // front and back
                     + 2*(Nb*(Nb+1)/2 + (Nc-Nb)*Nb) // sides
                     - 2*Na - 2*Nb - 4*Nc           // less edges
                     + 5;                           // plus vertices
            }
        }

        namespace StdPrismData
        {
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                        "than order in 'c' direction.");

                return Nb*StdTriData::getNumberOfCoefficients(Na,Nc);
            }

            inline int getNumberOfBndCoefficients( int Na, int Nb, int Nc)
            {
                ASSERTL1(Na > 1, "Order in 'a' direction must be > 1.");
                ASSERTL1(Nb > 1, "Order in 'b' direction must be > 1.");
                ASSERTL1(Nc > 1, "Order in 'c' direction must be > 1.");
                ASSERTL1(Na <= Nc, "Order in 'a' direction is higher "
                        "than order in 'c' direction.");

                return Na*Nb + 2*Nb*Nc              // rect faces
                   + 2*( Na*(Na+1)/2 + (Nc-Na)*Na ) // tri faces
                   - 2*Na - 3*Nb - 4*Nc             // less edges
                   + 6;                             // plus vertices
            }
        }

        inline int GetNumberOfCoefficients(ShapeType shape, std::vector<unsigned int> &modes, int offset)
        {
            int returnval = 0; 
            switch(shape)
            {
            case eSegment:
                returnval = modes[offset];
                break;
            case eTriangle:
                returnval = StdTriData::getNumberOfCoefficients(modes[offset],modes[offset+1]);
                break;
            case eQuadrilateral:
                returnval = modes[offset]*modes[offset+1];
                break;
            case eTetrahedron:
                returnval = StdTetData::getNumberOfCoefficients(modes[offset],modes[offset+1],modes[offset+2]);
                break;
            case ePyramid:
                returnval = StdPyrData::getNumberOfCoefficients(modes[offset],modes[offset+1],modes[offset+2]);
                break;
            case ePrism:
                returnval = StdPrismData::getNumberOfCoefficients(modes[offset],modes[offset+1],modes[offset+2]);
                break;
            case eHexahedron:
                returnval = modes[offset]*modes[offset+1]*modes[offset+2];
                break;
            default:
                NEKERROR(ErrorUtil::efatal,"Unknown Shape Type");
                break;
            }

            return returnval;
        }


        inline int GetNumberOfCoefficients(ShapeType shape, int na, int nb, int nc)
        {
            int returnval = 0; 
            switch(shape)
            {
            case eSegment:
                returnval = na;
                break;
            case eTriangle:
                returnval = StdTriData::getNumberOfCoefficients(na,nb);
                break;
            case eQuadrilateral:
                returnval = na*nb;
                break;
            case eTetrahedron:
                returnval = StdTetData::getNumberOfCoefficients(na,nb,nc);
                break;
            case ePyramid:
                returnval = StdPyrData::getNumberOfCoefficients(na,nb,nc);
                break;
            case ePrism:
                returnval = StdPrismData::getNumberOfCoefficients(na,nb,nc);
                break;
            case eHexahedron:
                returnval = na*nb*nc;
                break;
            default:
                NEKERROR(ErrorUtil::efatal,"Unknown Shape Type");
                break;
            }

            return returnval;
        }
    }
}

#endif

