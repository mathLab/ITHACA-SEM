///////////////////////////////////////////////////////////////////////////////
//
// File: NodalDemo.cpp
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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DEMOS_STDREGIONS_STDDEMOSUPPORT_HPP
#define DEMOS_STDREGIONS_STDDEMOSUPPORT_HPP

#include <string>
#include <vector>

#include <StdRegions/StdPointExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdNodalTetExp.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;
using namespace Nektar::StdRegions;
namespace po = boost::program_options;

class DemoSupport
{
public:
    DemoSupport() : m_desc("Available options")
    {
        m_desc.add_options()
            ("help,h",
             "Produce this help message and list basis and shape types.")
            ("nodal,n",
             po::value<string>(&m_ntype),
             "Optional nodal type, autofills shape and basis choices.")
            ("shape,s",
             po::value<string>(&m_shape),
             "Region shape to project function on.")
            ("basis,b",
             po::value<vector<string>>(&m_basis)->multitoken(),
             "Basis type, separate by spaces for higher dimensions.")
            ("order,o",
             po::value<vector<int>>(&m_order)->multitoken()->required(),
             "Order of basis sets, separate by spaces for higher dimensions.")
            ("points,p",
             po::value<vector<int>>(&m_points)->multitoken()->required(),
             "Number of quadrature points, separate by spaces for "
             "higher dimensions.")
            ("pointstype,P",
             po::value<vector<string>>(&m_pointstype)->multitoken(),
             "Optional points type, separate by spaces for higher dimensions.");
    }

    void ParseArguments(int argc, char *argv[])
    {
        try
        {
            po::store(po::parse_command_line(argc, argv, m_desc), m_vm);
            if (m_vm.count("help"))
            {
                cout << m_desc;
                cout << endl << "All nodal types, -n [ --nodal ], are:" << endl;
                for (int i = 22; i < SIZE_PointsType; ++i)
                {
                    cout << kPointsTypeStr[i] << endl;
                };
                cout << endl << "All shape types, -s [ --shape ], are:" << endl;
                for (int i = 1; i < SIZE_ShapeType; ++i)
                {
                    cout << ShapeTypeMap[i] << endl;
                };
                cout << endl << "All basis types, -b [ --basis ], are:" << endl;
                for (int i = 1; i < SIZE_BasisType; ++i)
                {
                    cout << BasisTypeMap[i] << endl;
                };
                cout << endl << "All points types, -P [ --pointstype ], are:"
                     << endl;
                for (int i = 1; i < SIZE_PointsType; ++i)
                {
                    cout << kPointsTypeStr[i] << endl;
                };
                exit(0);
            }
            po::notify(m_vm);
        }
        catch (const exception &e)
        {
            cerr << "Error: " << e.what() << endl << m_desc;
            exit(1);
        }
    }

    StdExpansion *CreateStdExpansion()
    {
        vector<PointsType> ptype;
        if (m_vm.count("pointstype"))
        {
            for (auto &p : m_pointstype)
            {
                PointsType tmp = eNoPointsType;

                // starts at nodal points
                for (int i = 1; i < SIZE_PointsType; ++i)
                {
                    if (boost::iequals(kPointsTypeStr[i], p))
                    {
                        tmp = static_cast<PointsType>(i);
                        break;
                    }
                    ASSERTL0(i != SIZE_PointsType - 1,
                             "The points type '" + p + "' does not exist");
                }

                ptype.push_back(tmp);
            }
        }

        // Convert string input argument to nodal type
        PointsType nodaltype = eNoPointsType;
        ShapeType stype = eNoShapeType;
        vector<BasisType> btype(3, eNoBasisType);
        if (m_vm.count("nodal"))
        {
            for (int i = 22; i < SIZE_PointsType; ++i) // starts at nodal points
            {
                if (boost::iequals(kPointsTypeStr[i], m_ntype))
                {
                    nodaltype = static_cast<PointsType>(i);
                    break;
                }
                ASSERTL0(i != SIZE_PointsType - 1,
                         "The nodal type '" + m_ntype + "' does not exist");
            }
            switch (nodaltype)
            {
                case eNodalTriElec:
                case eNodalTriFekete:
                case eNodalTriSPI:
                case eNodalTriEvenlySpaced:
                    btype[0] = eOrtho_A;
                    btype[1] = eOrtho_B;
                    stype = eTriangle;
                    break;
                case eNodalQuadElec:
                    btype[0] = eOrtho_A;
                    btype[1] = eOrtho_B;
                    stype = eQuadrilateral;
                    break;
                case eNodalTetElec:
                case eNodalTetSPI:
                case eNodalTetEvenlySpaced:
                    btype[0] = eOrtho_A;
                    btype[1] = eOrtho_B;
                    btype[2] = eOrtho_C;
                    stype = eTetrahedron;
                    break;
                case eNodalPrismElec:
                case eNodalPrismSPI:
                case eNodalPrismEvenlySpaced:
                    btype[0] = eOrtho_A;
                    btype[1] = eOrtho_A;
                    btype[2] = eOrtho_B;
                    stype = ePrism;
                    break;
                case eNodalHexElec:
                    btype[0] = eOrtho_A;
                    btype[1] = eOrtho_A;
                    btype[2] = eOrtho_A;
                    stype = eHexahedron;
                    break;
                default:
                    ASSERTL0(!nodaltype, ("The nodal type '" + m_ntype +
                                          "' is invalid for StdProject."));
                    break;
            }
        }

        //Convert string input argument to shape type
        if (stype == eNoShapeType)
        {
            for (int i = 1; i < SIZE_ShapeType; ++i)
            {
                if (boost::iequals(ShapeTypeMap[i], m_shape))
                {
                    stype = static_cast<ShapeType>(i);
                    break;
                }
                ASSERTL0(i != SIZE_ShapeType - 1,
                         "The shape type '" + m_shape + "' does not exist");
            }
        }

        // Check arguments supplied equals dimension
        const int dimension = (stype == ePoint) ? 1 : ShapeTypeDimMap[stype];
        ASSERTL0(m_order.size() == dimension,
                 "Number of orders supplied should match shape dimension");
        ASSERTL0(m_points.size() == dimension,
                 "Number of points supplied should match shape dimension");
        ASSERTL0(ptype.size() == dimension || ptype.size() == 0,
                 "Number of points types should match shape dimension if "
                 "supplied.");

        if (!m_vm.count("nodal"))
        {
            ASSERTL0(m_basis.size() == dimension,
                     "Number of bases supplied should match shape dimension");

            // Convert string input argument to basis types
            for (int i = 0; i < dimension; ++i)
            {
                for (int j = 1; j < SIZE_BasisType; ++j)
                {
                    if (boost::iequals(BasisTypeMap[j], m_basis[i]))
                    {
                        btype[i] = static_cast<BasisType>(j);
                        break;
                    }
                    ASSERTL0(j != SIZE_BasisType - 1,
                             "Basis type '" + m_basis[i] + "' does not exist");
                }
            }
        }

        //check basis selection is permitted for chosen shape
        map<ShapeType, vector<vector<BasisType>>> allowableBasis;
        allowableBasis[ePoint] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
             eLegendre, eChebyshev, eMonomial, eFourierSingleMode,
             eFourierHalfModeRe, eFourierHalfModeIm}
        };
        allowableBasis[eSegment] = { allowableBasis[ePoint][0] };
        allowableBasis[eTriangle] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
        };
        allowableBasis[eQuadrilateral] = {
            allowableBasis[eSegment][0], allowableBasis[eSegment][0]
        };
        allowableBasis[eTetrahedron] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_C, eModified_C, eGLL_Lagrange, eGauss_Lagrange}
        };
        allowableBasis[ePyramid] = {
            {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_A,    eModified_A,    eGLL_Lagrange, eGauss_Lagrange},
            {eOrthoPyr_C, eModifiedPyr_C, eGLL_Lagrange, eGauss_Lagrange}
        };
        allowableBasis[ePrism] = {
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_A, eModified_A, eGLL_Lagrange, eGauss_Lagrange},
            {eOrtho_B, eModified_B, eGLL_Lagrange, eGauss_Lagrange}
        };
        allowableBasis[eHexahedron] = {
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
             eLegendre, eChebyshev, eMonomial},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
             eLegendre, eChebyshev, eMonomial},
            {eOrtho_A, eModified_A, eFourier, eGLL_Lagrange, eGauss_Lagrange,
             eLegendre, eChebyshev, eMonomial}
        };


        for (int i = 0; i < dimension; ++i)
        {
            const unsigned int basisListLength = allowableBasis[stype][i].size();
            for (int j = 0; j < basisListLength; ++j)
            {
                if (allowableBasis[stype][i][j] == btype[i])
                {
                    break;
                }
                ASSERTL0(j != basisListLength - 1,
                         ("The basis type '" +
                          static_cast<string>(BasisTypeMap[btype[i]]) +
                          "' is invalid for basis argument " + to_string(i + 1) +
                          " for shape '" + ShapeTypeMap[stype] + "'."))
                    }
        }

        //Declaration of other variables needed
        StdExpansion *E = nullptr;

        // Assign points type according to basis type selection, if not already
        // assigned.
        if (ptype.size() == 0)
        {
            ptype.resize(dimension);
            for (int i = 0; i < dimension; ++i)
            {
                if (btype[i] == eFourier)
                {
                    ptype[i] = eFourierEvenlySpaced;
                }
                else if (btype[i] == eFourierSingleMode ||
                         btype[i] == eFourierHalfModeRe ||
                         btype[i] == eFourierHalfModeIm)
                {
                    ptype[i] = eFourierSingleModeSpaced;
                }
                else
                {
                    if (i == 1 && (stype == eTriangle || stype == eTetrahedron))
                    {
                        ptype[i] = eGaussRadauMAlpha1Beta0;
                    }
                    else if (i == 2 && (stype == eTetrahedron || stype == ePyramid))
                    {
                        ptype[i] = eGaussRadauMAlpha2Beta0;
                    }
                    else if (i == 2 && stype == ePrism)
                    {
                        ptype[i] = eGaussRadauMAlpha1Beta0;
                    }
                    else
                    {
                        ptype[i] = eGaussLobattoLegendre;
                    }
                }
            }
        }

        vector<PointsKey> pkey;
        vector<BasisKey> bkey;
        for (int i = 0; i < dimension; ++i)
        {
            pkey.emplace_back(PointsKey(m_points[i], ptype[i]));
            bkey.emplace_back(BasisKey(btype[i], m_order[i], pkey[i]));
        }

        switch (stype)
        {
            case ePoint:
            {
                E = new StdPointExp(bkey[0]);
                break;
            }
            case eSegment:
            {
                E = new StdSegExp(bkey[0]);
                break;
            }
            case eTriangle:
            {
                E = nodaltype != eNoPointsType ? new StdNodalTriExp(bkey[0],
                                                                    bkey[1],
                                                                    nodaltype)
                    : new StdTriExp(bkey[0], bkey[1]);
                break;
            }
            case eQuadrilateral:
            {
                E = new StdQuadExp(bkey[0], bkey[1]);
                break;
            }
            case eTetrahedron:
            {
                E = nodaltype != eNoPointsType ? new StdNodalTetExp(bkey[0],
                                                                    bkey[1],
                                                                    bkey[2],
                                                                    nodaltype)
                    : new StdTetExp(bkey[0], bkey[1],
                                    bkey[2]);
                break;
            }
            case ePyramid:
            {
                E = new StdPyrExp(bkey[0], bkey[1], bkey[2]);
                break;
            }
            case ePrism:
            {
                E = nodaltype != eNoPointsType ? new StdNodalPrismExp(bkey[0],
                                                                      bkey[1],
                                                                      bkey[2],
                                                                      nodaltype)
                    : new StdPrismExp(bkey[0], bkey[1],
                                      bkey[2]);
                break;
            }
            case eHexahedron:
            {
                E = new StdHexExp(bkey[0], bkey[1], bkey[2]);
                break;
            }
            default:
                break;
        }

        return E;
    }

    po::options_description &GetOptions()
    {
        return m_desc;
    }

    po::variables_map &GetVariableMap()
    {
        return m_vm;
    }

    std::vector<string> &GetPointsType()
    {
        return m_pointstype;
    }

    Array<OneD, Array<OneD, NekDouble>> GetCoords(StdExpansion *E)
    {
        int dimension = E->GetShapeDimension();
        const auto totPoints = (unsigned) E->GetTotPoints();

        Array<OneD, NekDouble> x(totPoints), y(totPoints), z(totPoints);
        Array<OneD, Array<OneD, NekDouble>> coords(dimension);

        switch (dimension)
        {
            case 1:
            {
                E->GetCoords(x);
                coords[0] = x;
                break;
            }

            case 2:
            {
                E->GetCoords(x, y);
                coords[0] = x;
                coords[1] = y;
                break;
            }

            case 3:
            {
                E->GetCoords(x, y, z);
                coords[0] = x;
                coords[1] = y;
                coords[2] = z;
                break;
            }
            default:
                break;
        }

        return coords;
    }

protected:
    po::options_description m_desc;
    po::variables_map m_vm;

    std::string    m_shape;
    std::string    m_ntype;
    vector<string> m_basis{3, "NoBasisType"};
    vector<string> m_pointstype{3, "NoPointsType"};
    vector<int>    m_order;
    vector<int>    m_points;
};

#endif
