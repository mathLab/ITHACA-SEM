///////////////////////////////////////////////////////////////////////////////
//
// File NodalPrismElec.cpp
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
// Description: 3D Nodal Prism eletrostatic Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/NodalPrismElec.h>
#include <LibUtilities/Foundations/Points.h>
#include <vector>

namespace Nektar
{
namespace LibUtilities
{

bool NodalPrismElec::initPointsManager[] = {
    PointsManager().RegisterCreator(PointsKey(0, eNodalPrismElec),         NodalPrismElec::Create)
};

namespace
{
bool isVertex(int t, int y, int npts)
{
    return (t == 0 && y == 0) || (t == 1 && y == 0) || (t == 2 && y == 0) ||
           (t == 0 && y == (npts - 1)) || (t == 1 && y == (npts - 1)) ||
           (t == 2 && y == (npts - 1));
}

bool isEdge_01(int t, int y, int npts)
{
    return y == 0 && t > 2 && t <= npts;
}

bool isEdge_12(int t, int y, int npts)
{
    boost::ignore_unused(y, npts);
    return t == 1;
}

bool isEdge_23(int t, int y, int npts)
{
    return y == (npts - 1) && t > 2 && t <= npts;
}

bool isEdge_30(int t, int y, int npts)
{
    boost::ignore_unused(y, npts);
    return t == 0;
}

bool isEdge_04(int t, int y, int npts)
{
    return y == 0 && t >= 3 + 2 * (npts - 2) && t < 3 + 3 * (npts - 2);
}

bool isEdge_14(int t, int y, int npts)
{
    return y == 0 && t >= 3 + (npts - 2) && t < 3 + 2 * (npts - 2);
}

bool isEdge_25(int t, int y, int npts)
{
    return y == npts - 1 && t >= 3 + (npts - 2) && t < 3 + 2 * (npts - 2);
}

bool isEdge_35(int t, int y, int npts)
{
    return y == npts - 1 && t >= 3 + 2 * (npts - 2) && t < 3 + 3 * (npts - 2);
}

bool isEdge_45(int t, int y, int npts)
{
    boost::ignore_unused(y, npts);
    return t == 2;
}

bool isEdge(int t, int y, int npts)
{
    return isEdge_01(t, y, npts) || isEdge_12(t, y, npts) ||
           isEdge_23(t, y, npts) || isEdge_30(t, y, npts) ||
           isEdge_04(t, y, npts) || isEdge_14(t, y, npts) ||
           isEdge_25(t, y, npts) || isEdge_35(t, y, npts) ||
           isEdge_45(t, y, npts);
}

bool isFace_0123(int t, int y, int npts)
{
    boost::ignore_unused(y);
    return t < 3 + (npts - 2);
}

bool isFace_014(int t, int y, int npts)
{
    boost::ignore_unused(t, npts);
    return y == 0;
}

bool isFace_1254(int t, int y, int npts)
{
    boost::ignore_unused(y);
    return t < 3 + 2 * (npts - 2) && t >= 3 + (npts - 2);
}

bool isFace_325(int t, int y, int npts)
{
    boost::ignore_unused(t);
    return y == (npts - 1);
}

bool isFace_0354(int t, int y, int npts)
{
    boost::ignore_unused(y);
    return t < 3 + 3 * (npts - 2) && t >= 3 + 2 * (npts - 2);
}

bool isFace(int t, int y, int npts)
{
    return isFace_0123(t, y, npts) || isFace_014(t, y, npts) ||
           isFace_1254(t, y, npts) || isFace_325(t, y, npts) ||
           isFace_0354(t, y, npts);
}
}

// Calculate evenly spaced number of points
void NodalPrismElec::CalculatePoints()
{
    // Allocate the storage for points
    PointsBaseType::CalculatePoints();

    // Populate m_points
    unsigned int npts = GetNumPoints();

    LibUtilities::PointsKey pkey1(npts, LibUtilities::eNodalTriElec);
    Array<OneD, NekDouble> u1, v1;
    LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
    LibUtilities::PointsKey pkey2(npts, LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> u;
    LibUtilities::PointsManager()[pkey2]->GetPoints(u);

    for (unsigned int y = 0, index = 0; y < npts; y++)
    {
        for (size_t t = 0; t < u1.size(); t++, index++)
        {
            m_points[0][index] = u1[t];
            m_points[1][index] = u[y];
            m_points[2][index] = v1[t];
        }
    }

    NodalPointReorder3d();
    m_util = MemoryManager<NodalUtilPrism>::AllocateSharedPtr(
        npts - 1, m_points[0], m_points[1], m_points[2]);
}

void NodalPrismElec::NodalPointReorder3d()
{
    unsigned int npts = GetNumPoints();
    using std::vector;
    vector<int> vertex;
    vector<int> iEdge_01;             // interior edge 0
    vector<int> iEdge_12;             // interior edge 1
    vector<int> iEdge_23;             // interior edge 2
    vector<int> iEdge_30;             // interior edge 3
    vector<int> iEdge_04;             // interior edge 4
    vector<int> iEdge_14;             // interior edge 5
    vector<int> iEdge_25;             // interior edge 6
    vector<int> iEdge_35;             // interior edge 7
    vector<int> iEdge_45;             // interior edge 8
    vector<int> iFace_0123;           // interior face 0
    vector<int> iFace_014;            // interior face 1
    vector<int> iFace_1254;           // interior face 2
    vector<int> iFace_325;            // interior face 3
    vector<int> iFace_0354;           // interior face 4
    vector<int> interiorVolumePoints; // interior volume points
    vector<int> map;

    // Build the lattice prism left to right - bottom to top
    for (unsigned int y = 0, index = 0; y < npts; y++)
    {
        for (unsigned int t = 0; t < npts * (npts + 1) / 2; t++, index++)
        {
            if (isVertex(t, y, npts))
            {
                vertex.push_back(index);
            }
            else if (isEdge(t, y, npts))
            {
                if (isEdge_01(t, y, npts))
                {
                    iEdge_01.push_back(index);
                }
                else if (isEdge_12(t, y, npts))
                {
                    iEdge_12.push_back(index);
                }
                else if (isEdge_23(t, y, npts))
                {
                    iEdge_23.push_back(index);
                }
                else if (isEdge_30(t, y, npts))
                {
                    iEdge_30.push_back(index);
                }
                else if (isEdge_04(t, y, npts))
                {
                    iEdge_04.push_back(index);
                }
                else if (isEdge_14(t, y, npts))
                {
                    iEdge_14.push_back(index);
                }
                else if (isEdge_25(t, y, npts))
                {
                    iEdge_25.push_back(index);
                }
                else if (isEdge_35(t, y, npts))
                {
                    iEdge_35.push_back(index);
                }
                else if (isEdge_45(t, y, npts))
                {
                    iEdge_45.push_back(index);
                }
            }
            else if (isFace(t, y, npts))
            {
                if (isFace_0123(t, y, npts))
                {
                    iFace_0123.push_back(index);
                }
                else if (isFace_014(t, y, npts))
                {
                    iFace_014.push_back(index);
                }
                else if (isFace_1254(t, y, npts))
                {
                    iFace_1254.push_back(index);
                }
                else if (isFace_325(t, y, npts))
                {
                    iFace_325.push_back(index);
                }
                else if (isFace_0354(t, y, npts))
                {
                    iFace_0354.push_back(index);
                }
            }
            else
            {
                interiorVolumePoints.push_back(index);
            }
        }
    }

    // sort vertices
    std::swap(vertex[2], vertex[4]);
    // sort edges
    std::reverse(iEdge_23.begin(), iEdge_23.end());
    std::reverse(iEdge_30.begin(), iEdge_30.end());
    std::reverse(iEdge_04.begin(), iEdge_04.end());
    std::reverse(iEdge_35.begin(), iEdge_35.end());

    // faces
    for (unsigned int i = 0; i < npts - 2; i++)
    {
        for (unsigned int j = i + 1; j < npts - 2; j++)
        {
            std::swap(iFace_1254[i * (npts - 2) + j],
                      iFace_1254[j * (npts - 2) + i]);
        }
    }
    for (int i = 0; i < npts - 2; i++)
    {
        std::reverse(iFace_0354.begin() + (i * (npts - 2)),
                     iFace_0354.begin() + (i * (npts - 2) + npts - 2));
    }
    for (unsigned int i = 0; i < npts - 2; i++)
    {
        for (unsigned int j = i + 1; j < npts - 2; j++)
        {
            std::swap(iFace_0354[i * (npts - 2) + j],
                      iFace_0354[j * (npts - 2) + i]);
        }
    }

    for (unsigned int n = 0; n < vertex.size(); ++n)
    {
        map.push_back(vertex[n]);
    }

    for (unsigned int n = 0; n < iEdge_01.size(); ++n)
    {
        map.push_back(iEdge_01[n]);
    }

    for (unsigned int n = 0; n < iEdge_12.size(); ++n)
    {
        map.push_back(iEdge_12[n]);
    }

    for (unsigned int n = 0; n < iEdge_23.size(); ++n)
    {
        map.push_back(iEdge_23[n]);
    }

    for (unsigned int n = 0; n < iEdge_30.size(); ++n)
    {
        map.push_back(iEdge_30[n]);
    }

    for (unsigned int n = 0; n < iEdge_04.size(); ++n)
    {
        map.push_back(iEdge_04[n]);
    }

    for (unsigned int n = 0; n < iEdge_14.size(); ++n)
    {
        map.push_back(iEdge_14[n]);
    }

    for (unsigned int n = 0; n < iEdge_25.size(); ++n)
    {
        map.push_back(iEdge_25[n]);
    }

    for (unsigned int n = 0; n < iEdge_35.size(); ++n)
    {
        map.push_back(iEdge_35[n]);
    }

    for (unsigned int n = 0; n < iEdge_45.size(); ++n)
    {
        map.push_back(iEdge_45[n]);
    }

    for (unsigned int n = 0; n < iFace_0123.size(); ++n)
    {
        map.push_back(iFace_0123[n]);
    }

    for (unsigned int n = 0; n < iFace_014.size(); ++n)
    {
        map.push_back(iFace_014[n]);
    }

    for (unsigned int n = 0; n < iFace_1254.size(); ++n)
    {
        map.push_back(iFace_1254[n]);
    }

    for (unsigned int n = 0; n < iFace_325.size(); ++n)
    {
        map.push_back(iFace_325[n]);
    }

    for (unsigned int n = 0; n < iFace_0354.size(); ++n)
    {
        map.push_back(iFace_0354[n]);
    }

    for (unsigned int n = 0; n < interiorVolumePoints.size(); ++n)
    {
        map.push_back(interiorVolumePoints[n]);
    }

    Array<OneD, NekDouble> points[3];
    points[0] = Array<OneD, NekDouble>(GetTotNumPoints());
    points[1] = Array<OneD, NekDouble>(GetTotNumPoints());
    points[2] = Array<OneD, NekDouble>(GetTotNumPoints());

    for (unsigned int index = 0; index < map.size(); ++index)
    {
        points[0][index] = m_points[0][index];
        points[1][index] = m_points[1][index];
        points[2][index] = m_points[2][index];
    }

    for (unsigned int index = 0; index < map.size(); ++index)
    {
        m_points[0][index] = points[0][map[index]];
        m_points[1][index] = points[1][map[index]];
        m_points[2][index] = points[2][map[index]];
    }
}

void NodalPrismElec::CalculateWeights()
{
    // Allocate the storage for points
    PointsBaseType::CalculateWeights();

    typedef DataType T;

    // Solve the Vandermonde system of integrals for the weight vector
    NekVector<T> w = m_util->GetWeights();

    m_weights = Array<OneD, T>(w.GetRows(), w.GetPtr());
}

// ////////////////////////////////////////
//        CalculateInterpMatrix()
void NodalPrismElec::CalculateInterpMatrix(
    const Array<OneD, const NekDouble> &xia,
    const Array<OneD, const NekDouble> &yia,
    const Array<OneD, const NekDouble> &zia, Array<OneD, NekDouble> &interp)
{
    Array<OneD, Array<OneD, NekDouble> > xi(3);
    xi[0] = xia;
    xi[1] = yia;
    xi[1] = zia;

    std::shared_ptr<NekMatrix<NekDouble> > mat =
        m_util->GetInterpolationMatrix(xi);
    Vmath::Vcopy(mat->GetRows() * mat->GetColumns(), mat->GetRawPtr(), 1,
                 &interp[0], 1);
}

// ////////////////////////////////////////
//        CalculateDerivMatrix()
void NodalPrismElec::CalculateDerivMatrix()
{
    // Allocate the derivative matrix.
    PointsBaseType::CalculateDerivMatrix();

    m_derivmatrix[0] = m_util->GetDerivMatrix(0);
    m_derivmatrix[1] = m_util->GetDerivMatrix(1);
    m_derivmatrix[2] = m_util->GetDerivMatrix(2);
}

std::shared_ptr<PointsBaseType> NodalPrismElec::Create(const PointsKey &key)
{
    std::shared_ptr<PointsBaseType> returnval(
        MemoryManager<NodalPrismElec>::AllocateSharedPtr(key));

    returnval->Initialize();

    return returnval;
}

}
}
