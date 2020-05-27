////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCurve.cpp
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
//  Description: create curved edges using a custom expression y = f(x)
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Interpolator.h>

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessCurve.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessCurve::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "curve"), 
    ProcessCurve::create,
    "Creates curved edge information using a function y=f(x) (2D only).");

/**
 * @brief Default constructor.
 */
ProcessCurve::ProcessCurve(MeshSharedPtr m) : ProcessCurvedEdges(m)
{
    m_config["function"] = ConfigOption(false, "NotSet",
        "Expression of the curve: y = f(x).");
    m_config["file"] = ConfigOption(false, "NotSet",
        "Pts file containing coordinates (x,y).");
    m_config["niter"] = ConfigOption(false, "50",
        "Number of iterations to perform to obtain evenly distribution of points.");
    m_config["gamma"] = ConfigOption(false, "0.1",
        "Relaxation parameter.");
}

/**
 * @brief Destructor.
 */
ProcessCurve::~ProcessCurve()
{
}

/**
 * This function generates curved edge information in 2D, by including
 * equispaced nodes to an edge following the prescribed function \f$y = f(x)\f$.
 *
 * An approximately equispaced node distribution is obtained by initially
 * considering a uniform distribution in the \f$ x\f$ direction, and then
 * adjusting the interior nodes using the iteration
 *  \f[ x^{n}_{i+1} = x^{n-1}_{i+1} + \Delta x^{n}_{i+1}   \f]
 * \f[ \Delta x^{n}_{i+1} = \Delta x^{n}_{i} +
 *    \gamma (\overline{\Delta s} - \Delta s^{n}_{i}) sign(\Delta x^{n}_{i}) \f]
 *
 * where \f$ x^{n}\f$ is the \f$ x\f$ coordinate of the \f$ n^{th}\f$ node,
 *  \f$ i\f$ is the iteration counter, \f$ \Delta s^{n}\f$ is the Cartesian
 * distance between nodes \f$ n\f$ and \f$ n-1\f$ and \f$ \overline{\Delta s}\f$
 * is the average of these distances. The relaxation factor \f$ \gamma\f$ and
 * the fixed number of iterations \f$ Niter\f$ are parameters of the module.
 *
 * In case the correction to \f$ \Delta x\f$ leads to an invalid distribution
 * (i.e. reversing directions), the relaxation parameter is successively halved
 * until it leads to a valid distribution.
 *
 * @param edge Edge which will be modified
 */
void ProcessCurve::v_GenerateEdgeNodes(EdgeSharedPtr edge)
{
    NodeSharedPtr n1 = edge->m_n1;
    NodeSharedPtr n2 = edge->m_n2;

    int       nq    = m_config["N"].as<int>();
    int       niter = m_config["niter"].as<int>();
    NekDouble gamma = m_config["gamma"].as<double>();

    edge->m_edgeNodes.resize(nq - 2);

    // Read function defining the curve
    if (m_config["function"].as<string>().compare("NotSet") != 0)
    {
        ASSERTL0(m_config["file"].as<string>().compare("NotSet") == 0,
             "Function and file cannot be defined at the same time.");

        std::string fstr = m_config["function"].as<string>();
        m_fExprId = m_fEval.DefineFunction("x y z", fstr);
        m_fromFile = false;
    }
    else
    {
        ASSERTL0(m_config["file"].as<string>().compare("NotSet") != 0,
             "Need to define either function or file.");
        std::string inFile = m_config["file"].as<string>();

        LibUtilities::CommSharedPtr  c =
                LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
        LibUtilities::PtsIOSharedPtr ptsIO =
                MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(c);
        ptsIO->Import(inFile, m_fieldPts);

        m_fromFile = true;
    }

    // Coordinates of points
    Array<OneD, NekDouble> x(nq,0.0);
    Array<OneD, NekDouble> y(nq,0.0);
    // Distances between points
    Array<OneD, NekDouble> dx(nq-1,0.0);
    Array<OneD, NekDouble> dy(nq-1,0.0);
    Array<OneD, NekDouble> ds(nq-1,0.0);
    // Average distance and deviation from average
    Array<OneD, NekDouble> s_deviation(nq-1,0.0);
    NekDouble s_average;

    // Fix start point
    x[0] =  n1->m_x;
    y[0] =  EvaluateCoordinate(x[0]);
    // Start with uniform distribution along x-axis
    Vmath::Sadd(nq-1, (n2->m_x - n1->m_x) / (nq-1), dx, 1, dx, 1);

    // Iterate a few times to make points more evenly distributed
    for (int s = 0; s < niter; ++s)
    {
        s_average = 0.0;
        for (int k = 1; k < nq; ++k)
        {
            x[k] = x[k-1] + dx[k-1] + gamma*s_deviation[k-1];
            y[k] =  EvaluateCoordinate(x[k]);

            dx[k-1] = x[k] - x[k-1];
            dy[k-1] = y[k] - y[k-1];
            ds[k-1] = sqrt(dx[k-1]*dx[k-1] + dy[k-1]*dy[k-1]);

            s_average = s_average + ds[k-1]/(nq-1);
        }
        // Calculate correction for next iteration
        for (int k = 0; k < nq-1; ++k)
        {
            s_deviation[k] =  (s_average - ds[k])* 
                                ( dx[k]/abs(dx[k]));
        }
        // Adjust gama to make sure next partition is valid
        // (no reversals in dx)
        bool valid = false;
        gamma = m_config["gamma"].as<double>();
        while( !valid)
        {
            valid = true;
            for (int k = 0; k < nq-1; ++k)
            {
                if( dx[k]*(dx[k]+gamma*s_deviation[k]) < 0.0)
                {
                    gamma = gamma/2;
                    valid = false;
                    continue;
                }
            }
        }
    }

    // Write interior nodes to edge
    for (int k = 1; k < nq-1; ++k)
    {
        edge->m_edgeNodes[k-1] = NodeSharedPtr(
                    new Node(0, x[k], y[k], n1->m_z));
    }
    edge->m_curveType = LibUtilities::ePolyEvenlySpaced;
}

NekDouble ProcessCurve::EvaluateCoordinate(NekDouble xCoord)
{
    if (m_fromFile)
    {
        Array<OneD, Array<OneD, NekDouble> > tmp(2);
        tmp[0] = Array<OneD, NekDouble>(1, xCoord);
        tmp[1] = Array<OneD, NekDouble>(1, 0.0);

        LibUtilities::PtsFieldSharedPtr toPts =
            MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(1, tmp);

        LibUtilities::Interpolator interp;
        interp.Interpolate(m_fieldPts, toPts);

        return tmp[1][0];
    }
    else
    {
        return m_fEval.Evaluate(m_fExprId, xCoord, 0.0, 0.0, 0.0);
    }
}

}
}
