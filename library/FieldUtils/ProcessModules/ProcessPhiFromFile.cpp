///////////////////////////////////////////////////////////////////////////////
//
// File: ProcessPhiFromFile.cpp
//
// For more information, please see: http://www.nektar.info/
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
// Description: Reads an STL file.
//
///////////////////////////////////////////////////////////////////////////////

#include "ProcessPhiFromFile.h"
#include <boost/core/ignore_unused.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar;
using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessPhiFromFile::m_className = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "phifile"),
        ProcessPhiFromFile::create,
        "Computes the Phi function from a file, used in IB methods.")
};

/**
 * @brief Set up ProcessPhiFromFile object.
 *
 */
ProcessPhiFromFile::ProcessPhiFromFile(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scale"] = ConfigOption(false, "NotSet",
                                     "Scale coefficient for Phi.");
    m_config["file"] = ConfigOption(false, "NotSet",
                                    "File with the IB definition.");
}

/**
 *
 */
ProcessPhiFromFile::~ProcessPhiFromFile()
{
}

/**
 *
 */
void ProcessPhiFromFile::Process(po::variables_map &vm)
{
    // Ignore warnings due to 'vm'
    boost::ignore_unused(vm);

    // Do not run in parallel
    ASSERTL0(m_f->m_session->GetComm()->IsSerial(), "Parallel execution is "
                                                    "not supported yet in "
                                                    "this module.")

    // Check if required params are defined
    ASSERTL0(m_f->m_graph, "A session file file must be provided before the "
                           "STL file.")

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Read Phi function from the session file...
    if (m_config["file"].as<string>().compare("NotSet") == 0)
    {
        WARNINGL0(m_config["scale"].as<string>().compare("NotSet") == 0,
                  "Reading Phi function from session file, the provided scale "
                  "value will not be used");
        GetPhifromSession();
    }
    // ...or Read STL file and append Phi values to the existing expansions
    else
    {
        ASSERTL0(m_config["scale"].as<string>().compare("NotSet") != 0,
                 "Need to specify a scale coefficient, scale=value");

        STLobject phiFile = ReadSTL(m_config["file"].as<string>());
        GetPhifromSTL(phiFile);
    }
}

/**
 * @brief Read one 3D vector from a STL file, starting from the next line
 * of the input 'ifstream'. Numbers in ifstream are defined as 'float'
 *
 * @param in
 * @return Array<OneD, NekDouble>
 */
Array<OneD, NekDouble> ProcessPhiFromFile::ReadVector(ifstream &in)
{
    Array<OneD, NekDouble> out(3, 0.0);
    float tmp;
    char buf[4];

    in.read(buf, sizeof(buf));
    memcpy(&tmp, buf, sizeof(buf)); out[0] = tmp;
    in.read(buf, sizeof(buf));
    memcpy(&tmp, buf, sizeof(buf)); out[1] = tmp;
    in.read(buf, sizeof(buf));
    memcpy(&tmp, buf, sizeof(buf)); out[2] = tmp;

    return out;
}

/**
 * @brief Read an STL binary file and returns a struct of type 'STLobject'
 * containing the parsed data
 *
 * @param filename
 * @return ProcessPhiFromFile::STLobject
 */
ProcessPhiFromFile::STLobject ProcessPhiFromFile::ReadSTL(string filename)
{
    STLobject out;

    // Open file
    ifstream fileStl(filename.c_str(), ios::binary);
    ASSERTL0(fileStl, "An error occurred while trying to open the STL file.")

    // Buffers
    char headerBuf[80];
    char numTriBuf[4];
    char dumpBuf[2];

    // Read header and num of triangles
    fileStl.read(headerBuf, sizeof(headerBuf));
    fileStl.read(numTriBuf, sizeof(numTriBuf));
    memcpy(&out.numTri, numTriBuf, sizeof(numTriBuf));

    // Read triangle data
    out.triangles  = Array<OneD, triangle>(out.numTri);
    for (NekUInt32 i = 0; i < out.numTri; ++i)
    {
        // Read normal vector
        triangle tmpTri;
        tmpTri.normal = ReadVector(fileStl);

        // Read three vertices
        tmpTri.v0 = ReadVector(fileStl);
        tmpTri.v1 = ReadVector(fileStl);
        tmpTri.v2 = ReadVector(fileStl);

        // Add centroid to the triangle object
        Array<OneD, NekDouble> centre(3);
        centre[0] = (tmpTri.v0[0]+tmpTri.v1[0]+tmpTri.v2[0]) / 3.0;
        centre[1] = (tmpTri.v0[1]+tmpTri.v1[1]+tmpTri.v2[1]) / 3.0;
        centre[2] = (tmpTri.v0[2]+tmpTri.v1[2]+tmpTri.v2[2]) / 3.0;
        tmpTri.centroid = centre;

        out.triangles[i] = tmpTri;

        // Dump triangle type
        fileStl.read(dumpBuf, sizeof(dumpBuf));
    }

    // Close the file
    fileStl.close();

    return out;
}

/**
 * @brief Smoothing function for the SPM method given a distance value
 * and a scaling coefficient
 *
 * @param dist
 * @param coeff
 * @return NekDouble
 */
NekDouble ProcessPhiFromFile::PhiFunction(double dist, double coeff)
{
    return -0.5*(std::tanh(dist/coeff)-1.0);
}

/**
 * @brief Assigns to 'phi' the values indicated by 'ShapeFunction'
 *
 */
void ProcessPhiFromFile::GetPhifromSession()
{
    // Check that 'ShapeFunction' appears in the session file
    ASSERTL0(m_f->m_session->DefinesFunction("ShapeFunction"),
             "If file=file.stl is not supplied as an argument, a "
             "'ShapeFunction' block must be provided in the session "
             "file.")

    // Phi function in session file
    LibUtilities::EquationSharedPtr phiFunction =
        m_f->m_session->GetFunction("ShapeFunction", "Phi");

    // Get info about the domain
    int nPts  = m_f->m_exp[0]->GetNpoints();
    int nVars = m_f->m_variables.size();
    int nStrips;
    m_f->m_session->LoadParameter("Strip_Z", nStrips, 1);

    // Add new variable
    m_f->m_variables.push_back("phi");

    for (int s = 0; s < nStrips; ++s)
    {
        // Get current coords of the point
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        for (int i = 0; i < 3; ++i)
        {
            coords[i] = Array<OneD, NekDouble>(nPts, 0.0);
        }
        m_f->m_exp[s*nVars]->GetCoords(coords[0], coords[1], coords[2]);

        // Append Phi expansion to 'm_f'
        MultiRegions::ExpListSharedPtr Exp;
        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        phiFunction->Evaluate(coords[0], coords[1], coords[2],
                              Exp->UpdatePhys());
        Exp->FwdTrans(Exp->GetPhys(), Exp->UpdateCoeffs());

        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, Exp);
    }
}

/**
 * @brief Assigns to 'phi' the corresponding values of Phi
 *
 * @param file
 */
void ProcessPhiFromFile::GetPhifromSTL(
        const ProcessPhiFromFile::STLobject &file)
{
    // Get info about the domain
    int nPts  = m_f->m_exp[0]->GetNpoints();
    int nVars = m_f->m_variables.size();

    // Get coordinates
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    for (int i = 0; i < 3; ++i)
    {
        coords[i] = Array<OneD, NekDouble>(nPts, 0.0);
    }
    m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

    // Add new variable
    m_f->m_variables.push_back("phi");

    // Number of homogeneous strips
    int nStrips;
    m_f->m_session->LoadParameter("Strip_Z", nStrips, 1);

    // Find bounds of the mesh
    Array<OneD, NekDouble> bounds(6);
    bounds[0] = coords[0][0];
    bounds[1] = coords[0][0];
    bounds[2] = coords[1][0];
    bounds[3] = coords[1][0];
    bounds[4] = coords[2][0];
    bounds[5] = coords[2][0];
    for (int i = 1; i < nPts; ++i)
    {
        bounds[0] = (bounds[0] < coords[0][i]) ? bounds[0] : coords[0][i];
        bounds[1] = (bounds[1] > coords[0][i]) ? bounds[1] : coords[0][i];
        bounds[2] = (bounds[2] < coords[1][i]) ? bounds[2] : coords[1][i];
        bounds[3] = (bounds[3] > coords[1][i]) ? bounds[3] : coords[1][i];
        bounds[4] = (bounds[4] < coords[2][i]) ? bounds[4] : coords[2][i];
        bounds[5] = (bounds[5] > coords[2][i]) ? bounds[5] : coords[2][i];
    }
    // and add a margin to avoid rounding errors
    for (int i = 0; i < 6; ++i)
    {
        bounds[i] -= pow(-1,i) * 1e-10;
    }

    // Array of centroids of triangles in the STL object
    Array<OneD, Array<OneD, NekDouble> > centroids(file.numTri);
    for (int i = 0; i < file.numTri; ++i)
    {
        centroids[i] = file.triangles[i].centroid;
    }

    // Initialise octree
    m_tree = Octree(centroids, 10, bounds);

    // For each strip...
    for (int s = 0; s < nStrips; ++s)
    {
        // Append Phi expansion to 'm_f'
        MultiRegions::ExpListSharedPtr phi;
        phi = m_f->AppendExpList(m_f->m_numHomogeneousDir);

        // Parallelisation is highly recommended here
        for (int i = 0; i < nPts; ++i)
        {
            // Get coordinates of each point
            Array<OneD, NekDouble> tmpCoords(3);
            tmpCoords[2] = coords[2][i];
            tmpCoords[1] = coords[1][i];
            tmpCoords[0] = coords[0][i];

            // Find the shortest distance to the body(ies)
            double dist;
            FindShortestDist(file, tmpCoords, dist);

            // Get corresponding value of Phi
            phi->UpdatePhys()[i] = PhiFunction(dist, 
                                               m_config["scale"].as<double>());
        }

        // Update vector of expansions
        phi->FwdTrans(phi->GetPhys(), phi->UpdateCoeffs());
        auto it = m_f->m_exp.begin() + s * (nVars + 1) + nVars;
        m_f->m_exp.insert(it, phi);
    }
}

/**
 * @brief Checks if a ray traced from 'Origin' with direction 'Dvec' hits
 * the triangle defined by 'tri'. Returns the distance to the plane
 * defined by 'tri' in any case. A negative distance means that the hit
 * happened in the direction oposite that of the ray. Approach to calculate
 * the intersection point found in:
 *
 * Fast, minimum storage ray/triangle intersection,
 * Tomas Moeller, Ben Trumbore
 *
 * @param tri
 * @param Origin
 * @param Dvec
 * @param distance
 * @param u
 * @param v
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::CheckHit(const ProcessPhiFromFile::triangle &tri,
                                  const Array<OneD, NekDouble> &Origin,
                                  const Array<OneD, NekDouble> &Dvec,
                                  double &distance, double &u, double &v)
{
    // Edge vectors
    Array<OneD, NekDouble> E1(3);
    Array<OneD, NekDouble> E2(3);
    Vmath::Vsub(3, tri.v1, 1, tri.v0, 1, E1, 1);
    Vmath::Vsub(3, tri.v2, 1, tri.v0, 1, E2, 1);

    // If det == 0, ray parallel to triangle
    Array<OneD, NekDouble> Pvec = Cross(Dvec, E2);
    double det = Vmath::Dot(3, Pvec, E1);
    double inv_det = 1.0 / det;
    if (IsEqual(0.0, det, 1e-10))
    {
        distance = numeric_limits<double>::max();
        u        = numeric_limits<double>::max();
        v        = numeric_limits<double>::max();
        return false;
    }

    // Vector T and parameter u = (0.0, 1.0)
    Array<OneD, NekDouble> Tvec(3);
    Vmath::Vsub(3, Origin, 1, tri.v0, 1, Tvec, 1);
    u = Vmath::Dot(3, Pvec, Tvec) * inv_det;

    // Vector Q and parameter v = (0.0, 1.0)
    Array<OneD, NekDouble> Qvec = Cross(Tvec, E1);
    v = Vmath::Dot(3, Qvec, Dvec) * inv_det;

    // There is a hit if (u,v) coordinates are bounded
    distance = Vmath::Dot(3, Qvec, E2) * inv_det;
    if ((u < 0.0 || u > 1.0) || (v < 0.0 || u+v > 1.0))
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * @brief Calculates the shortest distance from a point \f$x\f$ to the closed
 * body contained in the STL file
 *
 * @param file
 * @param x
 * @param dist
 */
void ProcessPhiFromFile::FindShortestDist(
                                const ProcessPhiFromFile::STLobject &file,
                                const Array<OneD, NekDouble> &x,
                                double &dist)
{
    // Find closest triangles
    // First, use the ones in the node of 'x'
    int node = m_tree.QueryNode(x);

    // Set 'dist' to an unreal value
    dist = numeric_limits<double>::max();

    // If the node's depth is less than 3 the point is far from the object
    int depth = m_tree.QueryDepth(node);

    if (depth > 2)
    {
        vector<int> treeTriangles;
        vector<int> tmpTriangles = m_tree.QueryPoints(node);
        for (int point : tmpTriangles)
        {
            treeTriangles.push_back(point);
        }

        // Then, save the ones in the neighbouring nodes
        vector<int> neighbours = m_tree.QueryNeighbours(node);
        for (int neighNode : neighbours)
        {
            tmpTriangles = m_tree.QueryPoints(neighNode);
            for (int point : tmpTriangles)
            {
                treeTriangles.push_back(point);
            }
        }

        // Keep the sign (interior or exterior),
        int distSign = 1;
        // the normal vector of closest triangle so far,
        Array<OneD, NekDouble> triNormal =
            file.triangles[treeTriangles[0]].normal;
        // and the distance to the triangle PLANE
        double tParam = dist;

        for (int i = 0; i < treeTriangles.size(); ++i)
        {
            // Find distance to triangle
            triangle tri = file.triangles[treeTriangles[i]];
            double currentTparam;
            double u, v;
            bool hit = CheckHit(tri, x, tri.normal, currentTparam, u, v);

            // Save "sign" of 'currentTparam'
            int currentDistSign = (currentTparam >= 0) - (currentTparam < 0);
            // and distance to the triangle
            double currentDist;

            // Vector linking the hit point with the node
            Array<OneD, NekDouble> distVector(3);

            if (hit)
            {
                Vmath::Smul(3, -currentTparam, tri.normal, 1, distVector, 1);
                currentDist = fabs(currentTparam);
            }
            else
            {
                // The minimum has to be in one of the edges
                if (v < 0)   // Edge V0-V1
                {
                    distVector = Vector2edge(x, tri.v0, tri.v1);
                }
                else if (u < 0)   // Edge V0-V2
                {
                    distVector = Vector2edge(x, tri.v0, tri.v2);
                }
                else   // Edge V1-V2
                {
                    distVector = Vector2edge(x, tri.v1, tri.v2);
                }
                currentDist = sqrt(Vmath::Dot(3, distVector, distVector));
            }

            // Update 'dist', MAGIC CONSTANT AHEAD!
            // In the case of corners, check that the new triangle is not
            // giving contradictory information and, if so, use the one
            // closer to 'x'. Otherwise, some exterior points will be treated
            // as interior and viceversa
            if (dist-currentDist > 1e-5*currentDist ||
                (IsEqual(dist, currentDist, 1e-5) &&
                 IsNegative(Vmath::Dot(3, triNormal, tri.normal), 1e-5) &&
                 fabs(currentTparam) > fabs(tParam)))
            {
                dist      = currentDist;
                tParam    = currentTparam;
                distSign  = currentDistSign;
                triNormal = tri.normal;
            }
        }

        // Update distance sign
        dist *= -distSign;
    }
}

/**
 * @brief Returns true if \f$x=y\f$ within the relative tolerance 'relTol'
 * (relative to 'y')
 *
 * @param x
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::IsEqual(double x, double y, double relTol)
{
    return (fabs(x-y) <= relTol*y);
}

/**
 * @brief Returns true if \f$x<tol\f$
 *
 * @param x
 * @param relTol
 * @return true
 * @return false
 */
bool ProcessPhiFromFile::IsNegative(double x, double tol)
{
    return (x < tol);
}

/**
 * @brief Returns the cross product of vectors 'v0' y 'v1'
 *
 * @param v0
 * @param v1
 * @return Array<OneD, NekDouble>
 */
Array<OneD, NekDouble> ProcessPhiFromFile::Cross(
                                const Array<OneD, NekDouble> &v0,
                                const Array<OneD, NekDouble> &v1)
{
    Array<OneD, NekDouble> out(3);

    out[0] = v0[1]*v1[2] - v0[2]*v1[1];
    out[1] = v0[2]*v1[0] - v0[0]*v1[2];
    out[2] = v0[0]*v1[1] - v0[1]*v1[0];

    return out;
}

/**
 * @brief Determines the shortest distance from a point 'x' to the segment
 * defined by the points 'e1' and 'e2'. Note that this distance may be
 * equal to that to one of the end points. The vector returned points towards
 * the point 'x'
 *
 * @param x
 * @param e1
 * @param e2
 * @return Array<OneD, NekDouble>
 */
Array<OneD, NekDouble> ProcessPhiFromFile::Vector2edge(
                                        const Array<OneD, NekDouble> &x,
                                        const Array<OneD, NekDouble> &e1,
                                        const Array<OneD, NekDouble> &e2)
{
    size_t n = x.size();
    Array<OneD, NekDouble> e1x(n);
    Array<OneD, NekDouble> e1e2(n);
    Vmath::Vsub(n, x, 1, e1, 1, e1x, 1);
    Vmath::Vsub(n, e2, 1, e1, 1, e1e2, 1);
    double norm = sqrt(Vmath::Dot(n, e1e2, e1e2));
    for (size_t i = 0; i < n; ++i)
    {
        e1e2[i] /= norm;
    }

    Array<OneD, NekDouble> out(n);
    double proj = Vmath::Dot(n, e1x, e1e2);
    if (proj < 0.0)
    {
        Vmath::Vsub(n, x, 1, e1, 1, out, 1);
    }
    else if (proj > norm)
    {
        Vmath::Vsub(n, x, 1, e2, 1, out, 1);
    }
    else
    {
        Vmath::Svtsvtp(n, 1.0, e1x, 1, -proj, e1e2, 1, out, 1);
    }

    return out;
}
}
}
