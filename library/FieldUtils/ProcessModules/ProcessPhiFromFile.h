///////////////////////////////////////////////////////////////////////////////
//
// File: ProcessPhiFromFile.h
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
// Description: geometry file converter.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSPHIFROMFILE
#define FIELDUTILS_PROCESSPHIFROMFILE

#include "../Module.h"
#include "../Octree.h"
#include <random>

namespace Nektar
{
namespace FieldUtils
{

/**
 * Converter for STL files.
 */
class ProcessPhiFromFile : public ProcessModule
{
public:
    /// Creates an instance of this class
    static ModuleSharedPtr create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessPhiFromFile>::AllocateSharedPtr(f);
    }
    /// ModuleKey for class.
    static ModuleKey m_className;
    
    ProcessPhiFromFile(FieldSharedPtr f);
    virtual ~ProcessPhiFromFile();
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessPhiFromFile";
    }

    virtual std::string GetModuleDescription()
    {
        return "Processing input STL file to calculate Phi";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eModifyExp;
    }

protected:
    /// Object representing a 3D triangle
    struct triangle
    {
        Array<OneD, NekDouble> normal;
        Array<OneD, NekDouble> v0;
        Array<OneD, NekDouble> v1;
        Array<OneD, NekDouble> v2;
        Array<OneD, NekDouble> centroid;

        triangle() : normal(3, 0.0), v0(3, 0.0), v1(3, 0.0), v2(3, 0.0) {}
    };
    /// STL file object
    struct STLobject
    {
        // Number of triangles
        NekUInt32 numTri;
        // Triangles definition
        Array<OneD, triangle> triangles;
    };
    /// Octree object
    Octree m_tree;

    // Reads one vector from a binary STL file
    Array<OneD, NekDouble> ReadVector(std::ifstream &in);
    // Reads an STL file and returns an 'STLobject' struct
    STLobject ReadSTL(std::string filename);
    // Smoothing function
    NekDouble PhiFunction(double dist, double coeff);
    // Calculates Phi from 'ShapeFunction' in the session file
    void GetPhifromSession();
    // Calculates Phi from an external STL binary file
    void GetPhifromSTL(const STLobject &file);
    // Checks if a ray hits a specific 3D triangle
    bool CheckHit(const triangle &tri,
                  const Array<OneD, NekDouble> &Origin,
                  const Array<OneD, NekDouble> &Dvec,
                  double &distance, double &u, double &v);
    // Shortest distance from a point to a 3D geometry
    void FindShortestDist(const STLobject &file,
                          const Array<OneD, NekDouble> &x,
                          double &dist);
    // Utility to find if two doubles are equal
    bool IsEqual(double x, double y, double relTol);
    // Utility to find if a double is negative with a certain tolerance
    bool IsNegative(double x, double tol);
    // Utility to calculate the cross-product of two 3D vectors
    Array<OneD, NekDouble> Cross(const Array<OneD, NekDouble> &v0,
                                 const Array<OneD, NekDouble> &v1);
    // Utility to measure shortest distance to a segment
    Array<OneD, NekDouble> Vector2edge(const Array<OneD, NekDouble> &x,
                                       const Array<OneD, NekDouble> &e1,
                                       const Array<OneD, NekDouble> &e2);

private:
};
}
}

#endif   // FIELDUTILS_PROCESSPHIFROMFILE
