////////////////////////////////////////////////////////////////////////////////
//
//  File: InputGmsh.h
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_INPUTGMSH
#define UTILITIES_NEKMESH_INPUTGMSH

#include <NekMeshUtils/Module/Module.h>

namespace Nektar
{
namespace Utilities
{

/**
 * Converter for Gmsh files.
 */
class InputGmsh : public NekMeshUtils::InputModule
{
public:
    InputGmsh(NekMeshUtils::MeshSharedPtr m);
    virtual ~InputGmsh();
    virtual void Process();

    /// Creates an instance of this class
    static NekMeshUtils::ModuleSharedPtr create(NekMeshUtils::MeshSharedPtr m)
    {
        return MemoryManager<InputGmsh>::AllocateSharedPtr(m);
    }
    /// %ModuleKey for class.
    static NekMeshUtils::ModuleKey className;
    static std::map<unsigned int, NekMeshUtils::ElmtConfig> GenElmMap();

    /**
     * Element map; takes a msh id to an %ElmtConfig object.
     */
    static std::map<unsigned int, NekMeshUtils::ElmtConfig> elmMap;
    static std::vector<int> CreateReordering(unsigned int InputGmshEntity);

private:
    int GetNnodes(unsigned int InputGmshEntity);
    static std::vector<int> TriReordering  (NekMeshUtils::ElmtConfig conf);
    static std::vector<int> QuadReordering (NekMeshUtils::ElmtConfig conf);
    static std::vector<int> HexReordering  (NekMeshUtils::ElmtConfig conf);
    static std::vector<int> PrismReordering(NekMeshUtils::ElmtConfig conf);
    static std::vector<int> TetReordering  (NekMeshUtils::ElmtConfig conf);
    static std::vector<int> LineReordering (NekMeshUtils::ElmtConfig conf);

    // Gmsh file version
    NekDouble m_version;
    // Previous id for contiguousness
    int m_prevId;
    // Id map if non-contiguous
    std::map<int, int> m_idMap;
    // Highest tag number
    int m_maxTagId;
    // This map takes each element ID and maps it to a permutation map
    // that is required to take Gmsh element node orderings and map them
    // to Nektar++ orderings.
    std::unordered_map<int, std::vector<int>> m_orderingMap;

    void ReadNextNode();
    void ReadNextNodeBlock(int nVertices = 0);
    void SaveNode(int id, NekDouble x = 0, NekDouble y = 0, NekDouble z = 0);
    void ReadNextElement(int tag = 0, int elm_type = 0);
};
}
}

#endif
