////////////////////////////////////////////////////////////////////////////////
//
//  File: InputStarTec.h
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
//  Description: Tecplot (ascii .dat) converter.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_INPUTSTAR
#define UTILITIES_NEKMESH_INPUTSTAR

#include <NekMesh/Module/Module.h>

#include <libccmio/ccmio.h>

namespace Nektar
{
namespace NekMesh
{
/// Converter for VTK files.
class InputStar : public NekMesh::InputModule
{
public:
    /// Creates an instance of this class
    static NekMesh::ModuleSharedPtr create(NekMesh::MeshSharedPtr m)
    {
        return MemoryManager<InputStar>::AllocateSharedPtr(m);
    }
    static NekMesh::ModuleKey className;

    InputStar(NekMesh::MeshSharedPtr m);
    virtual ~InputStar();

    /// Populate and validate required data structures.
    virtual void Process();

    void ReadZone(int &nComposite);

    virtual std::string GetModuleName()
    {
        return "InputStar";
    }

protected:
    void GenElement3D(std::vector<NekMesh::NodeSharedPtr> &Nodes,
                      int i,
                      std::vector<int> &ElementFaces,
                      std::unordered_map<int, std::vector<int> > &FaceNodes,
                      int ncomposite,
                      bool DoOrient);

    void GenElement2D(std::vector<NekMesh::NodeSharedPtr> &Nodes,
                      int i,
                      std::vector<int> &FaceNodes,
                      int ncomposite);

    Array<OneD, int> SortEdgeNodes(
                    std::vector<NekMesh::NodeSharedPtr> &Nodes,
                    std::vector<int> &FaceNodes);

    Array<OneD, int> SortFaceNodes(
                    std::vector<NekMesh::NodeSharedPtr> &Nodes,
                    std::vector<int> &ElementFaces,
                    std::unordered_map<int, std::vector<int> > &FaceNodes);

    void ResetNodes(std::vector<NekMesh::NodeSharedPtr> &Nodes,
                    Array<OneD, std::vector<int> > &ElementFaces,
                    std::unordered_map<int, std::vector<int> > &FaceNodes);

private:
    CCMIOError m_ccmErr; // Star CCM error flag
    CCMIOID m_ccmTopology; // Star CCM mesh topology
    CCMIOID m_ccmProcessor;
    std::map<int, std::string> m_faceLabels; // label from CCM into composite

    void InitCCM(void);

    void ReadNodes(std::vector<NekMesh::NodeSharedPtr> &Nodes);

    void ReadInternalFaces(std::unordered_map<int, std::vector<int> > &FacesNodes,
                           Array<OneD, std::vector<int> > &ElementFaces);

    void ReadBoundaryFaces(std::vector<std::vector<int> > &BndElementFaces,
                           std::unordered_map<int, std::vector<int> > &FacesNodes,
                           Array<OneD, std::vector<int> > &ElementFaces,
                           std::vector<std::string> &facelabels);

    void SetupElements(void);
};
}
}

#endif
