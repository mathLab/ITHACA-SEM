////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.h
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_PROCESSVAROPTI
#define UTILITIES_NEKMESH_PROCESSVAROPTI

#include "../../Module.h"

namespace Nektar
{
namespace Utilities
{

struct ElData
{
    ElData(vector<NodeSharedPtr> innodes, int spaceDim)
    {
        nodes.resize(innodes.size());
        for (int i = 0; i < innodes.size(); ++i)
        {
            nodes[i].resize(spaceDim);
            nodes[i][0] = &innodes[i]->m_x;

            if (spaceDim >= 2)
            {
                nodes[i][1] = &innodes[i]->m_y;
            }

            if (spaceDim >= 3)
            {
                nodes[i][2] = &innodes[i]->m_z;
            }
        }
    }

    ElementSharedPtr el;
    vector<vector<NekDouble *> > nodes;
    vector<Array<OneD, NekDouble> > maps;
};
typedef boost::shared_ptr<ElData> ElDataSharedPtr;

struct DerivUtil
{
    NekMatrix<NekDouble> VdmD[3];
    NekVector<NekDouble> quadW;
};
typedef boost::shared_ptr<DerivUtil> DerivUtilSharedPtr;

struct PtsHelper
{
    int ptsHigh;
    int ptsLow;
};
typedef boost::shared_ptr<PtsHelper> PtsHelperSharedPtr;

enum optimiser
{
    eLinEl,
    eWins,
    eRoca,
    eHypEl
};

struct Residual
{
    NekDouble val;
    int n;
    int nDirc;
    int nDoF;
    int startInv;
    NekDouble worstJac;
};

typedef boost::shared_ptr<Residual> ResidualSharedPtr;

/**
 * @brief This processing module calculates the Jacobian of elements
 * using %SpatialDomains::GeomFactors and the %Element::GetGeom
 * method. For now it simply prints a list of elements which have
 * negative Jacobian.
 */
class ProcessVarOpti : public ProcessModule
{
public:
    /// Creates an instance of this class
    static boost::shared_ptr<Module> create(MeshSharedPtr m)
    {
        return MemoryManager<ProcessVarOpti>::AllocateSharedPtr(m);
    }
    static ModuleKey className;

    ProcessVarOpti(MeshSharedPtr m);
    virtual ~ProcessVarOpti();

    /// Write mesh to output file.
    virtual void Process();
private:
    typedef map<int, vector<ElDataSharedPtr> > NodeElMap;

    void FillQuadPoints();
    void GetElementMap();
    vector<Array<OneD, NekDouble> > MappingIdealToRef(ElementSharedPtr el);
    vector<vector<NodeSharedPtr> > GetColouredNodes();
    void WriteStats(string file);
    void EvaluateMesh();

    NodeElMap nodeElMap;
    vector<ElDataSharedPtr> dataSet;

    ResidualSharedPtr res;
    DerivUtilSharedPtr derivUtil;
    PtsHelperSharedPtr ptsHelp;
    optimiser opti;
};

}
}

#endif
