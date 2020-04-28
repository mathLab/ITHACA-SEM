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

#include "ProcessVarOpti.h"

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
namespace Utilities
{

map<LibUtilities::ShapeType, DerivUtilSharedPtr> ProcessVarOpti::BuildDerivUtil(
    int o)
{
    // build Vandermonde information
    map<LibUtilities::ShapeType, DerivUtilSharedPtr> ret;

    // Typedef for points types used in the variational optimser. First entry is
    // the evaluation points; second entry is the distributions used for the
    // full element (includes element boundary)
    typedef std::pair<LibUtilities::PointsType, LibUtilities::PointsType>
        PTypes;

    map<LibUtilities::ShapeType, PTypes> typeMap;

    if (m_mesh->m_nummode + o <= 11)
    {
        typeMap[LibUtilities::eTriangle] =
            PTypes(LibUtilities::eNodalTriSPI, LibUtilities::eNodalTriElec);
        typeMap[LibUtilities::eTetrahedron] =
            PTypes(LibUtilities::eNodalTetSPI, LibUtilities::eNodalTetElec);
        typeMap[LibUtilities::ePrism] =
            PTypes(LibUtilities::eNodalPrismSPI, LibUtilities::eNodalPrismElec);
    }

    typeMap[LibUtilities::eQuadrilateral] =
        PTypes(LibUtilities::eNodalQuadElec, LibUtilities::eNodalQuadElec);
    // typeMap[LibUtilities::eHexahedron] =
    //    PTypes(LibUtilities::eNodalHexElec, LibUtilities::eNodalHexElec);

    for (auto &it : typeMap)
    {
        PTypes pType           = it.second;
        DerivUtilSharedPtr der = std::shared_ptr<DerivUtil>(new DerivUtil());

        LibUtilities::PointsKey pkey1(m_mesh->m_nummode, pType.second);
        LibUtilities::PointsKey pkey2(m_mesh->m_nummode + o, pType.first);

        const int pDim  = pkey1.GetPointsDim();
        const int order = m_mesh->m_nummode - 1;

        Array<OneD, Array<OneD, NekDouble> > u1(pDim), u2(pDim);

        switch (pDim)
        {
            case 2:
            {
                LibUtilities::PointsManager()[pkey1]->GetPoints(u1[0], u1[1]);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2[0], u2[1]);
                break;
            }
            case 3:
            {
                LibUtilities::PointsManager()[pkey1]->GetPoints(u1[0], u1[1],
                                                                u1[2]);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2[0], u2[1],
                                                                u2[2]);
                break;
            }
        }

        der->ptsStd = u1[0].size();
        der->pts    = u2[0].size();

        LibUtilities::NodalUtil *nodalUtil = NULL;

        if (it.first == LibUtilities::eTriangle)
        {
            nodalUtil =
                new LibUtilities::NodalUtilTriangle(order, u1[0], u1[1]);
        }
        else if (it.first == LibUtilities::eQuadrilateral)
        {
            nodalUtil = new LibUtilities::NodalUtilQuad(order, u1[0], u1[1]);
        }
        else if (it.first == LibUtilities::eTetrahedron)
        {
            nodalUtil = new LibUtilities::NodalUtilTetrahedron(order, u1[0],
                                                               u1[1], u1[2]);
        }
        else if (it.first == LibUtilities::ePrism)
        {
            nodalUtil =
                new LibUtilities::NodalUtilPrism(order, u1[0], u1[1], u1[2]);
        }
        else if (it.first == LibUtilities::eHexahedron)
        {
            nodalUtil =
                new LibUtilities::NodalUtilHex(order, u1[0], u1[1], u1[2]);
        }
        else
        {
            ASSERTL0(false,
                     "Unknown element type for derivative utility setup");
        }

        NekMatrix<NekDouble> interp = *nodalUtil->GetInterpolationMatrix(u2);
        NekMatrix<NekDouble> Vandermonde  = *nodalUtil->GetVandermonde();
        NekMatrix<NekDouble> VandermondeI = Vandermonde;
        VandermondeI.Invert();

        for (int i = 0; i < pDim; ++i)
        {
            der->VdmDStd[i] =
                *nodalUtil->GetVandermondeForDeriv(i) * VandermondeI;
            der->VdmD[i] = interp * der->VdmDStd[i];
        }

        Array<OneD, NekDouble> qds =
            LibUtilities::PointsManager()[pkey2]->GetW();
        NekVector<NekDouble> quadWi(qds);
        der->quadW = quadWi;

        ret[it.first] = der;
        delete nodalUtil;
    }

    return ret;
}



struct NodeComparator
{
    const vector<int> & value_vector;

    NodeComparator(const vector<int> & val_vec):
        value_vector(val_vec) {}

    bool operator()(int i1, int i2)
    {
        return value_vector[i1] > value_vector[i2];
    }
};



vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes(
    vector<ElUtilSharedPtr> elLock)
{

    // create set of nodes to be ignored and hence not included in the coloursets
    NodeSet ignoredNodes;
    for (int i = 0; i < elLock.size(); i++)
    {
        vector<NodeSharedPtr> nodes;
        elLock[i]->GetEl()->GetCurvedNodes(nodes);
        for (int j = 0; j < nodes.size(); j++)
        {
            ignoredNodes.insert(nodes[j]);
        }
    }

    // create set of nodes which are at the boundary and hence not included in
    // the colourset
    NodeSet boundaryNodes;

    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            for(auto &edge : m_mesh->m_edgeSet)
            {
                if(edge->m_elLink.size() == 2)
                {
                    continue;
                }

                boundaryNodes.insert(edge->m_n1);
                boundaryNodes.insert(edge->m_n2);
                for(int i = 0; i < edge->m_edgeNodes.size(); i++)
                {
                    boundaryNodes.insert(edge->m_edgeNodes[i]);
                }
            }
            break;
        }
        case 3:
        {
            if(!m_mesh->m_cad)
            {
                for (auto &face : m_mesh->m_faceSet)
                {
                    if (face->m_elLink.size() == 2)
                    {
                        continue;
                    }

                    vector<NodeSharedPtr> vs = face->m_vertexList;
                    for (int j = 0; j < vs.size(); j++)
                    {
                        boundaryNodes.insert(vs[j]);
                    }

                    vector<EdgeSharedPtr> es = face->m_edgeList;
                    for (int j = 0; j < es.size(); j++)
                    {
                        for (int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                        {
                            boundaryNodes.insert(es[j]->m_edgeNodes[k]);
                        }
                    }

                    for (int i = 0; i < face->m_faceNodes.size(); i++)
                    {
                        boundaryNodes.insert(face->m_faceNodes[i]);
                    }
                }
            }
            else
            {
                // If we have CAD therefore the only fixed nodes exist on
                // vertices only
                for (auto &node : m_mesh->m_vertexSet)
                {
                    if(node->GetNumCadCurve() > 1)
                    {
                        boundaryNodes.insert(node);
                    }
                }
            }
            break;
        }
        default:
            ASSERTL0(false,"space dim issue");
    }

    // Create vector of free nodes which "remain", hence will be included in the
    // coloursets
    vector<NodeSharedPtr> remainEdgeVertex;
    vector<NodeSharedPtr> remainFace;
    vector<NodeSharedPtr> remainVolume;
    m_res->nDoF = 0;

    // check if vertex nodes are in boundary or ignored nodes, otherwise add to
    // EDGE-VERTEX remain nodes
    for (auto &node : m_mesh->m_vertexSet)
    {
        if (boundaryNodes.find(node) == boundaryNodes.end() &&
            ignoredNodes.find(node) == ignoredNodes.end())
        {
            remainEdgeVertex.push_back(node);
            if (node->GetNumCadCurve() == 1)
            {
                m_res->nDoF++;
            }
            else if (node->GetNumCADSurf() == 1)
            {
                m_res->nDoF += 2;
            }
            else
            {
                m_res->nDoF += m_mesh->m_spaceDim;
            }
        }
    }

    // check if edge nodes are in boundary or ignored nodes, otherwise add to
    // EDGE-VERTEX remain nodes
    for (auto &edge : m_mesh->m_edgeSet)
    {
        vector<NodeSharedPtr> &n = edge->m_edgeNodes;
        for (int j = 0; j < n.size(); j++)
        {
            if (boundaryNodes.find(n[j]) == boundaryNodes.end() &&
                ignoredNodes.find(n[j]) == ignoredNodes.end())
            {
                remainEdgeVertex.push_back(n[j]);
                if (n[j]->GetNumCadCurve() == 1)
                {
                    m_res->nDoF++;
                }
                else if (n[j]->GetNumCADSurf() == 1)
                {
                    m_res->nDoF += 2;
                }
                else
                {
                    m_res->nDoF += m_mesh->m_spaceDim;
                }
            }
        }
    }

    // check if face nodes are in boundary or ignored nodes, otherwise add to
    // FACE remain nodes
    for (auto &face : m_mesh->m_faceSet)
    {
        vector<NodeSharedPtr> &n = face->m_faceNodes;
        for (int j = 0; j < n.size(); j++)
        {
            if (boundaryNodes.find(n[j]) == boundaryNodes.end() &&
                ignoredNodes.find(n[j]) == ignoredNodes.end())
            {
                remainFace.push_back(n[j]);
                if (n[j]->GetNumCADSurf() == 1)
                {
                    m_res->nDoF += 2;
                }
                else
                {
                    m_res->nDoF += m_mesh->m_spaceDim;
                }
            }
        }
    }

    // check if volume nodes are in boundary or ignored nodes, otherwise add to
    // VOLUME remain nodes
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> ns =
            m_mesh->m_element[m_mesh->m_expDim][i]->GetVolumeNodes();
        for (int j = 0; j < ns.size(); j++)
        {
            if (boundaryNodes.find(ns[j]) == boundaryNodes.end() &&
                ignoredNodes.find(ns[j]) == ignoredNodes.end())
            {
                remainVolume.push_back(ns[j]);
                m_res->nDoF += m_mesh->m_spaceDim;
            }
        }
    }

    // size of all free nodes to be included in the coloursets
    m_res->n = remainEdgeVertex.size()
                + remainFace.size() + remainVolume.size();

    // data structure for coloursets, that will ultimately contain all free
    // nodes
    vector<vector<NodeSharedPtr> > ret;
    vector<vector<NodeSharedPtr> > retPart;

    // edge and vertex nodes
    // create vector num_el of number of associated elements of each node
    vector<int> num_el(remainEdgeVertex.size());
    for (int i = 0; i < remainEdgeVertex.size(); i++)
    {
        //try to find node within all elements
        auto it = m_nodeElMap.find(remainEdgeVertex[i]->m_id);
        vector<ElUtilSharedPtr> &elUtils = it->second;
        num_el[i] = elUtils.size();
    }
    // finding the permutation according to num_el
    vector<int> permNode(remainEdgeVertex.size());
    for (int i = 0; i < remainEdgeVertex.size(); ++i)
    {
        permNode[i] = i;
    }
    std::sort(permNode.begin(), permNode.end(), NodeComparator(num_el));
    // applying the permutation to remainEdgeVertex
    vector<NodeSharedPtr> remainEdgeVertexSort(remainEdgeVertex.size());
    for (int i = 0; i < remainEdgeVertex.size(); ++i)
    {
        int j = permNode[i];
        remainEdgeVertexSort[i] = remainEdgeVertex[j];
    }

    retPart = CreateColoursets(remainEdgeVertexSort);
    if(m_mesh->m_verbose)
    {
        cout << "Number of Edge/Vertex Coloursets: " << retPart.size() << endl;
    }
    for (int i = 0; i < retPart.size(); i++)
    {
        ret.push_back(retPart[i]);
    }

    // face nodes
    retPart = CreateColoursets(remainFace);
    if(m_mesh->m_verbose)
    {
        cout << "Number of Face Coloursets: " << retPart.size() << endl;
    }
    for (int i = 0; i < retPart.size(); i++)
    {
        ret.push_back(retPart[i]);
    }

    // volume nodes
    retPart = CreateColoursets(remainVolume);
    if(m_mesh->m_verbose)
    {
        cout << "Number of Volume Coloursets: " << retPart.size() << endl;
    }
    for (int i = 0; i < retPart.size(); i++)
    {
        ret.push_back(retPart[i]);
    }


    if(m_mesh->m_verbose)
    {
        cout << endl;
    }

    return ret;
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::CreateColoursets(
         vector<NodeSharedPtr> remain)
{
    vector<vector<NodeSharedPtr> > retPart;

    // loop until all free nodes have been sorted
    while (remain.size() > 0)
    {
        vector<NodeSharedPtr> layer; // one colourset
        set<int> locked;
        set<int> completed;
        for (int i = 0; i < remain.size(); i++)
        {
            // Try to find node within all elements
            auto it = m_nodeElMap.find(remain[i]->m_id);
            ASSERTL0(it != m_nodeElMap.end(), "could not find node");

            // identify the vector of all associated elements of the node
            vector<ElUtilSharedPtr> &elUtils = it->second;

            // suppose node is not locked
            bool islocked = false;

            // loop over all associated elements of the node
            for (int j = 0; j < elUtils.size(); j++)
            {
                // check all nodes of the element. if node is within the set of
                // locked nodes then lock node and go to the next node
                if (locked.find(elUtils[j]->GetId()) != locked.end())
                {
                    islocked = true;
                    break;
                }
            }

            // if the node is not locked, insert it into the colourset and
            // insert sorted node into the completed list. Then, loop over all
            // other nodes of the same element and mark them as locked.
            if (!islocked)
            {
                layer.push_back(remain[i]);
                completed.insert(remain[i]->m_id);
                for (int j = 0; j < elUtils.size(); j++)
                {
                    locked.insert(elUtils[j]->GetId());
                }
            }
        }

        // identify nodes which are not sorted, yet and create new "remain"
        // vector
        vector<NodeSharedPtr> tmp = remain;
        remain.clear();
        for (int i = 0; i < tmp.size(); i++)
        {
            if (completed.find(tmp[i]->m_id) == completed.end())
            {
                remain.push_back(tmp[i]);
            }
        }

        // include layer or colourset into vector of coloursets
        retPart.push_back(layer);

        // print out progress
        if(m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                m_res->n - remain.size(), m_res->n, "Node Coloring");
        }

    }
    return retPart;
}


void ProcessVarOpti::GetElementMap(
    int o, map<LibUtilities::ShapeType, DerivUtilSharedPtr> derMap)
{
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);
        ElUtilSharedPtr d = std::shared_ptr<ElUtil>(new ElUtil(
            el, derMap[el->GetShapeType()], m_res, m_mesh->m_nummode, o));
        m_dataSet.push_back(d);
    }

    if (m_config["scalingfile"].beenSet)
    {
        LibUtilities::Interpolator interp =
            GetScalingFieldFromFile(
                m_config["scalingfile"].as<string>().c_str());

        for (int i = 0; i < m_dataSet.size(); ++i)
        {
            m_dataSet[i]->SetScaling(interp);
        }
    }

    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);

        for (int j = 0; j < ns.size(); j++)
        {
            m_nodeElMap[ns[j]->m_id].push_back(m_dataSet[i]);
        }

        ASSERTL0(derMap[el->GetShapeType()]->ptsStd == ns.size(),
                 "mismatch node count");
    }
}

vector<ElUtilSharedPtr> ProcessVarOpti::GetLockedElements(NekDouble thres)
{
    vector<ElUtilSharedPtr> elBelowThres;
    for (int i = 0; i < m_dataSet.size(); ++i)
    {
        if (m_dataSet[i]->GetScaledJac() < thres)
        {
            elBelowThres.push_back(m_dataSet[i]);
        }
    }

    std::unordered_set<int> inmesh;
    vector<ElUtilSharedPtr> totest;

    for (int i = 0; i < elBelowThres.size(); i++)
    {
        auto t = inmesh.insert(elBelowThres[i]->GetId());

        vector<FaceSharedPtr> f = elBelowThres[i]->GetEl()->GetFaceList();
        for (int j = 0; j < f.size(); j++)
        {
            for (int k = 0; k < f[j]->m_elLink.size(); k++)
            {
                if (f[j]->m_elLink[k].first.lock()->GetId() ==
                    elBelowThres[i]->GetId())
                {
                    continue;
                }

                t = inmesh.insert(f[j]->m_elLink[k].first.lock()->GetId());
                if (t.second)
                {
                    totest.push_back(
                        m_dataSet[f[j]->m_elLink[k].first.lock()->GetId()]);
                }
            }
        }
    }

    for (int i = 0; i < 6; i++)
    {
        vector<ElUtilSharedPtr> tmp = totest;
        totest.clear();
        for (int j = 0; j < tmp.size(); j++)
        {
            vector<FaceSharedPtr> f = tmp[j]->GetEl()->GetFaceList();
            for (int k = 0; k < f.size(); k++)
            {
                for (int l = 0; l < f[k]->m_elLink.size(); l++)
                {
                    if (f[k]->m_elLink[l].first.lock()->GetId() ==
                        tmp[j]->GetId())
                    {
                        continue;
                    }

                    auto t = inmesh.insert(
                        f[k]->m_elLink[l].first.lock()->GetId());
                    if (t.second)
                    {
                        totest.push_back(
                            m_dataSet[f[k]->m_elLink[l].first.lock()->GetId()]);
                    }
                }
            }
        }
    }

    // now need to invert the list
    vector<ElUtilSharedPtr> ret;
    for (int i = 0; i < m_dataSet.size(); ++i)
    {
        if (inmesh.find(m_dataSet[i]->GetId()) == inmesh.end())
        {
            ret.push_back(m_dataSet[i]);
        }
    }

    return ret;
}

void ProcessVarOpti::RemoveLinearCurvature()
{
    for(int i = 0; i < m_dataSet.size(); i++)
    {
        if(m_dataSet[i]->GetScaledJac() > 0.999)
        {
            ElementSharedPtr el = m_dataSet[i]->GetEl();
            vector<NodeSharedPtr> ns;
            el->SetVolumeNodes(ns);
        }
    }

    map<int, vector<FaceSharedPtr> > edgeToFace;

    for(auto &face : m_mesh->m_faceSet)
    {
        bool rm = true;
        for(int i = 0; i < face->m_elLink.size(); i++)
        {
            int id = face->m_elLink[i].first.lock()->GetId();
            if(m_dataSet[id]->GetScaledJac() <= 0.999)
            {
                rm = false;
                break;
            }
        }
        if(rm)
        {
            face->m_faceNodes.clear();
        }

        vector<EdgeSharedPtr> es = face->m_edgeList;
        for(int i = 0; i < es.size(); i++)
        {
            edgeToFace[es[i]->m_id].push_back(face);
        }
    }

    for(auto &edge : m_mesh->m_edgeSet)
    {
        auto it = edgeToFace.find(edge->m_id);
        ASSERTL0(it != edgeToFace.end(),"not found");
        bool rm = true;
        for(int i = 0; i < it->second.size(); i++)
        {
            if(it->second[i]->m_faceNodes.size() > 0)
            {
                rm = false;
                break;
            }
        }
        if(rm)
        {
            edge->m_edgeNodes.clear();
        }
    }
}

LibUtilities::Interpolator ProcessVarOpti::GetScalingFieldFromFile(string file)
{
    vector<vector<NekDouble> > data;

    ifstream f;
    f.open(file);
    ASSERTL0(f.is_open(), "No such scaling file")

    string fline;

    while (!f.eof())
    {
        getline(f, fline);

        vector<string> tmp;
        boost::split(tmp, fline, boost::is_any_of(" "));

        int i = 0;
        while (i < tmp.size())
        {
            if (tmp[i].size() == 0)
            {
                tmp.erase(tmp.begin() + i);
            }
            else
            {
                ++i;
            }
        }

        if (tmp.size() < 4)
        {
            continue;
        }

        vector<NekDouble> tmpD;
        tmpD.push_back(boost::lexical_cast<NekDouble>(tmp[0])); // x
        tmpD.push_back(boost::lexical_cast<NekDouble>(tmp[1])); // y
        tmpD.push_back(boost::lexical_cast<NekDouble>(tmp[3])); // scaling

        data.push_back(tmpD);
    }

    int dim = m_mesh->m_expDim;

    Array<OneD, Array<OneD, NekDouble> > inPts(dim + 1);
    for (int i = 0; i < dim + 1; ++i)
    {
        inPts[i] = Array<OneD, NekDouble>(data.size());

        for (int j = 0; j < data.size(); ++j)
        {
            inPts[i][j] = data[j][i];
        }
    }

    return GetField(inPts);
}

LibUtilities::Interpolator ProcessVarOpti::GetField(
    Array<OneD, Array<OneD, NekDouble> > inPts)
{
    int dim = m_mesh->m_expDim;

    vector<string> fieldNames;
    fieldNames.push_back("");

    map<LibUtilities::PtsInfo, int> ptsInfo = LibUtilities::NullPtsInfoMap;

    PtsFieldSharedPtr inField = MemoryManager<LibUtilities::PtsField>
        ::AllocateSharedPtr(dim, fieldNames, inPts, ptsInfo);

    Array<OneD, Array<OneD, NekDouble> > dummyPts(dim + 1);
    for (int i = 0; i < dim + 1; ++i)
    {
        dummyPts[i] = Array<OneD, NekDouble>(0);
    }

    PtsFieldSharedPtr dummyField = MemoryManager<LibUtilities::PtsField>
        ::AllocateSharedPtr(dim, fieldNames, dummyPts, ptsInfo);

    LibUtilities::Interpolator ret;
    ret.Interpolate(inField, dummyField);

    return ret;
}
}
}
