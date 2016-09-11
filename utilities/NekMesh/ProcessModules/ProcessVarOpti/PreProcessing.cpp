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

#include "ProcessVarOpti.h"

#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace Utilities
{

void ProcessVarOpti::BuildDerivUtil()
{
    //build Vandermonde information
    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            {
                LibUtilities::ShapeType st = LibUtilities::eTriangle;
                derivUtil[st] = boost::shared_ptr<DerivUtil>(new DerivUtil);
                derivUtil[st]->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;

                LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                              LibUtilities::eNodalTriElec);
                LibUtilities::PointsKey pkey2(m_mesh->m_nummode + 4,
                                              LibUtilities::eNodalTriSPI);
                Array<OneD, NekDouble> u1, v1, u2, v2;

                LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);

                derivUtil[st]->ptsHigh = u2.num_elements();

                LibUtilities::NodalUtilTriangle nodalTri(
                    m_mesh->m_nummode - 1, u1, v1);

                Array<OneD, Array<OneD, NekDouble> > uv2(2), uv1(2);
                uv1[0] = u1;
                uv1[1] = v1;
                uv2[0] = u2;
                uv2[1] = v2;

                NekMatrix<NekDouble> interp = *nodalTri.GetInterpolationMatrix(uv2);

                NekMatrix<NekDouble> Vandermonde = *nodalTri.GetVandermonde();
                NekMatrix<NekDouble> VandermondeI = Vandermonde;
                VandermondeI.Invert();

                derivUtil[st]->VdmDL[0] = *nodalTri.GetVandermondeForDeriv(0) * VandermondeI;
                derivUtil[st]->VdmDL[1] = *nodalTri.GetVandermondeForDeriv(1) * VandermondeI;
                derivUtil[st]->VdmD[0] = interp * derivUtil[st]->VdmDL[0];
                derivUtil[st]->VdmD[1] = interp * derivUtil[st]->VdmDL[1];
                //derivUtil->quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
                Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
                NekVector<NekDouble> quadWi(qds);
                derivUtil[st]->quadW = quadWi;
            }
            {
                LibUtilities::ShapeType st = LibUtilities::eQuadrilateral;
                derivUtil[st] = boost::shared_ptr<DerivUtil>(new DerivUtil);
                derivUtil[st]->ptsLow  = m_mesh->m_nummode*m_mesh->m_nummode;

                LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                              LibUtilities::eGaussLobattoLegendre);
                LibUtilities::PointsKey pkey2(m_mesh->m_nummode + 4,
                                              LibUtilities::eNodalQuadSPI);
                Array<OneD, NekDouble> u1(derivUtil[st]->ptsLow), v1(derivUtil[st]->ptsLow), u2, v2, tmp;

                LibUtilities::PointsManager()[pkey1]->GetPoints(tmp);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);

                //Get curved nodes ordering for a Quadrilateral is left to right
                //ordering so the we need to use gll and redistrobute it.
                for(int i = 0, ct = 0; i < m_mesh->m_nummode; i++)
                {
                    for(int j = 0; j < m_mesh->m_nummode; j++, ct++)
                    {
                        u1[ct] = tmp[j];
                        v1[ct] = tmp[i];
                    }
                }

                derivUtil[st]->ptsHigh = u2.num_elements();

                LibUtilities::NodalUtilQuad nodalQuad(
                    m_mesh->m_nummode - 1, u1, v1);

                Array<OneD, Array<OneD, NekDouble> > uv2(2), uv1(2);
                uv1[0] = u1;
                uv1[1] = v1;
                uv2[0] = u2;
                uv2[1] = v2;

                NekMatrix<NekDouble> interp = *nodalQuad.GetInterpolationMatrix(uv2);

                NekMatrix<NekDouble> Vandermonde = *nodalQuad.GetVandermonde();
                NekMatrix<NekDouble> VandermondeI = Vandermonde;
                VandermondeI.Invert();

                derivUtil[st]->VdmDL[0] = *nodalQuad.GetVandermondeForDeriv(0) * VandermondeI;
                derivUtil[st]->VdmDL[1] = *nodalQuad.GetVandermondeForDeriv(1) * VandermondeI;
                derivUtil[st]->VdmD[0] = interp * derivUtil[st]->VdmDL[0];
                derivUtil[st]->VdmD[1] = interp * derivUtil[st]->VdmDL[1];
                //derivUtil->quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
                Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
                NekVector<NekDouble> quadWi(qds);
                derivUtil[st]->quadW = quadWi;

                derivUtil[st]->ptx = u2;
                derivUtil[st]->pty = v2;
            }
        }
        break;
        case 3:
        {
            {
                LibUtilities::ShapeType st = LibUtilities::eTetrahedron;
                derivUtil[st] = boost::shared_ptr<DerivUtil>(new DerivUtil);
                derivUtil[st]->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)*(m_mesh->m_nummode+2)/6;
                LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                              LibUtilities::eNodalTetElec);
                LibUtilities::PointsKey pkey2(m_mesh->m_nummode+4,
                                              LibUtilities::eNodalTetSPI);
                Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
                LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

                derivUtil[st]->ptsHigh = u2.num_elements();

                LibUtilities::NodalUtilTetrahedron nodalTet(
                    m_mesh->m_nummode - 1, u1, v1, w1);

                Array<OneD, Array<OneD, NekDouble> > uv2(3), uv1(3);
                uv2[0] = u2;
                uv2[1] = v2;
                uv2[2] = w2;
                uv1[0] = u1;
                uv1[1] = v1;
                uv1[2] = w1;

                NekMatrix<NekDouble> interp = *nodalTet.GetInterpolationMatrix(uv2);
                NekMatrix<NekDouble> Vandermonde = *nodalTet.GetVandermonde();
                NekMatrix<NekDouble> VandermondeI = Vandermonde;
                VandermondeI.Invert();

                derivUtil[st]->VdmDL[0] = *nodalTet.GetVandermondeForDeriv(0) * VandermondeI;
                derivUtil[st]->VdmDL[1] = *nodalTet.GetVandermondeForDeriv(1) * VandermondeI;
                derivUtil[st]->VdmDL[2] = *nodalTet.GetVandermondeForDeriv(2) * VandermondeI;

                derivUtil[st]->VdmD[0] = interp * derivUtil[st]->VdmDL[0];
                derivUtil[st]->VdmD[1] = interp * derivUtil[st]->VdmDL[1];
                derivUtil[st]->VdmD[2] = interp * derivUtil[st]->VdmDL[2];
                Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
                NekVector<NekDouble> quadWi(qds);
                derivUtil[st]->quadW = quadWi;
            }

            {
                LibUtilities::ShapeType st = LibUtilities::ePrism;
                derivUtil[st] = boost::shared_ptr<DerivUtil>(new DerivUtil);
                derivUtil[st]->ptsLow  = m_mesh->m_nummode*m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;
                LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                              LibUtilities::eNodalPrismElec);
                LibUtilities::PointsKey pkey2(m_mesh->m_nummode+4,
                                              LibUtilities::eNodalPrismSPI);
                Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
                LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

                derivUtil[st]->ptsHigh = u2.num_elements();

                LibUtilities::NodalUtilPrism nodalPrism(
                    m_mesh->m_nummode - 1, u1, v1, w1);

                Array<OneD, Array<OneD, NekDouble> > uv2(3), uv1(3);
                uv2[0] = u2;
                uv2[1] = v2;
                uv2[2] = w2;
                uv1[0] = u1;
                uv1[1] = v1;
                uv1[2] = w1;

                NekMatrix<NekDouble> interp = *nodalPrism.GetInterpolationMatrix(uv2);
                NekMatrix<NekDouble> Vandermonde = *nodalPrism.GetVandermonde();
                NekMatrix<NekDouble> VandermondeI = Vandermonde;
                VandermondeI.Invert();

                derivUtil[st]->VdmDL[0] = *nodalPrism.GetVandermondeForDeriv(0) * VandermondeI;
                derivUtil[st]->VdmDL[1] = *nodalPrism.GetVandermondeForDeriv(1) * VandermondeI;
                derivUtil[st]->VdmDL[2] = *nodalPrism.GetVandermondeForDeriv(2) * VandermondeI;

                derivUtil[st]->VdmD[0] = interp * derivUtil[st]->VdmDL[0];
                derivUtil[st]->VdmD[1] = interp * derivUtil[st]->VdmDL[1];
                derivUtil[st]->VdmD[2] = interp * derivUtil[st]->VdmDL[2];
                Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
                NekVector<NekDouble> quadWi(qds);
                derivUtil[st]->quadW = quadWi;

                derivUtil[st]->ptx = u2;
                derivUtil[st]->pty = v2;
                derivUtil[st]->ptz = w2;
            }

            {
                LibUtilities::ShapeType st = LibUtilities::eHexahedron;
                derivUtil[st] = boost::shared_ptr<DerivUtil>(new DerivUtil);
                derivUtil[st]->ptsLow  = m_mesh->m_nummode*m_mesh->m_nummode*m_mesh->m_nummode;
                LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                              LibUtilities::eGaussLobattoLegendre);
                LibUtilities::PointsKey pkey2(m_mesh->m_nummode+4,
                                              LibUtilities::eGaussLobattoLegendre);
                Array<OneD, NekDouble> u1(m_mesh->m_nummode), v1(m_mesh->m_nummode), u2(m_mesh->m_nummode), v2, w1, w2, tmp;
                LibUtilities::PointsManager()[pkey1]->GetPoints(tmp);
                LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

                derivUtil[st]->ptsHigh = u2.num_elements();

                for(int i = 0, ct = 0; i < m_mesh->m_nummode; i++)
                {
                    for(int j = 0; j < m_mesh->m_nummode; j++)
                    {
                        for(int k = 0; k < m_mesh->m_nummode; k++, ct++)
                        {
                            u1[ct] = tmp[k];
                            v1[ct] = tmp[j];
                            w1[ct] = tmp[i];
                        }
                    }
                }

                LibUtilities::NodalUtilHex nodalHex(
                    m_mesh->m_nummode - 1, u1, v1, w1);

                Array<OneD, Array<OneD, NekDouble> > uv2(3), uv1(3);
                uv2[0] = u2;
                uv2[1] = v2;
                uv2[2] = w2;
                uv1[0] = u1;
                uv1[1] = v1;
                uv1[2] = w1;

                NekMatrix<NekDouble> interp = *nodalHex.GetInterpolationMatrix(uv2);
                NekMatrix<NekDouble> Vandermonde = *nodalHex.GetVandermonde();
                NekMatrix<NekDouble> VandermondeI = Vandermonde;
                VandermondeI.Invert();

                derivUtil[st]->VdmDL[0] = *nodalHex.GetVandermondeForDeriv(0) * VandermondeI;
                derivUtil[st]->VdmDL[1] = *nodalHex.GetVandermondeForDeriv(1) * VandermondeI;
                derivUtil[st]->VdmDL[2] = *nodalHex.GetVandermondeForDeriv(2) * VandermondeI;

                derivUtil[st]->VdmD[0] = interp * derivUtil[st]->VdmDL[0];
                derivUtil[st]->VdmD[1] = interp * derivUtil[st]->VdmDL[1];
                derivUtil[st]->VdmD[2] = interp * derivUtil[st]->VdmDL[2];
                Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
                NekVector<NekDouble> quadWi(qds);
                derivUtil[st]->quadW = quadWi;

                derivUtil[st]->ptx = u2;
                derivUtil[st]->pty = v2;
                derivUtil[st]->ptz = w2;
            }
        }
    }
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes(vector<ElementSharedPtr> elLock)
{
    NodeSet ignoredNodes;
    for(int i = 0; i < elLock.size(); i++)
    {
        vector<NodeSharedPtr> nodes;
        elLock[i]->GetCurvedNodes(nodes);
        for(int j = 0; j < nodes.size(); j++)
        {
            ignoredNodes.insert(nodes[j]);
        }
    }

    //this figures out the dirclet nodes and colors the others into paralell sets
    NodeSet boundaryNodes;

    if(!m_mesh->m_hasCAD)
    {
        switch (m_mesh->m_spaceDim)
        {
            case 2:
            {
                EdgeSet::iterator it;
                for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
                {
                    if((*it)->m_elLink.size() == 2)
                    {
                        continue;
                    }

                    boundaryNodes.insert((*it)->m_n1);
                    boundaryNodes.insert((*it)->m_n2);
                    for(int i = 0; i < (*it)->m_edgeNodes.size(); i++)
                    {
                        boundaryNodes.insert((*it)->m_edgeNodes[i]);
                    }
                }
                break;
            }
            case 3:
            {
                FaceSet::iterator it;
                for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
                {
                    if((*it)->m_elLink.size() == 2)
                    {
                        continue;
                    }

                    vector<NodeSharedPtr> vs = (*it)->m_vertexList;
                    for(int j = 0; j < vs.size(); j++)
                    {
                        boundaryNodes.insert(vs[j]);
                    }

                    vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                    for(int j = 0; j < es.size(); j++)
                    {
                        for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                        {
                            boundaryNodes.insert(es[j]->m_edgeNodes[k]);
                        }
                    }

                    for(int i = 0; i < (*it)->m_faceNodes.size(); i++)
                    {
                        boundaryNodes.insert((*it)->m_faceNodes[i]);
                    }
                }
                break;
            }
            default:
                ASSERTL0(false,"space dim issue");
        }
    }
    else
    {
        //if we have CAD we are 3D and therefore the only fixed nodes exist on vertices only
        NodeSet::iterator nit;
        for (nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); ++nit)
        {
            if((*nit)->GetNumCadCurve() > 1)
            {
                boundaryNodes.insert((*nit));
            }
        }
    }

    vector<NodeSharedPtr> remain;
    res->nDoF = 0;

    NodeSet::iterator nit;
    for (nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); ++nit)
    {
        NodeSet::iterator nit2 = boundaryNodes.find(*nit);
        NodeSet::iterator nit3 = ignoredNodes.find(*nit);
        if(nit2 == boundaryNodes.end() && nit3 == ignoredNodes.end())
        {
            remain.push_back(*nit);
            if((*nit)->GetNumCadCurve() == 1)
            {
                res->nDoF++;
            }
            else if((*nit)->GetNumCADSurf() == 1)
            {
                res->nDoF += 2;
            }
            else
            {
                res->nDoF += 3;
            }
        }
    }

    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        vector<NodeSharedPtr> n = (*eit)->m_edgeNodes;
        for(int j = 0; j < n.size(); j++)
        {
            NodeSet::iterator nit2 = boundaryNodes.find(n[j]);
            NodeSet::iterator nit3 = ignoredNodes.find(n[j]);
            if(nit2 == boundaryNodes.end() && nit3 == ignoredNodes.end())
            {
                remain.push_back(n[j]);
                if(n[j]->GetNumCadCurve() == 1)
                {
                    res->nDoF++;
                }
                else if(n[j]->GetNumCADSurf() == 1)
                {
                    res->nDoF += 2;
                }
                else
                {
                    res->nDoF += 3;
                }
            }
        }
    }

    FaceSet::iterator fit;
    for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
    {
        for(int j = 0; j < (*fit)->m_faceNodes.size(); j++)
        {
            NodeSet::iterator nit2 = boundaryNodes.find((*fit)->m_faceNodes[j]);
            NodeSet::iterator nit3 = ignoredNodes.find((*fit)->m_faceNodes[j]);
            if(nit2 == boundaryNodes.end() && nit3 == ignoredNodes.end())
            {
                remain.push_back((*fit)->m_faceNodes[j]);
                if((*fit)->m_faceNodes[j]->GetNumCADSurf() == 1)
                {
                    res->nDoF += 2;
                }
                else
                {
                    res->nDoF += 3;
                }
            }
        }
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> ns =
            m_mesh->m_element[m_mesh->m_expDim][i]->GetVolumeNodes();
        for(int j = 0; j < ns.size(); j++)
        {
            NodeSet::iterator nit2 = boundaryNodes.find(ns[j]);
            NodeSet::iterator nit3 = ignoredNodes.find(ns[j]);
            if(nit2 == boundaryNodes.end() && nit3 == ignoredNodes.end())
            {
                remain.push_back(ns[j]);
                res->nDoF += 3;
            }
        }
    }

    res->n = remain.size();

    vector<vector<NodeSharedPtr> > ret;

    while (remain.size() > 0)
    {
        vector<NodeSharedPtr> layer;
        set<int> locked;
        set<int> completed;
        for(int i = 0; i < remain.size(); i++)
        {
            NodeElMap::iterator it = nodeElMap.find(remain[i]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");

            vector<ElUtilSharedPtr> &elUtils = it->second;

            bool islocked = false;
            for(int j = 0; j < elUtils.size(); j++)
            {
                set<int>::iterator sit = locked.find(elUtils[j]->GetId());
                if(sit != locked.end())
                {
                    islocked = true;
                    break;
                }
            }
            if(!islocked)
            {
                layer.push_back(remain[i]);
                completed.insert(remain[i]->m_id);
                for(int j = 0; j < elUtils.size(); j++)
                {
                    locked.insert(elUtils[j]->GetId());
                }
            }
        }

        vector<NodeSharedPtr> tmp = remain;
        remain.clear();
        for(int i = 0; i < tmp.size(); i++)
        {
            set<int>::iterator sit = completed.find(tmp[i]->m_id);
            if(sit == completed.end())
            {
                remain.push_back(tmp[i]);
            }
        }
        ret.push_back(layer);
    }
    return ret;
}

void ProcessVarOpti::GetElementMap()
{
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);
        ElUtilSharedPtr d = boost::shared_ptr<ElUtil>(new ElUtil(el,
                                    derivUtil[el->GetShapeType()],
                                    res, m_mesh->m_nummode));
        dataSet.push_back(d);
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);

        for(int j = 0; j < ns.size(); j++)
        {
            nodeElMap[ns[j]->m_id].push_back(dataSet[i]);
        }

        ASSERTL0(derivUtil[dataSet[i]->GetEl()->GetShapeType()]->ptsLow == ns.size(), "mismatch node count");
    }
}

vector<ElementSharedPtr> ProcessVarOpti::GetLockedElements(NekDouble thres)
{
    vector<ElementSharedPtr> elBelowThres;
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> nodes;
        el->GetCurvedNodes(nodes);
        NekDouble mx = -1.0 * numeric_limits<double>::max();
        NekDouble mn =  numeric_limits<double>::max();

        if(m_mesh->m_expDim == 2)
        {
            NekVector<NekDouble> X(nodes.size()),Y(nodes.size());
            for(int j = 0; j < nodes.size(); j++)
            {
                X(j) = nodes[j]->m_x;
                Y(j) = nodes[j]->m_y;
            }

            NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),
                                 x2i(nodes.size()),y2i(nodes.size());

            x1i = derivUtil[el->GetShapeType()]->VdmDL[0]*X;
            y1i = derivUtil[el->GetShapeType()]->VdmDL[0]*Y;
            x2i = derivUtil[el->GetShapeType()]->VdmDL[1]*X;
            y2i = derivUtil[el->GetShapeType()]->VdmDL[1]*Y;

            for(int j = 0; j < nodes.size(); j++)
            {
                NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
                mx = max(mx,jacDet);
                mn = min(mn,jacDet);
            }
        }
        else if(m_mesh->m_expDim == 3)
        {
            NekVector<NekDouble> X(nodes.size()),Y(nodes.size()),Z(nodes.size());
            for(int j = 0; j < nodes.size(); j++)
            {
                X(j) = nodes[j]->m_x;
                Y(j) = nodes[j]->m_y;
                Z(j) = nodes[j]->m_z;
            }

            NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),z1i(nodes.size()),
                                 x2i(nodes.size()),y2i(nodes.size()),z2i(nodes.size()),
                                 x3i(nodes.size()),y3i(nodes.size()),z3i(nodes.size());

            x1i = derivUtil[el->GetShapeType()]->VdmDL[0]*X;
            y1i = derivUtil[el->GetShapeType()]->VdmDL[0]*Y;
            z1i = derivUtil[el->GetShapeType()]->VdmDL[0]*Z;
            x2i = derivUtil[el->GetShapeType()]->VdmDL[1]*X;
            y2i = derivUtil[el->GetShapeType()]->VdmDL[1]*Y;
            z2i = derivUtil[el->GetShapeType()]->VdmDL[1]*Z;
            x3i = derivUtil[el->GetShapeType()]->VdmDL[2]*X;
            y3i = derivUtil[el->GetShapeType()]->VdmDL[2]*Y;
            z3i = derivUtil[el->GetShapeType()]->VdmDL[2]*Z;

            for(int j = 0; j < nodes.size(); j++)
            {
                DNekMat dxdz(3,3,1.0,eFULL);
                dxdz(0,0) = x1i(j);
                dxdz(0,1) = x2i(j);
                dxdz(0,2) = x3i(j);
                dxdz(1,0) = y1i(j);
                dxdz(1,1) = y2i(j);
                dxdz(1,2) = y3i(j);
                dxdz(2,0) = z1i(j);
                dxdz(2,1) = z2i(j);
                dxdz(2,2) = z3i(j);

                NekDouble jacDet = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                                  -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                                  +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

                mx = max(mx,jacDet);
                mn = min(mn,jacDet);
            }
        }

        if(mn/mx < thres)
        {
            elBelowThres.push_back(el);
        }
    }

    boost::unordered_set<int> inmesh;
    pair<boost::unordered_set<int>::iterator, bool> t;
    vector<ElementSharedPtr> totest;

    for (int i = 0; i < elBelowThres.size(); i++)
    {
        t = inmesh.insert(elBelowThres[i]->GetId());

        vector<FaceSharedPtr> f = elBelowThres[i]->GetFaceList();
        for (int j = 0; j < f.size(); j++)
        {
            for (int k = 0; k < f[j]->m_elLink.size(); k++)
            {
                if (f[j]->m_elLink[k].first->GetId() == elBelowThres[i]->GetId())
                    continue;

                t = inmesh.insert(f[j]->m_elLink[k].first->GetId());
                if (t.second)
                {
                    totest.push_back(f[j]->m_elLink[k].first);
                }
            }
        }
    }

    for (int i = 0; i < 12; i++)
    {
        vector<ElementSharedPtr> tmp = totest;
        totest.clear();
        for (int j = 0; j < tmp.size(); j++)
        {
            vector<FaceSharedPtr> f = tmp[j]->GetFaceList();
            for (int k = 0; k < f.size(); k++)
            {
                for (int l = 0; l < f[k]->m_elLink.size(); l++)
                {
                    if (f[k]->m_elLink[l].first->GetId() == tmp[j]->GetId())
                        continue;

                    t = inmesh.insert(f[k]->m_elLink[l].first->GetId());
                    if (t.second)
                    {
                        totest.push_back(f[k]->m_elLink[l].first);
                    }
                }
            }
        }
    }

    //now need to invert the list
    vector<ElementSharedPtr> ret;
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        boost::unordered_set<int>::iterator s = inmesh.find(el->GetId());
        if(s == inmesh.end())
        {
            ret.push_back(el);
        }
    }

    return ret;
}

}
}
