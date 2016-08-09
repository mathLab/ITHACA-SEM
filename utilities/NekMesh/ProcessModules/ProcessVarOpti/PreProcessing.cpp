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
            derivUtil->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;

            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTriElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTriSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2;

            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);

            derivUtil->ptsHigh =
                LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            LibUtilities::NodalUtilTriangle nodalTri(
                m_mesh->m_nummode - 1, u1, v1);

            Array<OneD, Array<OneD, NekDouble> > uv2(2);
            uv2[0] = u2;
            uv2[1] = v2;

            NekMatrix<NekDouble> interp = *nodalTri.GetInterpolationMatrix(uv2);

            NekMatrix<NekDouble> Vandermonde = *nodalTri.GetVandermonde();
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();

            derivUtil->VdmDL[0] = *nodalTri.GetVandermondeForDeriv(0) * VandermondeI;
            derivUtil->VdmDL[1] = *nodalTri.GetVandermondeForDeriv(1) * VandermondeI;
            derivUtil->VdmD[0] = interp * derivUtil->VdmDL[0];
            derivUtil->VdmD[1] = interp * derivUtil->VdmDL[1];
            //derivUtil->quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            derivUtil->quadW = quadWi;
        }
        break;
        case 3:
        {
            derivUtil->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)*(m_mesh->m_nummode+2)/6;
            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTetElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTetSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

            derivUtil->ptsHigh =
                LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

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

            derivUtil->VdmDL[0] = *nodalTet.GetVandermondeForDeriv(0) * VandermondeI;
            derivUtil->VdmDL[1] = *nodalTet.GetVandermondeForDeriv(1) * VandermondeI;
            derivUtil->VdmDL[2] = *nodalTet.GetVandermondeForDeriv(2) * VandermondeI;

            derivUtil->VdmD[0] = interp * derivUtil->VdmDL[0];
            derivUtil->VdmD[1] = interp * derivUtil->VdmDL[1];
            derivUtil->VdmD[2] = interp * derivUtil->VdmDL[2];
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            derivUtil->quadW = quadWi;

            // Set up derivatives
            derivUtil->basisDeriv = Array<OneD, Array<OneD, NekDouble> >(
                derivUtil->ptsHigh);

            NekVector<NekDouble> tmp(derivUtil->ptsLow);
            NekVector<NekDouble> derivout[3];

            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < derivUtil->ptsLow; ++j)
                {
                    tmp(j) = uv1[i][j];
                }

                derivout[i] = derivUtil->VdmD[0] * tmp;
            }

            for (int i = 0; i < derivUtil->ptsHigh; ++i)
            {
                derivUtil->basisDeriv[i] = Array<OneD, NekDouble>(3);
                for (int j = 0; j < 3; ++j)
                {
                    derivUtil->basisDeriv[i][j] = derivout[j](i);
                }
            }
        }
    }
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes()
{
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
        if(nit2 == boundaryNodes.end())
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
            NodeSet::iterator nit = boundaryNodes.find(n[j]);
            if(nit == boundaryNodes.end())
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
            NodeSet::iterator nit = boundaryNodes.find((*fit)->m_faceNodes[j]);
            if(nit == boundaryNodes.end())
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
            NodeSet::iterator nit = boundaryNodes.find(ns[j]);
            if(nit == boundaryNodes.end())
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
            bool islocked = false;
            for(int j = 0; j < it->second.size(); j++)
            {
                set<int>::iterator sit = locked.find(it->second[j]->GetId());
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
                for(int j = 0; j < it->second.size(); j++)
                {
                    locked.insert(it->second[j]->GetId());
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
        ElUtilSharedPtr d = boost::shared_ptr<ElUtil>(new ElUtil(el, derivUtil,
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
    }
}

void ProcessVarOpti::FillQuadPoints()
{
    int nq = m_mesh->m_nummode;
    int id = m_mesh->m_vertexSet.size();

    LibUtilities::PointsKey ekey(m_mesh->m_nummode,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    EdgeSet::iterator eit;
    FaceSet::iterator fit;

    boost::unordered_map<int, SpatialDomains::Geometry1DSharedPtr> edgeGeoms;
    boost::unordered_map<int, SpatialDomains::Geometry2DSharedPtr> faceGeoms;
    boost::unordered_map<int, SpatialDomains::GeometrySharedPtr> volGeoms;

    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        SpatialDomains::Geometry1DSharedPtr geom =
            (*eit)->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        edgeGeoms[(*eit)->m_id] = geom;
    }

    for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
    {
        SpatialDomains::Geometry2DSharedPtr geom =
            (*fit)->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        faceGeoms[(*fit)->m_id] = geom;
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        SpatialDomains::GeometrySharedPtr geom =
            el->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        volGeoms[el->GetId()] = geom;
    }

    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        // if((*eit)->m_edgeNodes.size() > 0 && (*eit)->m_curveType == LibUtilities::eGaussLobattoLegendre)
        // {
        //     //already high-order just need to Id and double check type
        //     for(int i = 0; i < (*eit)->m_edgeNodes.size(); i++)
        //     {
        //         (*eit)->m_edgeNodes[i]->m_id = id++;
        //     }
        //     continue;
        // }
        SpatialDomains::Geometry1DSharedPtr geom = edgeGeoms[(*eit)->m_id];
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        vector<NodeSharedPtr> hons;

        switch (m_mesh->m_spaceDim)
        {
            case 2:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),0.0)));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
            case 3:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());
                Array<OneD, NekDouble> zc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);
                xmap->BwdTrans(coeffs2,zc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),
                                 xmap->PhysEvaluate(xp,zc))));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
        }
    }

    if(m_mesh->m_expDim == 2)
    {
        //for faces need to do volume nodes
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom = volGeoms[el->GetId()];
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);

            vector<NodeSharedPtr> hons;

            for(int j = 3 + 3*(nq-2); j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),0.0)));
            }

            el->SetVolumeNodes(hons);
            el->SetCurveType(LibUtilities::eNodalTriElec);
        }
    }
    else
    {
        FaceSet::iterator it;
        for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
        {
            if((*it)->m_faceNodes.size() > 0 && (*it)->m_curveType == LibUtilities::eNodalTriElec)
            {
                for(int i = 0; i < (*it)->m_faceNodes.size(); i++)
                {
                    (*it)->m_faceNodes[i]->m_id = id++;
                }
                continue;
            }
            SpatialDomains::Geometry2DSharedPtr geom = faceGeoms[(*it)->m_id];
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            vector<NodeSharedPtr> hons;

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            for(int j = 3 + 3*(nq-2); j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }
            (*it)->m_faceNodes.clear();
            (*it)->m_faceNodes = hons;
            (*it)->m_curveType = LibUtilities::eNodalTriElec;
        }
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom = volGeoms[el->GetId()];
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTetElec);
            Array<OneD, NekDouble> u, v, w;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            vector<NodeSharedPtr> hons;

            //need to finish for tet
            for(int j = 4 + 6*(nq-2) + 4 * ((nq-2)*(nq-3) / 2);
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(3);
                xp[0] = u[j];
                xp[1] = v[j];
                xp[2] = w[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }

            el->SetVolumeNodes(hons);
            el->SetCurveType(LibUtilities::eNodalTetElec);
        }
    }
}

}
}
