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

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace Utilities
{

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
                set<int>::iterator sit = locked.find(it->second[j]->el->GetId());
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
                    locked.insert(it->second[j]->el->GetId());
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
        ElUtilSharedPtr d = boost::shared_ptr<ElUtil>(new ElUtil(el));
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
        if((*eit)->m_edgeNodes.size() > 0)
        {
            //already high-order just need to Id
            for(int i = 0; i < (*eit)->m_edgeNodes.size(); i++)
            {
                (*eit)->m_edgeNodes[i]->m_id = id++;
            }
            continue;
        }
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
            if((*it)->m_faceNodes.size() > 0)
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

    EvaluateMesh();
}

void ProcessVarOpti::EvaluateMesh()
{
    res->startInv =0;
    res->worstJac = numeric_limits<double>::max();

    if(m_mesh->m_expDim == 2)
    {
        LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                      LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u1, v1;

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);

        NekVector<NekDouble> U1(u1), V1(v1);

        NekMatrix<NekDouble> Vandermonde = LibUtilities::GetVandermonde(U1,V1);
        NekMatrix<NekDouble> VandermondeI = Vandermonde;
        VandermondeI.Invert();
        NekMatrix<NekDouble> VdmDxt =  (
          LibUtilities::GetVandermondeForXDerivative(U1,V1) * VandermondeI);
        NekMatrix<NekDouble> VdmDyt =  (
          LibUtilities::GetVandermondeForYDerivative(U1,V1) * VandermondeI);



        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            vector<NodeSharedPtr> ns;
            el->GetCurvedNodes(ns);

            NekDouble mx = -1.0 * numeric_limits<double>::max();
            NekDouble mn =  numeric_limits<double>::max();


            NekVector<NekDouble> X(u1.num_elements()),Y(u1.num_elements());
            for(int j = 0; j < u1.num_elements(); j++)
            {
                X(j) = ns[j]->m_x;
                Y(j) = ns[j]->m_y;
            }

            NekVector<NekDouble> x1i(u1.num_elements()),y1i(u1.num_elements()),
                                 x2i(u1.num_elements()),y2i(u1.num_elements());

            x1i = VdmDxt*X;
            y1i = VdmDxt*Y;
            x2i = VdmDyt*X;
            y2i = VdmDyt*Y;

            for(int j = 0; j < u1.num_elements(); j++)
            {
                NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
                mx = max(mx,jacDet);
                mn = min(mn,jacDet);
            }

            if(mn < 0)
            {
                res->startInv++;
            }
            res->worstJac = min(res->worstJac,mn/mx);
        }
    }
    else
    {
        LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                      LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u1, v1,w1;

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1,w1);

        NekVector<NekDouble> U1(u1), V1(v1), W1(w1);

        NekMatrix<NekDouble> Vandermonde = LibUtilities::GetTetVandermonde(U1,V1,W1);
        NekMatrix<NekDouble> VandermondeI = Vandermonde;
        VandermondeI.Invert();
        NekMatrix<NekDouble> VdmDxt =  (
          LibUtilities::GetVandermondeForTetXDerivative(U1,V1,W1) * VandermondeI);
        NekMatrix<NekDouble> VdmDyt =  (
          LibUtilities::GetVandermondeForTetYDerivative(U1,V1,W1) * VandermondeI);
        NekMatrix<NekDouble> VdmDzt =  (
          LibUtilities::GetVandermondeForTetZDerivative(U1,V1,W1) * VandermondeI);

        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            vector<NodeSharedPtr> ns;
            el->GetCurvedNodes(ns);

            NekDouble mx = -1.0 * numeric_limits<double>::max();
            NekDouble mn =  numeric_limits<double>::max();

            NekVector<NekDouble> X(u1.num_elements()),Y(u1.num_elements()),Z(u1.num_elements());
            for(int j = 0; j < u1.num_elements(); j++)
            {
                X(j) = ns[j]->m_x;
                Y(j) = ns[j]->m_y;
                Z(j) = ns[j]->m_z;
            }

            NekVector<NekDouble> x1i(u1.num_elements()),y1i(u1.num_elements()),z1i(u1.num_elements()),
                                 x2i(u1.num_elements()),y2i(u1.num_elements()),z2i(u1.num_elements()),
                                 x3i(u1.num_elements()),y3i(u1.num_elements()),z3i(u1.num_elements());

            x1i = VdmDxt*X;
            y1i = VdmDxt*Y;
            z1i = VdmDxt*Z;
            x2i = VdmDyt*X;
            y2i = VdmDyt*Y;
            z2i = VdmDyt*Z;
            x3i = VdmDzt*X;
            y3i = VdmDzt*Y;
            z3i = VdmDzt*Z;

            for(int j = 0; j < u1.num_elements(); j++)
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



            if(mn < 0)
            {
                res->startInv++;
            }

            res->worstJac = min(res->worstJac,mn/mx);
        }
    }
}

void ProcessVarOpti::WriteStats(string file)
{
    ASSERTL0(file != "", "no file name given");

    ofstream out;
    out.open(file.c_str());
    out << scientific;

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> ns;
        m_mesh->m_element[m_mesh->m_expDim][i]->GetCurvedNodes(ns);
        int pts = ns.size();
        int dim = m_mesh->m_element[m_mesh->m_expDim][i]->GetDim();

        NekDouble mx = -1.0 * numeric_limits<double>::max();
        NekDouble mn = numeric_limits<double>::max();

        if(dim == 2)
        {
            NekVector<NekDouble> X(pts),Y(pts),Z(pts),
                                 x1(pts),y1(pts),
                                 x2(pts),y2(pts);
            for(int i = 0; i < pts; i++)
            {
                X(i) = ns[i]->m_x;
                Y(i) = ns[i]->m_y;
            }

            x1 = derivUtil->VdmD[0]*X;
            y1 = derivUtil->VdmD[0]*Y;
            x2 = derivUtil->VdmD[0]*X;
            y2 = derivUtil->VdmD[0]*Y;

            for(int i = 0; i < pts; i++)
            {
                mx = max(mx, x1(i)*y2(i)-y1(i)*x2(i));
                mn = min(mn, x1(i)*y2(i)-y1(i)*x2(i));
            }

        }
        else
        {
            ASSERTL0(false,"not coded");
        }
        out << mn / mx << endl;
    }
}

}
}
