////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include "Octree.h"
#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/Module/Module.h>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <boost/algorithm/string.hpp>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void Octree::Process()
{
    Array<OneD, NekDouble> boundingBox = m_mesh->m_cad->GetBoundingBox();

    // build curvature samples
    CompileSourcePointList();

    if (m_mesh->m_verbose)
    {
        cout << "\tCurvature samples: " << m_SPList.size() << endl;
    }

    // make master octant based on the bounding box of the domain
    m_dim = max((boundingBox[1] - boundingBox[0]) / 2.0,
                (boundingBox[3] - boundingBox[2]) / 2.0);

    m_dim = max(m_dim, (boundingBox[5] - boundingBox[4]) / 2.0);

    m_centroid    = Array<OneD, NekDouble>(3);
    m_centroid[0] = (boundingBox[1] + boundingBox[0]) / 2.0;
    m_centroid[1] = (boundingBox[3] + boundingBox[2]) / 2.0;
    m_centroid[2] = (boundingBox[5] + boundingBox[4]) / 2.0;

    m_masteroct = MemoryManager<Octant>::AllocateSharedPtr(
        0, m_centroid[0], m_centroid[1], m_centroid[2], m_dim, m_SPList);

    // begin recersive subdivision
    SubDivide();

    m_octants.clear();
    m_masteroct->CompileLeaves(m_octants);

    if (m_mesh->m_verbose)
    {
        cout << "\tOctants: " << m_octants.size() << endl;
    }

    SmoothSurfaceOctants();

    PropagateDomain();

    SmoothAllOctants();

    if (m_mesh->m_verbose)
    {
        int elem = CountElemt();
        cout << "\tPredicted mesh: " << elem << endl;
    }
}

NekDouble Octree::Query(Array<OneD, NekDouble> loc)
{
    // starting at master octant 0 move through succsesive m_octants which
    // contain the point loc until a leaf is found
    // first search through sourcepoints

    NekDouble tmp = numeric_limits<double>::max();

    for (int i = 0; i < m_lsources.size(); i++)
    {
        if (m_lsources[i].withinRange(loc))
        {
            tmp = min(m_lsources[i].delta, tmp);
        }
    }

    OctantSharedPtr n = m_masteroct;
    int quad = 0;

    bool found = false;

    while (!found)
    {
        Array<OneD, NekDouble> octloc = n->GetLoc();

        if (!(loc[0] < octloc[0]) && // forward
            !(loc[1] < octloc[1]) && // forward
            !(loc[2] < octloc[2]))   // forward
        {
            quad = 0;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] > octloc[2]))   // back
        {
            quad = 1;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] < octloc[2]))   // forward
        {
            quad = 2;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] > octloc[2]))   // back
        {
            quad = 3;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] < octloc[2]))   // forward
        {
            quad = 4;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] > octloc[2]))   // back
        {
            quad = 5;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] < octloc[2]))   // forward
        {
            quad = 6;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] > octloc[2]))   // back
        {
            quad = 7;
        }
        else
        {
            ASSERTL0(false, "Cannot locate quadrant");
        }

        n = n->GetChild(quad);

        if (n->IsLeaf())
        {
            found = true;
        }
    }

    return min(n->GetDelta(), tmp);
}

NekDouble Octree::GetMinDelta()
{
    NekDouble tmp = numeric_limits<double>::max();

    for (int i = 0; i < m_lsources.size(); i++)
    {
        tmp = min(m_lsources[i].delta, tmp);
    }
    return min(m_minDelta, tmp);
}

void Octree::WriteOctree(string nm)
{
    MeshSharedPtr oct = std::shared_ptr<Mesh>(new Mesh());
    oct->m_expDim     = 3;
    oct->m_spaceDim   = 3;
    oct->m_nummode    = 2;

    for (int i = 0; i < m_octants.size(); i++)
    {
        /*if(m_octants[i]->GetLocation() != eOnBoundary)
        {
            continue;
        }*/

        vector<NodeSharedPtr> ns(8);

        ns[0] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eBack),
                                               m_octants[i]->FX(eDown),
                                               m_octants[i]->FX(eRight)));

        ns[1] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eForward),
                                               m_octants[i]->FX(eDown),
                                               m_octants[i]->FX(eRight)));

        ns[2] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eForward),
                                               m_octants[i]->FX(eUp),
                                               m_octants[i]->FX(eRight)));

        ns[3] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eBack),
                                               m_octants[i]->FX(eUp),
                                               m_octants[i]->FX(eRight)));

        ns[4] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eBack),
                                               m_octants[i]->FX(eDown),
                                               m_octants[i]->FX(eLeft)));

        ns[5] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eForward),
                                               m_octants[i]->FX(eDown),
                                               m_octants[i]->FX(eLeft)));

        ns[6] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eForward),
                                               m_octants[i]->FX(eUp),
                                               m_octants[i]->FX(eLeft)));

        ns[7] = std::shared_ptr<Node>(new Node(0, m_octants[i]->FX(eBack),
                                               m_octants[i]->FX(eUp),
                                               m_octants[i]->FX(eLeft)));

        vector<int> tags;
        tags.push_back(0);
        ElmtConfig conf(LibUtilities::eHexahedron, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eHexahedron, conf, ns, tags);
        oct->m_element[3].push_back(E);
    }

    ModuleSharedPtr mod =
        GetModuleFactory().CreateInstance(ModuleKey(eOutputModule, "xml"), oct);
    mod->RegisterConfig("outfile", nm);
    mod->ProcessVertices();
    mod->ProcessEdges();
    mod->ProcessFaces();
    mod->ProcessElements();
    mod->ProcessComposites();
    mod->Process();
}

void Octree::SubDivide()
{
    bool repeat;
    int ct = 1;

    m_numoct = 1;
    m_masteroct->Subdivide(m_masteroct, m_numoct);

    if (m_mesh->m_verbose)
    {
        cout << "\tSubdivide iteration: ";
    }

    do
    {
        if (m_mesh->m_verbose)
        {
            cout << "\r                                                       ";
            cout << "\r";
            cout << "\tSubdivide iteration: " << ct;
            cout.flush();
        }
        ct++;
        repeat = false;
        m_octants.clear();
        // grab a list of the leaves curently in the octree
        m_masteroct->CompileLeaves(m_octants);

        VerifyNeigbours();

        // neeed to create a divide list, in the first list will be m_octants
        // which need to
        // subdivide based on curvature,
        // in the next list will be ocants which need to subdivide to make sure
        // level criteria is statisified for the previous list and so on
        // the list will keep building till no more lists are required.
        // the list will then be iterated through backwards to subdivide all the
        // sublists in turn.
        vector<vector<OctantSharedPtr> > dividelist;
        set<int> inlist;
        // build initial list
        {
            vector<OctantSharedPtr> sublist;
            for (int i = 0; i < m_octants.size(); i++)
            {
                if (m_octants[i]->NeedDivide() &&
                    m_octants[i]->DX() / 4.0 > m_minDelta)
                {
                    sublist.push_back(m_octants[i]);
                    inlist.insert(m_octants[i]->GetId());
                    repeat = true; // if there is a subdivision this whole
                                   // process needs to be repeated
                }
            }
            dividelist.push_back(sublist);
        }
        // then loop over building sublists until no more are required
        int ct2 = 0;
        while (true)
        {
            ct2++;
            vector<OctantSharedPtr> newsublist,
                previouslist = dividelist.back();
            for (int i = 0; i < previouslist.size(); i++)
            {
                map<OctantFace, vector<OctantSharedPtr> > nlist =
                    previouslist[i]->GetNeigbours();
                map<OctantFace, vector<OctantSharedPtr> >::iterator it;
                for (it = nlist.begin(); it != nlist.end(); it++)
                {
                    for (int j = 0; j < it->second.size(); j++)
                    {
                        if (previouslist[i]->DX() < it->second[j]->DX())
                        {
                            set<int>::iterator s =
                                inlist.find(it->second[j]->GetId());
                            if (s == inlist.end())
                            {
                                inlist.insert(it->second[j]->GetId());
                                newsublist.push_back(it->second[j]);
                            }
                        }
                    }
                }
            }
            if (newsublist.size() == 0)
            {
                break;
            }
            else
            {
                dividelist.push_back(newsublist);
            }
        }

        vector<vector<OctantSharedPtr> >::reverse_iterator rit;
        for (rit = dividelist.rbegin(); rit != dividelist.rend(); rit++)
        {
            vector<OctantSharedPtr> currentlist = *rit;
            for (int i = 0; i < currentlist.size(); i++)
            {
                currentlist[i]->Subdivide(currentlist[i], m_numoct);
            }
        }
    } while (repeat);

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }
}

bool Octree::VerifyNeigbours()
{
    // check all octant links to their neighbours
    // at all times in the subdivision the set of neigbours must
    // conform to a set of criteria such as smoothness in size
    // this checks that
    bool error = false;
    for (int i = 0; i < m_octants.size(); i++)
    {
        bool valid = true;
        map<OctantFace, vector<OctantSharedPtr> > nlist =
            m_octants[i]->GetNeigbours();
        map<OctantFace, vector<OctantSharedPtr> >::iterator it;
        for (it = nlist.begin(); it != nlist.end(); it++)
        {
            if (it->second.size() == 0)
            {
                NekDouble expectedfx = 0.0;
                switch (it->first)
                {
                    case eUp:
                        expectedfx = m_centroid[1] + m_dim;
                        break;
                    case eDown:
                        expectedfx = m_centroid[1] - m_dim;
                        break;
                    case eLeft:
                        expectedfx = m_centroid[2] + m_dim;
                        break;
                    case eRight:
                        expectedfx = m_centroid[2] - m_dim;
                        break;
                    case eForward:
                        expectedfx = m_centroid[0] + m_dim;
                        break;
                    case eBack:
                        expectedfx = m_centroid[0] - m_dim;
                        break;
                }
                if (fabs(m_octants[i]->FX(it->first) - expectedfx) > 1E-6)
                {
                    valid = false;
                    cout << "wall neigbour error" << endl;
                    cout << expectedfx << " " << m_octants[i]->FX(it->first)
                         << " " << it->first << endl;
                }
            }
            else if (it->second.size() == 1)
            {
                if (!(m_octants[i]->DX() == it->second[0]->DX() ||
                      it->second[0]->DX() == 2.0 * m_octants[i]->DX()))
                {
                    valid = false;
                    cout << " 1 neigbour error" << endl;
                    cout << m_octants[i]->DX() << " " << it->second[0]->DX()
                         << endl;
                }
            }
            else if (it->second.size() == 4)
            {
                if (!(m_octants[i]->DX() / 2.0 == it->second[0]->DX()))
                {
                    valid = false;
                    cout << "4 neibour error" << endl;
                    cout << m_octants[i]->DX() << " " << it->second[0]->DX()
                         << endl;
                }
            }
        }
        if (!valid)
        {
            error = true;
            cout << "invalid neigbour config" << endl;
        }
    }
    return !error;
}

void Octree::SmoothSurfaceOctants()
{
    // for all the m_octants which are surface containing and know their delta
    // specification, look over all neighbours and ensure the specification
    // between them is smooth
    int ct = 0;

    do
    {
        ct = 0;
        for (int i = 0; i < m_octants.size(); i++)
        {
            OctantSharedPtr oct = m_octants[i];

            if (oct->IsDeltaKnown())
            {
                vector<OctantSharedPtr> check;
                map<OctantFace, vector<OctantSharedPtr> > nList =
                    oct->GetNeigbours();
                map<OctantFace, vector<OctantSharedPtr> >::iterator it;

                for (it = nList.begin(); it != nList.end(); it++)
                {
                    for (int j = 0; j < it->second.size(); j++)
                    {
                        if (it->second[j]->IsDeltaKnown() &&
                            it->second[j]->GetDelta() < oct->GetDelta() &&
                            ddx(oct, it->second[j]) > 0.2)
                        {
                            check.push_back(it->second[j]);
                        }
                    }
                }

                // for each neighbour listed in check_id, figure out the
                // smoothed
                // delta, and asign the miminum of these to nodes[i].GetDelta()
                if (check.size() > 0)
                {
                    NekDouble deltaSM = numeric_limits<double>::max();
                    for (int j = 0; j < check.size(); j++)
                    {
                        NekDouble r = oct->Distance(check[j]);

                        if (0.199 * r + check[j]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.199 * r + check[j]->GetDelta();
                        }
                    }
                    oct->SetDelta(deltaSM);
                    ct += 1;
                }
            }
        }
    } while (ct > 0);
}

void Octree::PropagateDomain()
{
    // until all m_octants know their delta specifcation and orientaion
    // look over all m_octants and if their neighours know either their
    // orientation
    // or specifcation calculate one for this octant
    int ct = 0;

    do
    {
        ct = 0;
        for (int i = 0; i < m_octants.size(); i++)
        {
            OctantSharedPtr oct = m_octants[i];

            if (!oct->IsDeltaKnown())
            { // if delta has not been asigned
                vector<OctantSharedPtr> known;
                map<OctantFace, vector<OctantSharedPtr> > nList =
                    oct->GetNeigbours();
                map<OctantFace, vector<OctantSharedPtr> >::iterator it;

                for (it = nList.begin(); it != nList.end(); it++)
                {
                    for (int j = 0; j < it->second.size(); j++)
                    {
                        if (it->second[j]->IsDeltaKnown())
                        {
                            known.push_back(it->second[j]);
                        }
                    }
                }

                if (known.size() > 0)
                {
                    vector<NekDouble> deltaPrime;
                    for (int j = 0; j < known.size(); j++)
                    {
                        NekDouble r = oct->Distance(known[j]);

                        if (0.199 * r + known[j]->GetDelta() < m_maxDelta)
                        {
                            deltaPrime.push_back(0.199 * r +
                                                 known[j]->GetDelta());
                        }
                        else
                        {
                            deltaPrime.push_back(m_maxDelta);
                        }
                    }
                    NekDouble min = numeric_limits<double>::max();
                    for (int j = 0; j < deltaPrime.size(); j++)
                    {
                        if (deltaPrime[j] < min)
                        {
                            min = deltaPrime[j];
                        }
                    }
                    oct->SetDelta(min);
                    ct += 1;
                    deltaPrime.clear();
                }
                known.clear();
            }

            if (oct->GetLocation() == eUnknown)
            { // if the node does not know its location
                vector<OctantSharedPtr> known;
                map<OctantFace, vector<OctantSharedPtr> > nList =
                    oct->GetNeigbours();
                map<OctantFace, vector<OctantSharedPtr> >::iterator it;

                for (it = nList.begin(); it != nList.end(); it++)
                {
                    for (int j = 0; j < it->second.size(); j++)
                    {
                        if (it->second[j]->GetLocation() != eUnknown)
                        {
                            known.push_back(it->second[j]);
                        }
                    }
                }

                if (known.size() > 0)
                {
                    vector<OctantSharedPtr> isNotOnBound;
                    for (int j = 0; j < known.size(); j++)
                    {
                        if (known[j]->GetLocation() != eOnBoundary)
                        {
                            isNotOnBound.push_back(known[j]);
                        }
                    }

                    if (isNotOnBound.size() > 0)
                    {
                        oct->SetLocation(isNotOnBound[0]->GetLocation());
                    }
                    else
                    {
                        NekDouble dist = numeric_limits<double>::max();

                        OctantSharedPtr closest;

                        bool f = false;
                        for (int j = 0; j < known.size(); j++)
                        {
                            if (oct->Distance(known[j]) < dist)
                            {
                                closest = known[j];
                                dist    = oct->Distance(known[j]);
                                f       = true;
                            }
                        }
                        ASSERTL0(f, "closest never set");

                        SPBaseSharedPtr sp = closest->GetABoundPoint();

                        Array<OneD, NekDouble> octloc, sploc, vec(3), uv, N;
                        int surf;
                        sp->GetCAD(surf, uv);
                        N = m_mesh->m_cad->GetSurf(surf)->N(uv);

                        octloc = oct->GetLoc();
                        sploc  = sp->GetLoc();

                        vec[0] = octloc[0] - sploc[0];
                        vec[1] = octloc[1] - sploc[1];
                        vec[2] = octloc[2] - sploc[2];

                        NekDouble dot =
                            vec[0] * N[0] + vec[1] * N[1] + vec[2] * N[2];

                        if (dot <= 0.0)
                        {
                            oct->SetLocation(eOutside);
                            ct += 1;
                        }
                        else
                        {
                            oct->SetLocation(eInside);
                            ct += 1;
                        }
                    }
                }
                known.clear();
            }
        }

    } while (ct > 0);

    for (int i = 0; i < m_octants.size(); i++)
    {
        if (!m_octants[i]->IsDeltaKnown())
        {
            m_octants[i]->SetDelta(m_maxDelta);
        }
    }
}

void Octree::SmoothAllOctants()
{
    // until no more changes occur smooth the mesh specification between all
    // m_octants not particualrly strictly
    int ct = 0;

    do
    {
        ct = 0;
        for (int i = 0; i < m_octants.size(); i++)
        {
            OctantSharedPtr oct = m_octants[i];

            vector<OctantSharedPtr> check;
            map<OctantFace, vector<OctantSharedPtr> > nList =
                oct->GetNeigbours();
            map<OctantFace, vector<OctantSharedPtr> >::iterator it;
            for (it = nList.begin(); it != nList.end(); it++)
            {
                for (int j = 0; j < it->second.size(); j++)
                {
                    if (it->second[j]->GetDelta() < oct->GetDelta() &&
                        ddx(oct, it->second[j]) > 0.2)
                    {
                        check.push_back(it->second[j]);
                    }
                }
            }

            if (check.size() > 0)
            {
                NekDouble deltaSM = numeric_limits<double>::max();
                for (int j = 0; j < check.size(); j++)
                {
                    NekDouble r = oct->Distance(check[j]);

                    if (0.199 * r + check[j]->GetDelta() < deltaSM)
                    {
                        deltaSM = 0.199 * r + check[j]->GetDelta();
                    }
                }
                oct->SetDelta(deltaSM);
                ct += 1;
            }
        }

    } while (ct > 0);
}

int Octree::CountElemt()
{
    // by considering the volume of a tet evaluate the number of elements in the
    // mesh

    NekDouble total = 0.0;

    Array<OneD, NekDouble> boundingBox = m_mesh->m_cad->GetBoundingBox();

    for (int i = 0; i < m_octants.size(); i++)
    {
        OctantSharedPtr oct = m_octants[i];
        if (oct->GetLocation() == eInside)
        {
            total += 8.0 * oct->DX() * oct->DX() * oct->DX() /
                     (oct->GetDelta() * oct->GetDelta() * oct->GetDelta() /
                      6.0 / sqrt(2));
        }
        else if (oct->GetLocation() == eOnBoundary)
        {
            NekDouble vol = 1.0;
            if (oct->FX(eBack) < boundingBox[1] &&
                oct->FX(eForward) > boundingBox[0])
            {
                // then there is some over lap in x
                NekDouble min, max;
                if (oct->FX(eBack) > boundingBox[0])
                {
                    min = oct->FX(eBack);
                }
                else
                {
                    min = boundingBox[0];
                }
                if (boundingBox[1] < oct->FX(eForward))
                {
                    max = boundingBox[1];
                }
                else
                {
                    max = oct->FX(eForward);
                }
                vol *= (max - min);
            }
            else
            {
                vol *= 0.0;
            }

            if (oct->FX(eDown) < boundingBox[3] &&
                oct->FX(eUp) > boundingBox[2])
            {
                // then there is some over lap in x
                NekDouble min, max;
                if (oct->FX(eDown) > boundingBox[2])
                {
                    min = oct->FX(eDown);
                }
                else
                {
                    min = boundingBox[2];
                }
                if (boundingBox[3] < oct->FX(eUp))
                {
                    max = boundingBox[3];
                }
                else
                {
                    max = oct->FX(eUp);
                }
                vol *= (max - min);
            }
            else
            {
                vol *= 0.0;
            }

            if (oct->FX(eRight) < boundingBox[5] &&
                oct->FX(eLeft) > boundingBox[4])
            {
                // then there is some over lap in x
                NekDouble min, max;
                if (oct->FX(eRight) > boundingBox[4])
                {
                    min = oct->FX(eRight);
                }
                else
                {
                    min = boundingBox[4];
                }
                if (boundingBox[5] < oct->FX(eLeft))
                {
                    max = boundingBox[5];
                }
                else
                {
                    max = oct->FX(eLeft);
                }
                vol *= (max - min);
            }
            else
            {
                vol *= 0.0;
            }
            total += vol / 2.0 / (oct->GetDelta() * oct->GetDelta() *
                                  oct->GetDelta() / 6.0 / sqrt(2));
        }
    }

    return int(total);
}

void Octree::CompileSourcePointList()
{
    int totalEnt = 0;
    if(m_mesh->m_cad->Is2D())
    {
        totalEnt += m_mesh->m_cad->GetNumCurve();
        for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
        {
            if (m_mesh->m_verbose)
            {
                LibUtilities::PrintProgressbar(i, totalEnt,
                                               "\tCompiling source points");
            }

            CADCurveSharedPtr curve = m_mesh->m_cad->GetCurve(i);
            Array<OneD, NekDouble> bds = curve->GetBounds();
            //this works assuming the curves are not distorted
            int samples  = ceil(curve->Length(bds[0],bds[1]) / m_minDelta) * 2;
            samples = max(40, samples);
            NekDouble dt = (bds[1] - bds[0]) / (samples + 1);
            for (int j = 1; j < samples - 1; j++) // dont want first and last point
            {
                NekDouble t = bds[0] + dt * j;
                NekDouble C = curve->Curvature(t);

                Array<OneD, NekDouble> loc = curve->P(t);

                vector<pair<weak_ptr<CADSurf>, CADOrientation::Orientation> > ss =
                    curve->GetAdjSurf();
                Array<OneD, NekDouble> uv = ss[0].first.lock()->locuv(loc);

                if (C != 0.0)
                {
                    NekDouble del = 2.0 * (1.0 / C) * sqrt(m_eps * (2.0 - m_eps));

                    if (del > m_maxDelta)
                    {
                        del = m_maxDelta;
                    }
                    if (del < m_minDelta)
                    {
                        del = m_minDelta;
                    }

                    CPointSharedPtr newCPoint =
                        MemoryManager<CPoint>::AllocateSharedPtr(
                            ss[0].first.lock()->GetId(), uv, loc, del);

                    m_SPList.push_back(newCPoint);
                }
                else
                {
                    BPointSharedPtr newBPoint =
                        MemoryManager<BPoint>::AllocateSharedPtr(
                            ss[0].first.lock()->GetId(), uv, loc);

                    m_SPList.push_back(newBPoint);
                }
            }
        }
    }
    else
    {
        totalEnt = m_mesh->m_cad->GetNumSurf();
        for (int i = 1; i <= totalEnt; i++)
        {
            if (m_mesh->m_verbose)
            {
                LibUtilities::PrintProgressbar(i, totalEnt,
                                               "\tCompiling source points");
            }

            CADSurfSharedPtr surf = m_mesh->m_cad->GetSurf(i);
            Array<OneD, NekDouble> bounds = surf->GetBounds();

            // to figure out the amount of curvature sampling to be conducted on
            // each parameter plane the surface is first sampled with a 40x40
            // grid the real space lengths of this grid are analyzed to find the
            // largest stretching in the u and v directions
            // this stretching is then considered with the mindelta user input
            // to find a number of sampling points in each direction which
            // ensures that in the final octree each surface octant will have at
            // least one sample point within its volume.
            // the 40x40 grid is used to ensure each surface has a minimum of
            // 40x40 samples.
            NekDouble du = (bounds[1] - bounds[0]) / (40 - 1);
            NekDouble dv = (bounds[3] - bounds[2]) / (40 - 1);

            NekDouble DeltaU = 0.0;
            NekDouble DeltaV = 0.0;

            Array<TwoD, Array<OneD, NekDouble> > samplepoints(40, 40);

            for (int j = 0; j < 40; j++)
            {
                for (int k = 0; k < 40; k++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = k * du + bounds[0];
                    uv[1] = j * dv + bounds[2];
                    if (j == 40 - 1)
                        uv[1] = bounds[3];
                    if (k == 40 - 1)
                        uv[0]          = bounds[1];
                    samplepoints[k][j] = surf->P(uv);
                }
            }

            for (int j = 0; j < 40 - 1; j++)
            {
                for (int k = 0; k < 40 - 1; k++)
                {
                    NekDouble deltau = sqrt(
                        (samplepoints[k][j][0] - samplepoints[k + 1][j][0]) *
                        (samplepoints[k][j][0] - samplepoints[k + 1][j][0]) +
                        (samplepoints[k][j][1] - samplepoints[k + 1][j][1]) *
                        (samplepoints[k][j][1] - samplepoints[k + 1][j][1]) +
                        (samplepoints[k][j][2] - samplepoints[k + 1][j][2]) *
                        (samplepoints[k][j][2] - samplepoints[k + 1][j][2]));
                    NekDouble deltav = sqrt(
                        (samplepoints[k][j][0] - samplepoints[k][j + 1][0]) *
                        (samplepoints[k][j][0] - samplepoints[k][j + 1][0]) +
                        (samplepoints[k][j][1] - samplepoints[k][j + 1][1]) *
                        (samplepoints[k][j][1] - samplepoints[k][j + 1][1]) +
                        (samplepoints[k][j][2] - samplepoints[k][j + 1][2]) *
                        (samplepoints[k][j][2] - samplepoints[k][j + 1][2]));

                    if (deltau > DeltaU)
                        DeltaU = deltau;
                    if (deltav > DeltaV)
                        DeltaV = deltav;
                }
            }

            // these are the acutal number of sample points in each parametric
            // direction
            int nu = ceil(DeltaU * 40 / m_minDelta) * 2;
            int nv = ceil(DeltaV * 40 / m_minDelta) * 2;
            nu = max(40, nu);
            nv = max(40, nv);

            for (int j = 0; j < nu; j++)
            {
                for (int k = 0; k < nv; k++)
                {
                    Array<OneD, NekDouble> uv(2);
                    uv[0] = (bounds[1] - bounds[0]) / (nu - 1) * j + bounds[0];
                    uv[1] = (bounds[3] - bounds[2]) / (nv - 1) * k + bounds[2];

                    // this prevents round off error at the end of the surface
                    // may not be neseercary but works
                    if (j == nu - 1)
                        uv[0] = bounds[1];
                    if (k == nv - 1)
                        uv[1] = bounds[3];

                    NekDouble C = surf->Curvature(uv);

                    // create new point based on smallest R, flat surfaces have
                    // k=0 but still need a point for element estimation
                    if (C == -1.0)
                    {
                        // Curvature not defined
                        continue;
                    }
                    if (C != 0.0)
                    {
                        NekDouble del =
                            2.0 * (1.0 / C) * sqrt(m_eps * (2.0 - m_eps));

                        if (del > m_maxDelta)
                        {
                            del = m_maxDelta;
                        }
                        if (del < m_minDelta)
                        {
                            del = m_minDelta;
                        }

                        CPointSharedPtr newCPoint =
                            MemoryManager<CPoint>::AllocateSharedPtr(
                                surf->GetId(), uv, surf->P(uv), del);

                        m_SPList.push_back(newCPoint);
                    }
                    else
                    {
                        BPointSharedPtr newBPoint =
                            MemoryManager<BPoint>::AllocateSharedPtr(
                                surf->GetId(), uv, surf->P(uv));

                        m_SPList.push_back(newBPoint);
                    }
                }
            }
        }
    }
    if (m_mesh->m_verbose)
    {
        cout << endl;
    }

    if (m_refinement.size() > 0)
    {
        if (m_mesh->m_verbose)
        {
            cout << "\t\tModifying based on refinement lines" << endl;
        }
        // now deal with the user defined spacing
        vector<string> lines;

        boost::split(lines, m_refinement, boost::is_any_of(":"));

        for (int i = 0; i < lines.size(); i++)
        {
            vector<NekDouble> data;
            ParseUtils::GenerateVector(lines[i], data);

            Array<OneD, NekDouble> x1(3), x2(3);

            x1[0] = data[0];
            x1[1] = data[1];
            x1[2] = data[2];
            x2[0] = data[3];
            x2[1] = data[4];
            x2[2] = data[5];

            m_lsources.push_back(linesource(x1, x2, data[6], data[7]));
        }

        // this takes any existing sourcepoints within the influence range
        // and modifies them
        /*for (int i = 0; i < m_SPList.size(); i++)
        {
            for (int j = 0; j < m_lsources.size(); j++)
            {
                if (m_lsources[j].withinRange(m_SPList[i]->GetLoc()))
                {
                    if(m_SPList[i]->GetType() == ePBoundary)
                    {
                        BPointSharedPtr bp =
                            std::dynamic_pointer_cast<BPoint>
                                                            (m_SPList[i]);

                        m_SPList[i] = bp->ChangeType();

                    }
                    m_SPList[i]->SetDelta(m_lsources[j].delta);
                }
            }
        }*/
        /// TODO add extra source points from the line souce to the octree
    }
}

NekDouble Octree::ddx(OctantSharedPtr i, OctantSharedPtr j)
{
    return fabs(i->GetDelta() - j->GetDelta()) / i->Distance(j);
}
}
}
