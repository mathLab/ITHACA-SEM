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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

//#include <algorithm>
//#include <limits>

#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

NekDouble Octree::Query(Array<OneD, NekDouble> loc)
{
    // starting at master octant 0 move through succsesive m_octants which
    // contain
    // the point loc until a leaf is found
    OctantSharedPtr n = m_masteroct;
    int quad;

    bool found = false;

    while (!found)
    {
        Array<OneD, NekDouble> octloc = n->GetLoc();

        if (!(loc[0] < octloc[0]) && // forward
            !(loc[1] < octloc[1]) && // forward
            !(loc[2] < octloc[2])) // forward
        {
            quad = 0;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] > octloc[2])) // back
        {
            quad = 1;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] < octloc[2])) // forward
        {
            quad = 2;
        }
        else if (!(loc[0] < octloc[0]) && // forward
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] > octloc[2])) // back
        {
            quad = 3;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] < octloc[2])) // forward
        {
            quad = 4;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] < octloc[1]) && // forward
                 !(loc[2] > octloc[2])) // back
        {
            quad = 5;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] < octloc[2])) // forward
        {
            quad = 6;
        }
        else if (!(loc[0] > octloc[0]) && // back
                 !(loc[1] > octloc[1]) && // back
                 !(loc[2] > octloc[2])) // back
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
    return n->GetDelta();
}

void Octree::GetOctreeMesh(MeshSharedPtr m)
{
    for (int i = 0; i < m_octants.size(); i++)
    {
        /*if(m_octants[i]->GetLocation() != eOnBoundary)
        {
            continue;
        }*/

        vector<NodeSharedPtr> ns(8);

        ns[0] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eBack),
                                                 m_octants[i]->FX(eDown),
                                                 m_octants[i]->FX(eRight)));

        ns[1] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eForward),
                                                 m_octants[i]->FX(eDown),
                                                 m_octants[i]->FX(eRight)));

        ns[2] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eForward),
                                                 m_octants[i]->FX(eUp),
                                                 m_octants[i]->FX(eRight)));

        ns[3] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eBack),
                                                 m_octants[i]->FX(eUp),
                                                 m_octants[i]->FX(eRight)));

        ns[4] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eBack),
                                                 m_octants[i]->FX(eDown),
                                                 m_octants[i]->FX(eLeft)));

        ns[5] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eForward),
                                                 m_octants[i]->FX(eDown),
                                                 m_octants[i]->FX(eLeft)));

        ns[6] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eForward),
                                                 m_octants[i]->FX(eUp),
                                                 m_octants[i]->FX(eLeft)));

        ns[7] = boost::shared_ptr<Node>(new Node(0,
                                                 m_octants[i]->FX(eBack),
                                                 m_octants[i]->FX(eUp),
                                                 m_octants[i]->FX(eLeft)));

        vector<int> tags;
        tags.push_back(0);
        ElmtConfig conf(LibUtilities::eHexahedron, 1, false, false);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eHexahedron, conf, ns, tags);
        m->m_element[3].push_back(E);
    }
}

void Octree::Build()
{
    Array<OneD, NekDouble> boundingBox = m_cad->GetBoundingBox();

    if (m_verbose)
        cout << endl << "Octree system" << endl;

    // build curvature samples
    CompileCuravturePointList();

    if (m_verbose)
        cout << "\tCurvature samples: " << m_cpList.size() << endl;

    m_dim = max((boundingBox[1] - boundingBox[0]) / 2.0,
                (boundingBox[3] - boundingBox[2]) / 2.0);

    m_dim = max(m_dim, (boundingBox[5] - boundingBox[4]) / 2.0);

    m_centroid    = Array<OneD, NekDouble>(3);
    m_centroid[0] = (boundingBox[1] + boundingBox[0]) / 2.0;
    m_centroid[1] = (boundingBox[3] + boundingBox[2]) / 2.0;
    m_centroid[2] = (boundingBox[5] + boundingBox[4]) / 2.0;

    // make master octant based on the bounding box of the domain
    m_masteroct = MemoryManager<Octant>::AllocateSharedPtr(
        0, m_centroid[0], m_centroid[1], m_centroid[2], m_dim, m_cpList);

    SubDivide();

    m_octants.clear();
    m_masteroct->CompileLeaves(m_octants);

    if (m_verbose)
    {
        cout << "\tOctants: " << m_octants.size() << endl;
    }

    SmoothSurfaceOctants();

    PropagateDomain();

    SmoothAllOctants();

    if (m_verbose)
    {
        int elem = CountElemt();
        cout << "\tPredicted mesh: " << elem << endl;
    }
}

void Octree::SubDivide()
{
    bool repeat;
    int ct = 1;

    m_numoct = 1;
    m_masteroct->Subdivide(m_masteroct, m_numoct);

    if (m_verbose)
    {
        cout << "\tSubdivide iteration: ";
    }

    do
    {
        ct++;
        if (m_verbose)
        {
            cout << ct << " ";
            cout.flush();
        }
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

    if (m_verbose)
    {
        cout << endl;
    }
}

bool Octree::VerifyNeigbours()
{
    // check all neibours
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
                NekDouble expectedfx;
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
                            ddx(oct, it->second[j]) > 0.1)
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

                        if (0.099 * r + check[j]->GetDelta() < deltaSM)
                        {
                            deltaSM = 0.099 * r + check[j]->GetDelta();
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

                        if (0.14 * r + known[j]->GetDelta() < m_maxDelta)
                        {
                            deltaPrime.push_back(0.14 * r +
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

                        for (int j = 0; j < known.size(); j++)
                        {
                            if (oct->Distance(known[j]) < dist)
                            {
                                closest = known[j];
                                dist    = oct->Distance(known[j]);
                            }
                        }

                        CurvaturePointSharedPtr cp = closest->GetCPPoint();

                        Array<OneD, NekDouble> octloc, cploc, vec(3), uv, N;
                        int surf;
                        cp->GetCAD(surf, uv);
                        N = m_cad->GetSurf(surf)->N(uv);

                        octloc = oct->GetLoc();
                        cploc  = cp->GetLoc();

                        vec[0] = octloc[0] - cploc[0];
                        vec[1] = octloc[1] - cploc[1];
                        vec[2] = octloc[2] - cploc[2];

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
        ASSERTL0(m_octants[i]->IsDeltaKnown(),
                 "does not know delta after propergation");
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

    Array<OneD, NekDouble> boundingBox = m_cad->GetBoundingBox();

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

struct linesource
{
    Array<OneD, NekDouble> x1, x2;
    NekDouble R, delta;
    linesource(Array<OneD, NekDouble> p1,
               Array<OneD, NekDouble> p2,
               NekDouble r,
               NekDouble d)
        : x1(p1), x2(p2), R(r), delta(d)
    {
    }

    bool withinRange(Array<OneD, NekDouble> p)
    {
        Array<OneD, NekDouble> Le(3), Re(3), s(3);
        for (int i = 0; i < 3; i++)
        {
            Le[i] = p[i] - x1[i];
            Re[i] = p[i] - x2[i];
            s[i]  = x2[i] - x1[i];
        }
        Array<OneD, NekDouble> dev(3);
        dev[0] = Le[1] * Re[2] - Re[1] * Le[2];
        dev[1] = Le[0] * Re[2] - Re[0] * Le[2];
        dev[2] = Le[0] * Re[1] - Re[0] * Le[1];

        NekDouble dist =
            sqrt(dev[0] * dev[0] + dev[1] * dev[1] + dev[2] * dev[2]) /
            sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);

        NekDouble t = -1.0 * ((x1[0] - p[0]) * s[0] + (x1[1] - p[1]) * s[1] +
                              (x1[1] - p[1]) * s[1]) /
                      Length() / Length();

        if (dist < R && !(t > 1) && !(t < 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    NekDouble Length()
    {
        return sqrt((x1[0] - x2[0]) * (x1[0] - x2[0]) +
                    (x1[1] - x2[1]) * (x1[1] - x2[1]) +
                    (x1[2] - x2[2]) * (x1[2] - x2[2]));
    }
};

void Octree::CompileCuravturePointList()
{
    for (int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        CADSurfSharedPtr surf         = m_cad->GetSurf(i);
        Array<OneD, NekDouble> bounds = surf->GetBounds();

        // to figure out the amount of curvature sampling to be conducted on
        // each parameter plane the surface is first sampled with a 40x40 grid
        // the real space lengths of this grid are analysed to find the largest
        // strecthing in the u and v directions
        // this stretching this then cosnidered with the mindelta user input
        // to find a number of sampling points in each direction which
        // enures that in the final octree each surface octant will have at
        // least
        // one sample point within its volume.
        // the 40x40 grid is used to ensure each surface has a minimum of 40x40
        // samples.
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
        int nu = ceil(DeltaU / m_minDelta) * 40;
        int nv = ceil(DeltaV / m_minDelta) * 40;

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

                // create new point based on smallest R, flat surfaces have k=0
                // but still need a point for element estimation
                if (C != 0.0)
                {
                    bool minlimited = false;
                    NekDouble ideal;

                    NekDouble del =
                        2.0 * (1.0 / C) * sqrt(m_eps * (2.0 - m_eps));

                    if (del > m_maxDelta)
                    {
                        del = m_maxDelta;
                    }
                    if (del < m_minDelta)
                    {
                        ideal      = del;
                        del        = m_minDelta;
                        minlimited = true;
                    }

                    if (minlimited)
                    {
                        CurvaturePointSharedPtr newCPoint =
                            MemoryManager<CurvaturePoint>::AllocateSharedPtr(
                                surf->GetId(), uv, surf->P(uv), del, ideal);

                        m_cpList.push_back(newCPoint);
                    }
                    else
                    {
                        CurvaturePointSharedPtr newCPoint =
                            MemoryManager<CurvaturePoint>::AllocateSharedPtr(
                                surf->GetId(), uv, surf->P(uv), del);

                        m_cpList.push_back(newCPoint);
                    }
                }
                else
                {
                    CurvaturePointSharedPtr newCPoint =
                        MemoryManager<CurvaturePoint>::AllocateSharedPtr(
                            surf->GetId(), uv, surf->P(uv));

                    m_cpList.push_back(newCPoint);
                }
            }
        }
    }

    if (m_udsfile == "N")
    {
        return;
    }

    // now deal with the user defined spacing
    vector<linesource> lsources;
    fstream fle;
    fle.open(m_udsfile.c_str());

    string fileline;

    while (!fle.eof())
    {
        getline(fle, fileline);
        stringstream s(fileline);
        string word;
        s >> word;
        if (word == "#")
        {
            continue;
        }

        Array<OneD, NekDouble> x1(3), x2(3);
        NekDouble r, d;
        x1[0] = boost::lexical_cast<double>(word);
        s >> x1[1] >> x1[2] >> x2[0] >> x2[1] >> x2[2] >> r >> d;

        lsources.push_back(linesource(x1, x2, r, d));
    }
    fle.close();

    for (int j = 0; j < lsources.size(); j++)
    {
        cout << lsources[j].x1[0] << " " << lsources[j].x1[1] << " "
             << lsources[j].x1[2] << endl;
        cout << lsources[j].x2[0] << " " << lsources[j].x2[1] << " "
             << lsources[j].x2[2] << endl;
        cout << lsources[j].Length() << endl;
    }

    int ct = 0;
    for (int i = 0; i < m_cpList.size(); i++)
    {
        for (int j = 0; j < lsources.size(); j++)
        {
            if (lsources[j].withinRange(m_cpList[i]->GetLoc()))
            {
                ct++;
                m_cpList[i]->SetDelta(lsources[j].delta);
            }
        }
    }
    cout << ct << endl;

    ///@TODO need to add curvature points with the false tag to make octree
    ///modification work
    // off surfaces

    /*for(int i = 0; i < lsources.size(); i++)
    {
        int nc; //number of point to add cicularly
        int nl; //number of point to add length

        nc = ceil(2.0*3.142*lsources[i].R / lsources[i].delta)*2;
        nl = ceil(lsources[i].Length() / lsources[i].delta)*2;

        NekDouble dr = lsources[i].Length() / nl;
        NekDouble dtheta = 2.0*3.142 / nc;

        for(int j = 0; j < nl; j++)
        {
            NekDouble len =
            for(int k = 0; k < nc; k++)
            {

            }
        }

    }*/
}

NekDouble Octree::ddx(OctantSharedPtr i, OctantSharedPtr j)
{
    return fabs(i->GetDelta() - j->GetDelta()) / i->Distance(j);
}
}
}
