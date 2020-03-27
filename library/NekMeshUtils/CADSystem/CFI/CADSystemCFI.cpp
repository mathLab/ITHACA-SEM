////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystemCFI.cpp
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

#include "CADSystemCFI.h"
#include "CADCurveCFI.h"
#include "CADSurfCFI.h"
#include "CADVertCFI.h"

#include <boost/lexical_cast.hpp>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSystemCFI::key = GetEngineFactory().RegisterCreatorFunction(
    "cfi", CADSystemCFI::create, "Uses CFI as cad engine");

bool CADSystemCFI::LoadCAD()
{
    // it is possible to get CFI to lock on to a open gui session
    // not sure it ever will with this code
    m_cfiHandle.startServer();
    if (m_verbose)
    {
        cout << "cfi loaded in mode: ";
        if (m_cfiHandle.info.mode == cfi::MODE_STANDALONE)
        {
            cout << "standalone" << endl;
        }
        else if (m_cfiHandle.info.mode == cfi::MODE_CLIENT)
        {
            cout << "client" << endl;
        }
        else if (m_cfiHandle.info.mode == cfi::MODE_SERVER)
        {
            cout << "server" << endl;
        }
        else if (m_cfiHandle.info.mode == cfi::MODE_BOTH)
        {
            cout << "both" << endl;
        }
        else if (m_cfiHandle.info.mode == cfi::MODE_PLUGIN)
        {
            cout << "plugin" << endl;
        }
        else
        {
            cout << "unknown" << endl;
        }

        cout << "\tVersion " << m_cfiHandle.info.version << endl
             << "\tfixno " << m_cfiHandle.info.fixno << endl
             << "\tubid " << m_cfiHandle.info.ubid << endl;
    }

    if (m_config.count("UseCFIMesh"))
    {
        m_useCFIMesh = boost::lexical_cast<bool>(m_config["UseCFIMesh"]);
    }

    m_model = m_cfiHandle.openModelFile(m_name.c_str());

    if (m_model->getEntityTotal(cfi::TYPE_BODY, cfi::SUBTYPE_ALL) != 1)
    {
        if (m_useCFIMesh)
        {
            if (m_verbose)
            {
                cout << "\tWill extract mesh and have multibodies" << endl;
            }
        }
        else
        {
            if (m_verbose)
            {
                cout << "\tHas multibodies and instructions to mesh, this is "
                        "not possible" << endl;
            }
            abort();
        }
    }

    vector<cfi::Entity *> *bds =
        m_model->getEntityList(cfi::TYPE_BODY, cfi::SUBTYPE_ALL);

    for (auto &i : *bds)
    {
        cfi::Body *b = static_cast<cfi::Body *>(i);
        m_bodies.push_back(b);
    }

    // cfi doesnt mind stupid units so this scales everything back to meters
    // which is what nekmesh assumes its in
    // the m_scal object is passed to all cad entities and scales any operation
    // before running it.
    m_scal = 1.0;
    if (m_model->getUnits() == cfi::UNIT_INCHES)
    {
        if (m_verbose)
        {
            cout << "\tModel is in inches, scaling accordingly" << endl;
        }
        m_scal = 0.0254;
    }
    else if (m_model->getUnits() == cfi::UNIT_MILLIMETERS ||
             m_model->getUnits() == cfi::UNIT_MILLIMETRES)
    {
        if (m_verbose)
        {
            cout << "\tModel is in mm, scaling accordingly" << endl;
        }
        m_scal = 1e-3;
    }

    // CFI does everything by string identifers
    // currently nekmesh cad system uses integer ids.
    // it really should use strings but doesnt currently
    map<string, cfi::Point *> mapOfVerts;
    map<string, cfi::Line *> mapOfEdges;
    map<string, cfi::Face *> mapOfFaces;

    // nothing is unique in cfi, there can list of verts that arnt used
    // or are listed twice. This block gets all the real unique vertices in the
    // cad by cascading down from the faces
    // also builds a list on unique edges in the process

    for (int i = 0; i < m_bodies.size(); i++)
    {
        // check that it is not a group of bodies
        if (m_bodies[i]->getTopoSubtype() == cfi::SUBTYPE_COMBINED)
        {
            continue;
        }

        vector<cfi::Oriented<cfi::TopoEntity *>> *faceList =
            m_bodies[i]->getChildList();

        vector<cfi::Oriented<cfi::TopoEntity *>>::iterator it, it2, it3;
        for (it = faceList->begin(); it != faceList->end(); it++)
        {
            cfi::Oriented<cfi::TopoEntity *> orientatedFace = *it;
            cfi::Face *face = static_cast<cfi::Face *>(orientatedFace.entity);

            try
            {
                boost::optional<cfi::Oriented<cfi::Surface *>> surf =
                    face->getTopoEmbedding();
            }
            catch (std::exception &e)
            {
                // dont want this face
                continue;
            }

            auto it = mapOfFaces.find(face->getName());
            if (it == mapOfFaces.end())
            {
                vector<cfi::Oriented<cfi::TopoEntity *>> *edgeList =
                    face->getChildList();

                // unrolling combined edges
                vector<cfi::Oriented<cfi::TopoEntity *>> fullEdgeList;
                for (it2 = edgeList->begin(); it2 != edgeList->end(); it2++)
                {
                    cfi::Line *edge = static_cast<cfi::Line *>(it2->entity);

                    if (edge->getTopoSubtype() == cfi::SUBTYPE_COMBINED)
                    {
                        vector<cfi::Oriented<cfi::TopoEntity *>> *subEdgeList =
                            edge->getChildList();

                        for (it3 = subEdgeList->begin();
                             it3 != subEdgeList->end(); it3++)
                        {
                            fullEdgeList.push_back(*it3);
                        }
                    }
                    else
                    {
                        fullEdgeList.push_back(*it2);
                    }
                }

                delete edgeList;

                for (it2 = fullEdgeList.begin(); it2 != fullEdgeList.end(); it2++)
                {
                    cfi::Oriented<cfi::TopoEntity *> orientatedEdge = *it2;
                    cfi::Line *edge =
                        static_cast<cfi::Line *>(orientatedEdge.entity);
                    mapOfEdges[edge->getName()] = edge;

                    vector<cfi::Oriented<cfi::TopoEntity *>> *vertList =
                        edge->getChildList();
                    for (it3 = vertList->begin(); it3 != vertList->end(); it3++)
                    {
                        cfi::Oriented<cfi::TopoEntity *> orientatedVert = *it3;
                        cfi::Point *vert =
                            static_cast<cfi::Point *>(orientatedVert.entity);
                        mapOfVerts[vert->getName()] = vert;
                        m_mapVertToListEdge[vert->getName()].push_back(
                            edge->getName());
                    }
                    delete vertList;
                }

                mapOfFaces[face->getName()] = face;
            }
        }
        delete faceList;
    }

    // make the vertices and build a map of name to id
    map<string, cfi::Point *>::iterator vit;
    int i = 1; // from one to be consistent with oce
    for (vit = mapOfVerts.begin(); vit != mapOfVerts.end(); vit++, i++)
    {
        AddVert(i, vit->second);
        m_nameToVertId[vit->second->getName()] = i;
    }

    // build curves
    map<string, cfi::Line *>::iterator eit;
    i = 1;
    for (eit = mapOfEdges.begin(); eit != mapOfEdges.end(); eit++, i++)
    {
        AddCurve(i, eit->second);
        m_nameToCurveId[eit->second->getName()] = i;
    }

    // build surfaces
    map<string, cfi::Face *>::iterator fit;
    i = 1;
    for (fit = mapOfFaces.begin(); fit != mapOfFaces.end(); fit++, i++)
    {
        m_nameToFaceId[fit->second->getName()] = i;

        AddSurf(i, fit->second);
    }

    // TODO identify Degenerated faces and setdegen on vertices accordinaly

    // This checks that all edges are bound by two surfaces, sanity check.
    if (!m_2d && !m_useCFIMesh)
    {
        map<int, CADCurveSharedPtr>::iterator it;
        for (it = m_curves.begin(); it != m_curves.end(); it++)
        {
            ASSERTL0(it->second->GetAdjSurf().size() == 2,
                     "curve is not joined to 2 surfaces");
        }
    }

    if (m_verbose)
    {
        Report();
    }

    // Tidy up
    delete bds;

    return true;
}

void CADSystemCFI::AddVert(int i, cfi::Point *in)
{
    CADVertSharedPtr newVert = GetCADVertFactory().CreateInstance(key);

    static_pointer_cast<CADVertCFI>(newVert)->Initialise(i, in, m_scal);

    m_verts[i] = newVert;
    m_verts[i]->SetName(in->getName());
}

void CADSystemCFI::AddCurve(int i, cfi::Line *in)
{
    CADCurveSharedPtr newCurve = GetCADCurveFactory().CreateInstance(key);
    static_pointer_cast<CADCurveCFI>(newCurve)->Initialise(i, in, m_scal);

    vector<cfi::Oriented<cfi::TopoEntity *>> *vertList = in->getChildList();

    ASSERTL0(vertList->size() == 2, "should be two ends");

    vector<cfi::Oriented<cfi::TopoEntity *>>::iterator it;
    vector<NekDouble> t;
    for (it = vertList->begin(); it != vertList->end(); it++)
    {
        cfi::Oriented<cfi::TopoEntity *> orientatedVert = *it;
        cfi::Point *vert = static_cast<cfi::Point *>(orientatedVert.entity);
        boost::optional<cfi::Projected<double>> pj =
            in->calcTFromXYZ(vert->getGeometry(), -1);
        t.push_back(pj.value().parameters);
    }
    ASSERTL0(t[0] < t[1], "weirdness");

    vector<CADVertSharedPtr> vs;
    vs.push_back(m_verts[m_nameToVertId[vertList->at(0).entity->getName()]]);
    vs.push_back(m_verts[m_nameToVertId[vertList->at(1).entity->getName()]]);
    m_curves[i] = newCurve;
    m_curves[i]->SetVert(vs);
    m_curves[i]->SetName(in->getName());
    m_verts[m_nameToVertId[vertList->at(0).entity->getName()]]->AddAdjCurve(
        m_curves[i]);
    m_verts[m_nameToVertId[vertList->at(1).entity->getName()]]->AddAdjCurve(
        m_curves[i]);
    delete vertList;
}

void CADSystemCFI::AddSurf(int i, cfi::Face *in)
{
    CADSurfSharedPtr newSurf = GetCADSurfFactory().CreateInstance(key);
    static_pointer_cast<CADSurfCFI>(newSurf)->Initialise(i, in, m_scal);

    vector<cfi::Oriented<cfi::TopoEntity *>> *edgeList = in->getChildList();

    // unrolling combined edges
    vector<cfi::Oriented<cfi::TopoEntity *>> fullEdgeList;
    vector<cfi::Oriented<cfi::TopoEntity *>>::iterator it2, it3;
    for (it2 = edgeList->begin(); it2 != edgeList->end(); it2++)
    {
        cfi::Line *edge = static_cast<cfi::Line *>(it2->entity);

        if (edge->getTopoSubtype() == cfi::SUBTYPE_COMBINED)
        {
            vector<cfi::Oriented<cfi::TopoEntity *>> *subEdgeList =
                edge->getChildList();

            for (it3 = subEdgeList->begin(); it3 != subEdgeList->end(); it3++)
            {
                fullEdgeList.push_back(*it3);
            }

            if (it2->orientation == cfi::ORIENT_NEGATIVE)
            {
                reverse(fullEdgeList.begin() + fullEdgeList.size() -
                            subEdgeList->size(),
                        fullEdgeList.end());
            }
        }
        else
        {
            fullEdgeList.push_back(*it2);
        }
    }
    delete edgeList;

    vector<EdgeLoopSharedPtr> edgeloops;
    int done = 0;
    while (done != fullEdgeList.size())
    {
        EdgeLoopSharedPtr edgeloop = EdgeLoopSharedPtr(new EdgeLoop);
        string firstVert;
        vector<cfi::Oriented<cfi::TopoEntity *>> *vertList =
            fullEdgeList.at(done).entity->getChildList();
        if (fullEdgeList.at(done).orientation == cfi::ORIENT_POSITIVE)
        {
            firstVert = vertList->at(0).entity->getName();
            edgeloop->edgeo.push_back(CADOrientation::eForwards);
        }
        else
        {
            firstVert = vertList->at(1).entity->getName();
            edgeloop->edgeo.push_back(CADOrientation::eBackwards);
        }

        delete vertList;

        edgeloop->edges.push_back(
            m_curves[m_nameToCurveId[fullEdgeList.at(done).entity->getName()]]);

        for (done++; done < fullEdgeList.size(); done++)
        {
            bool end = false;
            vertList = fullEdgeList.at(done).entity->getChildList();
            if (fullEdgeList.at(done).orientation == cfi::ORIENT_POSITIVE)
            {
                if (vertList->at(1).entity->getName() == firstVert)
                {
                    end = true;
                }
                edgeloop->edgeo.push_back(CADOrientation::eForwards);
            }
            else
            {
                if (vertList->at(0).entity->getName() == firstVert)
                {
                    end = true;
                }
                edgeloop->edgeo.push_back(CADOrientation::eBackwards);
            }

            delete vertList;

            edgeloop->edges.push_back(
                m_curves[m_nameToCurveId[
                        fullEdgeList.at(done).entity->getName()]]);

            if (end)
            {
                done++;
                break;
            }
        }
        edgeloops.push_back(edgeloop);
    }

    // TODO find if surface has reversed normal or not

    int tote = 0;
    for (int k = 0; k < edgeloops.size(); k++)
    {
        tote += edgeloops[k]->edges.size();
    }

    ASSERTL0(tote != 1, "cannot handle periodic curves");

    CADSurf::OrientateEdges(newSurf, edgeloops);
    newSurf->SetEdges(edgeloops);

    // now the loops are orientated, tell the curves how they are
    for (int k = 0; k < edgeloops.size(); k++)
    {
        for (int j = 0; j < edgeloops[k]->edges.size(); j++)
        {
            edgeloops[k]->edges[j]->SetAdjSurf(
                make_pair(newSurf, edgeloops[k]->edgeo[j]));
        }
    }

    m_surfs[i] = newSurf;
    m_surfs[i]->SetName(in->getName());
}

Array<OneD, NekDouble> CADSystemCFI::GetBoundingBox()
{
    cfi::BoundingBox box = m_model->calcBoundingBox();

    Array<OneD, NekDouble> ret(6);
    ret[0] = box.xLower;
    ret[1] = box.xUpper;
    ret[2] = box.yLower;
    ret[3] = box.yUpper;
    ret[4] = box.zLower;
    ret[5] = box.zUpper;

    return ret;
}
}
}
