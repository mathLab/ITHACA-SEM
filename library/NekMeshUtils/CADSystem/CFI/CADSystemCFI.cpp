////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
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

#include "CADSystemCFI.h"
#include "CADVertCFI.h"
#include "CADCurveCFI.h"
#include "CADSurfCFI.h"

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{


std::string CADSystemCFI::key = GetEngineFactory().RegisterCreatorFunction(
        "cfi", CADSystemCFI::create, "Uses CFI as cad engine");

bool CADSystemCFI::LoadCAD()
{
    cout << "trying " << m_name << endl;

    cfiHandel.startServer();
    cout << "cfi loaded in mode: ";
    if(cfiHandel.info.mode == cfi::MODE_STANDALONE)
    {
        cout << "standalone" << endl;
    }
    else if(cfiHandel.info.mode == cfi::MODE_CLIENT)
    {
        cout << "client" << endl;
    }
    else if(cfiHandel.info.mode == cfi::MODE_SERVER)
    {
        cout << "server" << endl;
    }
    else if(cfiHandel.info.mode == cfi::MODE_BOTH)
    {
        cout << "both" << endl;
    }
    else if(cfiHandel.info.mode == cfi::MODE_PLUGIN)
    {
        cout << "plugin" << endl;
    }
    else
    {
        cout << "unknown" << endl;
    }

    model = cfiHandel.openModelFile(m_name.c_str());

    ASSERTL0(model->getEntityTotal(cfi::TYPE_BODY,cfi::SUBTYPE_ALL) == 1,
                    "cannot deal with multibodies");

    body = static_cast<cfi::Body*>(model->getEntity("W1"));

    map<string,cfi::Point*> mapOfVerts;
    map<string,cfi::Line*> mapOfEdges;

    vector< cfi::Oriented<cfi::TopoEntity*> >* faceList = body->getChildList();

    vector< cfi::Oriented<cfi::TopoEntity*> >::iterator it,it2,it3;
    for(it = faceList->begin(); it != faceList->end(); it++)
    {
        cfi::Oriented<cfi::TopoEntity*> orientatedFace = *it;
        cfi::Face* face = static_cast<cfi::Face*>(orientatedFace.entity);

        vector< cfi::Oriented<cfi::TopoEntity*> >* edgeList = face->getChildList();
        for(it2 = edgeList->begin(); it2 != edgeList->end(); it2++)
        {
            cfi::Oriented<cfi::TopoEntity*> orientatedEdge = *it2;
            cfi::Line* edge = static_cast<cfi::Line*>(orientatedEdge.entity);
            mapOfEdges[edge->getName()] = edge;

            vector< cfi::Oriented<cfi::TopoEntity*> >* vertList = edge->getChildList();
            for(it3 = vertList->begin(); it3 != vertList->end(); it3++)
            {
                cfi::Oriented<cfi::TopoEntity*> orientatedVert = *it3;
                cfi::Point* vert = static_cast<cfi::Point*>(orientatedVert.entity);
                mapOfVerts[vert->getName()] = vert;
                mapVertToListEdge[vert->getName()].push_back(edge->getName());
            }
        }
    }

    map<string,int> nameToVertId;
    map<string,cfi::Point*>::iterator vit;
    int i = 1; //from one to be consistent with oce
    for(vit = mapOfVerts.begin(); vit != mapOfVerts.end(); vit++, i++)
    {
        AddVert(i, vit->second);
        nameToVertId[vit->second->getName()] = i;
    }

    map<string,cfi::Line*>::iterator eit;
    i = 1;
    for(eit = mapOfEdges.begin(); eit != mapOfEdges.end(); eit++, i++)
    {
        vector< cfi::Oriented<cfi::TopoEntity*> >* vertList = eit->second->getChildList();
        vector<int> ids;
        for(it3 = vertList->begin(); it3 != vertList->end(); it3++)
        {
            cfi::Oriented<cfi::TopoEntity*> orientatedVert = *it3;
            cfi::Point* vert = static_cast<cfi::Point*>(orientatedVert.entity);
            ids.push_back(nameToVertId[vert->getName()]);
        }
        ASSERTL0(ids.size()==2,"doesnt make sense");
        nameToCurveId[eit->second->getName()] = i;
        AddCurve(i, eit->second, ids[0], ids[1]);
    }

    map<int, vector<int> > adjsurfmap;
    i = 1;
    for(it = faceList->begin(); it != faceList->end(); it++, i++)
    {
        cfi::Oriented<cfi::TopoEntity*> orientatedFace = *it;
        cfi::Face* face = static_cast<cfi::Face*>(orientatedFace.entity);
        nameToFaceId[face->getName()] = i;
        vector< cfi::Oriented<cfi::TopoEntity*> >* edgeList = face->getChildList();

        vector<EdgeLoop> edgeloops;
        vector<vector<cfi::Oriented<cfi::TopoEntity*> > > cfiloops;
        int done = 0;
        while(done != edgeList->size())
        {
            EdgeLoop edgeloop;
            vector<cfi::Oriented<cfi::TopoEntity*> > cfiloop;
            string firstVert;
            vector< cfi::Oriented<cfi::TopoEntity*> >* vertList = edgeList->at(done).entity->getChildList();
            if(edgeList->at(done).orientation == 1)
            {
                firstVert = vertList->at(0).entity->getName();

            }
            else
            {
                firstVert = vertList->at(1).entity->getName();
            }
            edgeloop.edges.push_back(m_curves[nameToCurveId[edgeList->at(done).entity->getName()]]);
            cfiloop.push_back(edgeList->at(done));
            adjsurfmap[nameToCurveId[edgeList->at(done).entity->getName()]].push_back(i);
            edgeList->at(done).orientation == 1 ? edgeloop.edgeo.push_back(0) : edgeloop.edgeo.push_back(1);

            for(done++; done < edgeList->size(); done++)
            {
                bool end = false;
                vertList = edgeList->at(done).entity->getChildList();
                if(edgeList->at(done).orientation == 1)
                {
                    if(vertList->at(1).entity->getName() == firstVert)
                    {
                        end = true;
                    }
                }
                else
                {
                    if(vertList->at(0).entity->getName() == firstVert)
                    {
                        end = true;
                    }
                }

                edgeloop.edges.push_back(m_curves[nameToCurveId[edgeList->at(done).entity->getName()]]);
                cfiloop.push_back(edgeList->at(done));
                adjsurfmap[nameToCurveId[edgeList->at(done).entity->getName()]].push_back(i);
                edgeList->at(done).orientation == 1 ? edgeloop.edgeo.push_back(0) : edgeloop.edgeo.push_back(1);

                if(end)
                {
                    done++;
                    break;
                }
            }
            cfiloops.push_back(cfiloop);
            edgeloops.push_back(edgeloop);
        }

        cfi::FaceMassProperties* fmp = face->calcMassProperties(5.0,5.0);
        cfi::UVPosition uv = face->calcUVFromXYZ(fmp->centreOfMass);
        edgeloops[0].center = Array<OneD, NekDouble>(2);
        edgeloops[0].center[0] = uv.u;
        edgeloops[0].center[1] = uv.v;
        for(int i = 1; i < cfiloops.size(); i++)
        {
            for(int j = 0; j < cfiloops[i].size(); j++)
            {
                cfiloops[i][j].orientation == cfi::ORIENT_POSITIVE ?
                            cfiloops[i][j].orientation = cfi::ORIENT_NEGATIVE :
                            cfiloops[i][j].orientation = cfi::ORIENT_POSITIVE;
            }
            cfi::Face* tmpface = cfi::Face::createBasic(model, cfiloops[i],face->getTopoEmbedding().value());
            cfi::FaceMassProperties* fmp = tmpface->calcMassProperties(5.0,5.0);
            cfi::UVPosition uv = tmpface->calcUVFromXYZ(fmp->centreOfMass);
            edgeloops[i].center = Array<OneD, NekDouble>(2);
            edgeloops[i].center[0] = uv.u;
            edgeloops[i].center[1] = uv.v;
        }
        AddSurf(i, face, edgeloops);
    }

    //TODO identify Degenerated faces and setdegen on vertices accordinaly

    ASSERTL0(adjsurfmap.size() == m_curves.size(),"incorrect curve info");

    // This checks that all edges are bound by two surfaces, sanity check.
    for (map<int, vector<int> >::iterator it = adjsurfmap.begin();
         it != adjsurfmap.end(); it++)
    {
        ASSERTL0(it->second.size() == 2, "no three curve surfaces");
        vector<CADSurfSharedPtr> sfs;
        for (int i = 0; i < it->second.size(); i++)
        {
            sfs.push_back(m_surfs[it->second[i]]);
        }
        m_curves[it->first]->SetAdjSurf(sfs);
    }

    return true;
}

void CADSystemCFI::AddVert(int i, cfi::Point* in)
{
    CADVertSharedPtr newVert = GetCADVertFactory().CreateInstance(key);

    static_pointer_cast<CADVertCFI>(newVert)->Initialise(i,in);

    m_verts[i] = newVert;
}

void CADSystemCFI::AddCurve(int i, cfi::Line* in, int fv, int lv)
{
    CADCurveSharedPtr newCurve = GetCADCurveFactory().CreateInstance(key);
    static_pointer_cast<CADCurveCFI>(newCurve)->Initialise(i, in);

    vector<CADVertSharedPtr> vs;
    vs.push_back(m_verts[fv]);
    vs.push_back(m_verts[lv]);
    m_curves[i] = newCurve;
    m_curves[i]->SetVert(vs);
}

void CADSystemCFI::AddSurf(int i, cfi::Face* in, std::vector<EdgeLoop> ein)
{
    CADSurfSharedPtr newSurf = GetCADSurfFactory().CreateInstance(key);
    static_pointer_cast<CADSurfCFI>(newSurf)->Initialise(i, in, ein);
    m_surfs[i] = newSurf;

    //if (in.Orientation() == 0)
    //{
    //    m_surfs[i]->SetReverseNomral();
    //}

    int tote = 0;
    for (int i = 0; i < ein.size(); i++)
    {
        tote += ein[i].edges.size();
    }

    ASSERTL0(tote != 1, "cannot handle periodic curves");
}

Array<OneD, NekDouble> CADSystemCFI::GetBoundingBox()
{
    cfi::BoundingBox box = model->calcBoundingBox();

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
