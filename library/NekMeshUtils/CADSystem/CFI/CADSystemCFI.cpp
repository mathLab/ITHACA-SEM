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

    cfi::Cfi cfi;
    model = cfi.openModelFile(m_name.c_str());

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
            }
        }
    }

    cout << "verts " << mapOfVerts.size() << endl;

    map<string,int> nameToVertId;
    map<string,cfi::Point*>::iterator vit;
    int i = 1; //from one to be consistent with oce
    for(vit = mapOfVerts.begin(); vit != mapOfVerts.end(); vit++, i++)
    {
        AddVert(i, vit->second);
        nameToVertId[vit->second->getName()] = i;
    }

    cout << "curves " << mapOfEdges.size() << endl;

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
        AddCurve(i, eit->second, ids[0], ids[1]);
    }


    exit(-1);
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
}

}
}
