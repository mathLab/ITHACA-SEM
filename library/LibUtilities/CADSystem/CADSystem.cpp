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

#include <string>
#include <sstream>
#include <limits>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <LibUtilities/CADSystem/CADSystem.h>

using namespace std;

namespace Nektar {
namespace LibUtilities {

string CADSystem::GetName()
{
    return m_name;
}

void CADSystem::Report()
{
    cout << endl << "CAD report:" << endl;
    cout << "\tCAD has: " << m_curves.size() << " curves." << endl;
    cout << "\tCAD has: " << m_surfs.size() << " surfaces." << endl;
    cout << "\tCAD Euler-Poincaré characteristic: " << m_epc << endl;
}

Array<OneD, NekDouble> CADSystem::GetBoundingBox()
{
    Array<OneD, NekDouble> bound(6);
    bound[0] = numeric_limits<double>::max(); //xmin
    bound[1] = numeric_limits<double>::min(); //xmax
    bound[2] = numeric_limits<double>::max(); //ymin
    bound[3] = numeric_limits<double>::min(); //ymax
    bound[4] = numeric_limits<double>::max(); //zmin
    bound[5] = numeric_limits<double>::min(); //zmax

    for(int i = 1; i <= m_curves.size(); i++)
    {
        gp_Pnt start, end;
        CADCurveSharedPtr c = GetCurve(i);
        Array<OneD, NekDouble> ends = c->GetMinMax();

        bound[0] = min(bound[0], min(ends[0],ends[3]));
        bound[1] = max(bound[1], max(ends[0],ends[3]));

        bound[2] = min(bound[2], min(ends[1],ends[4]));
        bound[3] = max(bound[3], max(ends[1],ends[4]));

        bound[4] = min(bound[4], min(ends[2],ends[5]));
        bound[5] = max(bound[5], max(ends[2],ends[5]));
    }

    return bound;
}

bool CADSystem::LoadCAD()
{
    if (!boost::filesystem::exists(m_name.c_str()))
    {
        return false;
    }

    string ext;
    size_t pos = m_name.find(".");
    ext = m_name.substr(pos);

    if (boost::iequals(ext,".STEP") || boost::iequals(ext,".STP"))
    {
        // Takes step file and makes OpenCascade shape
        STEPControl_Reader reader;
        reader = STEPControl_Reader();
        reader.ReadFile(m_name.c_str());
        reader.NbRootsForTransfer();
        reader.TransferRoots();
        shape = reader.OneShape();
        if(shape.IsNull())
        {
            return false;
        }
    }
    else if(boost::iequals(ext,".IGES") || boost::iequals(ext,".IGS"))
    {
        // Takes IGES file and makes OpenCascade shape
        IGESControl_Reader reader;
        reader = IGESControl_Reader();
        reader.ReadFile(m_name.c_str());
        reader.NbRootsForTransfer();
        reader.TransferRoots();
        shape = reader.OneShape();
        if(shape.IsNull())
        {
            return false;
        }
    }
    else
    {
        return false;
    }

    Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
    sfs->Init(shape);
    sfs->Perform();

    if(sfs->Status(ShapeExtend_DONE) )
    {
        shape = sfs->Shape();
    }
    else if(sfs->Status(ShapeExtend_FAIL))
    {
        ASSERTL0(false,"Shape could not be fixed");
    }

    Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe(shape);
    sfwf->ModeDropSmallEdges() = Standard_True;
    sfwf->FixSmallEdges();
    sfwf->FixWireGaps();

    if(sfwf->StatusWireGaps(ShapeExtend_DONE) )
    {
        shape = sfwf->Shape();
    }
    else if(sfwf->StatusWireGaps(ShapeExtend_FAIL))
    {
        ASSERTL0(false,"Shape could not be fixed");
    }
    if(sfwf->StatusSmallEdges(ShapeExtend_DONE) )
    {
        shape = sfwf->Shape();
    }
    else if(sfwf->StatusSmallEdges(ShapeExtend_FAIL))
    {
        ASSERTL0(false,"Shape could not be fixed");
    }

    TopTools_IndexedMapOfShape s;
    TopExp::MapShapes(shape, TopAbs_SHELL, s);

    Handle(ShapeFix_Shell) sfsh = new ShapeFix_Shell;
    sfsh->FixFaceOrientation(TopoDS::Shell(s.FindKey(1)), false, false);

    if(sfsh->Status(ShapeExtend_DONE) )
    {
        shape = sfsh->Shape();
    }
    else if(sfsh->Status(ShapeExtend_FAIL))
    {
        ASSERTL0(false,"Shape could not be fixed");
    }



    // From OpenCascade maps calculate Euler-Poincare number.
    TopTools_IndexedMapOfShape mapOfVerts, ec;
    TopTools_IndexedMapOfShape mapOfFaces;
    TopExp::MapShapes(shape, TopAbs_FACE, mapOfFaces);
    TopExp::MapShapes(shape, TopAbs_EDGE, ec);
    TopExp::MapShapes(shape, TopAbs_VERTEX, mapOfVerts);

    m_epc = mapOfVerts.Extent() - ec.Extent() + mapOfFaces.Extent();

    TopTools_IndexedMapOfShape mapOfEdges; //empty map which is manually built from valid edges

    //standard mm to m conversion
    gp_Trsf transform;
    gp_Pnt ori(0.0, 0.0, 0.0);
    transform.SetScale(ori, 1.0 / 1000.0);
    TopLoc_Location mv(transform);

    for(int i = 1; i <= mapOfVerts.Extent(); i++)
    {
        TopoDS_Shape v = mapOfVerts.FindKey(i);
        v.Move(mv);
        gp_Pnt sp = BRep_Tool::Pnt(TopoDS::Vertex(v));
        Array<OneD, NekDouble> p(3);
        p[0] = sp.X(); p[1] = sp.Y(); p[2] = sp.Z();
        m_verts.push_back(p);
    }

    // For each face of the geometry, get the local edges which bound it. If
    // they are valid (their type != 7), then add them to an edge map. This
    // filters out the dummy edges which OCC uses.
    for(int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face= mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape localEdges;
        TopExp::MapShapes(face, TopAbs_EDGE, localEdges);

        for(int j = 1; j <= localEdges.Extent(); j++)
        {
            TopoDS_Shape edge = localEdges.FindKey(j);
            BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
            if(curve.GetType() != 7)
            {
                if(!(mapOfEdges.Contains(edge)))
                {
                    mapOfEdges.Add(edge);
                }
            }
        }
    }

    map<int, vector<int> > adjsurfmap;

    // Adds edges to our type and map
    for(int i = 1; i <= mapOfEdges.Extent(); i++)
    {
        TopoDS_Shape edge = mapOfEdges.FindKey(i);
        TopoDS_Vertex fv = TopExp::FirstVertex(TopoDS::Edge(edge), Standard_True);
        TopoDS_Vertex lv = TopExp::LastVertex (TopoDS::Edge(edge), Standard_True);

        if(edge.Orientation() == 0)
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(fv), mapOfVerts.FindIndex(lv));
        }
        else
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(lv), mapOfVerts.FindIndex(fv));
        }
    }

    // For each face, examine all the wires (i.e. bounding loops) and
    // investigates the loop. Using this information, connectivity is determined
    // and edges are associated with surfaces.
    for(int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        vector<vector<pair<int,int> > > edges;

        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape mapOfWires;
        TopExp::MapShapes(face, TopAbs_WIRE, mapOfWires);

        for(int j = 1; j <= mapOfWires.Extent(); j++)
        {
            vector<pair<int,int> > edgeloop;

            TopoDS_Shape wire = mapOfWires.FindKey(j);

            ShapeAnalysis_Wire wiretest(TopoDS::Wire(wire),
                                        TopoDS::Face(face),
                                        1E-6);

            BRepTools_WireExplorer exp;

            exp.Init(TopoDS::Wire(wire));

            while(exp.More())
            {
                TopoDS_Shape edge = exp.Current();

                if(mapOfEdges.Contains(edge))
                {
                    pair<int,int> e;
                    e.first = mapOfEdges.FindIndex(edge);
                    adjsurfmap[e.first].push_back(i);
                    e.second = exp.Orientation();
                    edgeloop.push_back(e);
                }

                exp.Next();
            }

            edges.push_back(edgeloop);
        }

        AddSurf(i, face, edges);
    }

    // This checks that all edges are bound by two surfaces, sanity check.
    for(map<int,vector<int> >::iterator it = adjsurfmap.begin();
        it != adjsurfmap.end(); it++)
    {
        ASSERTL0(it->second.size() == 2, "no three curve surfaces");
        m_curves[it->first]->SetAdjSurf(it->second);
    }

    return true;
}

void CADSystem::AddCurve(int i, TopoDS_Shape in, int fv, int lv)
{
    CADCurveSharedPtr newCurve = MemoryManager<CADCurve>::
                                            AllocateSharedPtr(i,in);

    vector<int> vs;
    vs.push_back(fv-1);
    vs.push_back(lv-1);
    m_curves[i] = newCurve;
    m_curves[i]->SetVert(vs);
}

void CADSystem::AddSurf(int i, TopoDS_Shape in,
                        std::vector<std::vector<std::pair<int,int> > > ein)
{
    CADSurfSharedPtr newSurf = MemoryManager<CADSurf>::
                                            AllocateSharedPtr(i,in,ein);
    m_surfs[i] = newSurf;

    //check the face normal is pointing interior
    Array<OneD, NekDouble> bounds = m_surfs[i]->GetBounds();
    Array<OneD, NekDouble> uv(2);
    uv[0] = bounds[1];
    uv[1] = bounds[2];
    Array<OneD, NekDouble> N = m_surfs[i]->N(uv);
    Array<OneD, NekDouble> P = m_surfs[i]->P(uv);

    //create a test point whihc is one unit normal from the center of the surface
    Array<OneD, NekDouble> loc(3);
    loc[0] = P[0] + 1E-3*N[0];
    loc[1] = P[1] + 1E-3*N[1];
    loc[2] = P[2] + 1E-3*N[2];

    ASSERTL0((in.Orientation()==1 && InsideShape(loc)) ||
             (in.Orientation()==0 && !InsideShape(loc)),
             "cannot determine normal");

    if(in.Orientation()==0)
    {
        m_surfs[i]->SetReverseNomral();
    }
}

bool CADSystem::InsideShape(Array<OneD, NekDouble> loc)
{
    gp_Pnt p(loc[0]*1000.0, loc[1]*1000.0, loc[2]*1000.0);

    BRepClass3d_SolidClassifier test(shape, p, 1E-7);
    if(test.State() == TopAbs_IN || test.State() == TopAbs_ON)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}
}
