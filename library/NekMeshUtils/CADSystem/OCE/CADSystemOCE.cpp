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
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/CADSystem/OCE/CADCurveOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADSurfOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADSystemOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADVertOCE.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSystemOCE::key = GetEngineFactory().RegisterCreatorFunction(
    "oce", CADSystemOCE::create, "Uses OCE as cad engine");

bool CADSystemOCE::LoadCAD()
{
    if (m_naca.size() == 0)
    {
        // not a naca profile behave normally
        // Takes step file and makes OpenCascade shape
        STEPControl_Reader reader;
        reader = STEPControl_Reader();
        reader.ReadFile(m_name.c_str());
        reader.NbRootsForTransfer();
        reader.TransferRoots();
        shape = reader.OneShape();
        if (shape.IsNull())
        {
            return false;
        }
    }
    else
    {
        shape = BuildNACA(m_name);
    }

    TopExp_Explorer explr;

    // build map of verticies
    for (explr.Init(shape, TopAbs_VERTEX); explr.More(); explr.Next())
    {
        TopoDS_Shape v = explr.Current();
        if (mapOfVerts.Contains(v))
        {
            continue;
        }
        int i = mapOfVerts.Add(v);
        AddVert(i, v);
    }

    // For each face of the geometry, get the local edges which bound it. If
    // they are valid (their type != 7), then add them to an edge map. This
    // filters out the dummy edges which OCC uses.
    for (explr.Init(shape, TopAbs_EDGE); explr.More(); explr.Next())
    {
        TopoDS_Shape e = explr.Current().Oriented(TopAbs_FORWARD);
        if (mapOfEdges.Contains(e))
        {
            continue;
        }
        BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(e));
        if (curve.GetType() != 7)
        {
            int i = mapOfEdges.Add(e);
            AddCurve(i, e);
        }
    }

    for (explr.Init(shape, TopAbs_FACE); explr.More(); explr.Next())
    {
        TopoDS_Shape f = explr.Current();
        ASSERTL0(!mapOfFaces.Contains(f), "duplicated faces");
        int i = mapOfFaces.Add(f);

        AddSurf(i, f);
    }

    // attempts to identify properties of the vertex on the degen edge
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i).Oriented(TopAbs_FORWARD);

        TopTools_IndexedMapOfShape localEdges;
        TopExp::MapShapes(face, TopAbs_EDGE, localEdges);

        for (int j = 1; j <= localEdges.Extent(); j++)
        {
            TopoDS_Shape edge = localEdges.FindKey(j);
            if (BRep_Tool::Degenerated(TopoDS::Edge(edge)))
            {
                gp_Pnt2d p1, p2;

                BRep_Tool::UVPoints(TopoDS::Edge(edge), TopoDS::Face(face), p1,
                                    p2);

                m_verts[mapOfVerts.FindIndex(TopExp::FirstVertex(
                            TopoDS::Edge(edge), Standard_True))]
                    ->SetDegen(i, m_surfs[i], (p1.X() + p2.X()) / 2.0,
                               (p1.Y() + p2.Y()) / 2.0);
            }
        }
    }

    // This checks that all edges are bound by two surfaces, sanity check.
    if (!m_2d)
    {
        map<int, CADCurveSharedPtr>::iterator it;
        for (it = m_curves.begin(); it != m_curves.end(); it++)
        {
            ASSERTL0(it->second->GetAdjSurf().size() == 2,
                     "curve is not joined to 2 surfaces");
        }
    }

    return true;
}

void CADSystemOCE::AddVert(int i, TopoDS_Shape in)
{
    CADVertSharedPtr newVert = GetCADVertFactory().CreateInstance(key);

    boost::static_pointer_cast<CADVertOCE>(newVert)->Initialise(i, in);

    m_verts[i] = newVert;
}

void CADSystemOCE::AddCurve(int i, TopoDS_Shape in)
{
    CADCurveSharedPtr newCurve = GetCADCurveFactory().CreateInstance(key);
    boost::static_pointer_cast<CADCurveOCE>(newCurve)->Initialise(i, in);

    TopoDS_Vertex fv = TopExp::FirstVertex(TopoDS::Edge(in));
    TopoDS_Vertex lv = TopExp::LastVertex(TopoDS::Edge(in));

    vector<CADVertSharedPtr> vs;
    vs.push_back(m_verts[mapOfVerts.FindIndex(fv)]);
    vs.push_back(m_verts[mapOfVerts.FindIndex(lv)]);
    newCurve->SetVert(vs);

    m_curves[i] = newCurve;
}

void CADSystemOCE::AddSurf(int i, TopoDS_Shape in)
{
    CADSurfSharedPtr newSurf = GetCADSurfFactory().CreateInstance(key);
    boost::static_pointer_cast<CADSurfOCE>(newSurf)->Initialise(i, in);

    // do the exploration on forward oriented
    TopoDS_Shape face = in.Oriented(TopAbs_FORWARD);
    TopTools_IndexedMapOfShape mapOfWires;
    TopExp::MapShapes(face, TopAbs_WIRE, mapOfWires);
    vector<EdgeLoopSharedPtr> edgeloops;
    // now we acutally analyse the loops for cad building
    for (int j = 1; j <= mapOfWires.Extent(); j++)
    {
        EdgeLoopSharedPtr edgeloop = EdgeLoopSharedPtr(new EdgeLoop);

        TopoDS_Shape wire = mapOfWires.FindKey(j);

        BRepTools_WireExplorer exp;

        exp.Init(TopoDS::Wire(wire));

        while (exp.More())
        {
            TopoDS_Shape edge = exp.Current();

            if (mapOfEdges.Contains(edge))
            {
                int e = mapOfEdges.FindIndex(edge);
                edgeloop->edges.push_back(m_curves[e]);
                edgeloop->edgeo.push_back(exp.Orientation() == TopAbs_FORWARD
                                              ? CADOrientation::eForwards
                                              : CADOrientation::eBackwards);
            }
            exp.Next();
        }
        edgeloops.push_back(edgeloop);
    }

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
}

Array<OneD, NekDouble> CADSystemOCE::GetBoundingBox()
{
    Array<OneD, NekDouble> bound(6);
    bound[0] = numeric_limits<double>::max(); // xmin
    bound[1] = numeric_limits<double>::min(); // xmax
    bound[2] = numeric_limits<double>::max(); // ymin
    bound[3] = numeric_limits<double>::min(); // ymax
    bound[4] = numeric_limits<double>::max(); // zmin
    bound[5] = numeric_limits<double>::min(); // zmax

    for (int i = 1; i <= m_curves.size(); i++)
    {
        CADCurveSharedPtr c = GetCurve(i);
        Array<OneD, NekDouble> ends = c->GetMinMax();

        bound[0] = min(bound[0], min(ends[0], ends[3]));
        bound[1] = max(bound[1], max(ends[0], ends[3]));

        bound[2] = min(bound[2], min(ends[1], ends[4]));
        bound[3] = max(bound[3], max(ends[1], ends[4]));

        bound[4] = min(bound[4], min(ends[2], ends[5]));
        bound[5] = max(bound[5], max(ends[2], ends[5]));
    }

    return bound;
}

TopoDS_Shape CADSystemOCE::BuildNACA(string naca)
{
    ASSERTL0(naca.length() == 4, "not a 4 digit code");
    vector<NekDouble> data;
    ParseUtils::GenerateUnOrderedVector(m_naca.c_str(), data);
    ASSERTL0(data.size() == 5, "not a vaild domain");

    int n       = boost::lexical_cast<int>(naca);
    NekDouble T = (n % 100) / 100.0;
    n /= 100;
    NekDouble P = (n % 10) / 10.0;
    n /= 10;
    NekDouble M = (n % 10) / 100.0;

    int np = 25;

    Array<OneD, NekDouble> xc(np);
    NekDouble dtheta = M_PI / (np - 1);
    for (int i = 0; i < np; i++)
    {
        xc[i] = (1.0 - cos(i * dtheta)) / 2.0;
    }

    Array<OneD, NekDouble> yc(np), dyc(np);
    for (int i = 0; i < np; i++)
    {
        if (xc[i] < P)
        {
            yc[i]  = M / P / P * (2.0 * P * xc[i] - xc[i] * xc[i]);
            dyc[i] = 2.0 * M / P / P * (P - xc[i]);
        }
        else
        {
            yc[i] = M / (1.0 - P) / (1.0 - P) *
                    (1.0 - 2.0 * P + 2.0 * P * xc[i] - xc[i] * xc[i]);
            dyc[i] = 2.0 * M / (1.0 - P) / (1.0 - P) * (P - xc[i]);
        }
    }

    Array<OneD, NekDouble> yt(np);
    for (int i = 0; i < np; i++)
    {
        yt[i] =
            T / 0.2 * (0.2969 * sqrt(xc[i]) - 0.1260 * xc[i] -
                       0.3516 * xc[i] * xc[i] + 0.2843 * xc[i] * xc[i] * xc[i] -
                       0.1015 * xc[i] * xc[i] * xc[i] * xc[i]);
    }

    Array<OneD, NekDouble> x(2 * np - 1), y(2 * np - 1);
    int l = 0;
    for (int i = np - 1; i >= 0; i--, l++)
    {
        NekDouble theta = atan(dyc[i]);

        x[l] = xc[i] - yt[i] * sin(theta);
        y[l] = yc[i] + yt[i] * cos(theta);
    }
    for (int i = 1; i < np; i++)
    {
        NekDouble theta = atan(dyc[i]);

        x[i + np - 1] = xc[i] + yt[i] * sin(theta);
        y[i + np - 1] = yc[i] - yt[i] * cos(theta);
    }

    TColgp_Array1OfPnt pointArray(0, 2 * np - 2);

    for (int i = 0; i < 2 * np - 1; i++)
    {
        pointArray.SetValue(i, gp_Pnt(x[i] * 1000.0, y[i] * 1000.0, 0.0));
    }

    GeomAPI_PointsToBSpline spline(pointArray);
    Handle(Geom_BSplineCurve) curve = spline.Curve();

    BRepBuilderAPI_MakeEdge areoEdgeBuilder(curve);
    TopoDS_Edge aeroEdge = areoEdgeBuilder.Edge();
    BRepBuilderAPI_MakeEdge aeroTEBuilder(
        gp_Pnt(x[0] * 1000.0, y[0] * 1000.0, 0.0),
        gp_Pnt(x[2 * np - 2] * 1000.0, y[2 * np - 2] * 1000.0, 0.0));
    TopoDS_Edge TeEdge = aeroTEBuilder.Edge();

    BRepBuilderAPI_MakeWire aeroWireBuilder(aeroEdge, TeEdge);
    TopoDS_Wire aeroWire = aeroWireBuilder.Wire();

    gp_Trsf transform;
    gp_Ax1 rotAx(gp_Pnt(500.0, 0.0, 0.0), gp_Dir(gp_Vec(0.0, 0.0, -1.0)));
    transform.SetRotation(rotAx, data[4] / 180.0 * M_PI);
    TopLoc_Location mv(transform);
    aeroWire.Move(mv);

    BRepBuilderAPI_MakeEdge domInlBuilder(
        gp_Pnt(data[0] * 1000.0, data[1] * 1000.0, 0.0),
        gp_Pnt(data[0] * 1000.0, data[3] * 1000.0, 0.0));
    TopoDS_Edge inlEdge = domInlBuilder.Edge();

    BRepBuilderAPI_MakeEdge domTopBuilder(
        gp_Pnt(data[0] * 1000.0, data[3] * 1000.0, 0.0),
        gp_Pnt(data[2] * 1000.0, data[3] * 1000.0, 0.0));
    TopoDS_Edge topEdge = domTopBuilder.Edge();

    BRepBuilderAPI_MakeEdge domOutBuilder(
        gp_Pnt(data[2] * 1000.0, data[3] * 1000.0, 0.0),
        gp_Pnt(data[2] * 1000.0, data[1] * 1000.0, 0.0));
    TopoDS_Edge outEdge = domOutBuilder.Edge();

    BRepBuilderAPI_MakeEdge domBotBuilder(
        gp_Pnt(data[2] * 1000.0, data[1] * 1000.0, 0.0),
        gp_Pnt(data[0] * 1000.0, data[1] * 1000.0, 0.0));
    TopoDS_Edge botEdge = domBotBuilder.Edge();

    BRepBuilderAPI_MakeWire domWireBuilder(inlEdge, topEdge, outEdge, botEdge);
    TopoDS_Wire domWire = domWireBuilder.Wire();

    BRepBuilderAPI_MakeFace face(domWire, true);
    face.Add(aeroWire);

    ShapeFix_Face sf(face.Face());
    sf.FixOrientation();


    return sf.Face();
}
}
}
