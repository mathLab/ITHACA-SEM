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

#include <LibUtilities/BasicUtils/ParseUtils.h>

#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/CADSystem/OCE/CADCurveOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADSurfOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADSystemOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADVertOCE.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <ElCLib.hxx>
#include <gce_MakeCirc.hxx>
#include <gce_MakePln.hxx>

#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <ShapeFix_Solid.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <StlAPI_Writer.hxx>
#include <TopoDS_Solid.hxx>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSystemOCE::key = GetEngineFactory().RegisterCreatorFunction(
    "oce", CADSystemOCE::create, "Uses OCE as cad engine");

void filterModShape(TopTools_DataMapOfShapeShape &modShape, TopoDS_Shape &S)
{
    bool repeat = true;
    while (repeat)
    {
        repeat = false;
        if (modShape.IsBound(S))
        {
            repeat = true;
            S      = modShape.Find(S);
        }
    }
}

bool CADSystemOCE::LoadCAD()
{
    Handle(XSControl_WorkSession) WS;
    Handle(Interface_InterfaceModel) Model;
    Handle(XSControl_TransferReader) TR;
    Handle(Transfer_TransientProcess) TP;
    Handle(XCAFDoc_ShapeTool) STool;

    bool fromStep = false;

    if (m_naca.size() == 0)
    {
        // not a naca profile behave normally
        // could be a geo
        string ext = boost::filesystem::extension(m_name);

        if (boost::iequals(ext, ".geo"))
        {
            shape = BuildGeo(m_name);
        }
        else
        {
            // Takes step file and makes OpenCascade shape
            STEPCAFControl_Reader readerCAF;
            readerCAF.SetNameMode(true);
            readerCAF.SetLayerMode(true);
            readerCAF.SetColorMode(true);

            Handle(TDocStd_Document) document =
                new TDocStd_Document(Storage::Version());
            readerCAF.ReadFile(m_name.c_str());
            readerCAF.Transfer(document);

            STEPControl_Reader reader = readerCAF.Reader();

            WS    = reader.WS();
            Model = WS->Model();
            TR    = WS->TransferReader();
            TP    = TR->TransientProcess();
            XCAFDoc_DocumentTool::ShapeTool(document->Main());

            reader.NbRootsForTransfer();
            reader.TransferRoots();
            shape = reader.OneShape();
            if (shape.IsNull())
            {
                return false;
            }
            fromStep = true;
        }
    }
    else
    {
        shape = BuildNACA(m_name);
    }

    TopExp_Explorer explr;

    TopTools_DataMapOfShapeShape modShape;

    if(!m_2d)
    {
        BRepBuilderAPI_Sewing sew(1e-1);

        for (explr.Init(shape, TopAbs_FACE); explr.More(); explr.Next())
        {
            sew.Add(explr.Current());
        }

        sew.Perform();

        for (explr.Init(shape, TopAbs_FACE); explr.More(); explr.Next())
        {
            if (sew.IsModified(explr.Current()))
            {
                modShape.Bind(explr.Current(), sew.Modified(explr.Current()));
            }
        }

        shape = sew.SewedShape();

        int shell = 0;
        for (explr.Init(shape, TopAbs_SHELL); explr.More(); explr.Next())
        {
            shell++;
        }

        /*ASSERTL0(shell == 1,
          "Was not able to form a topological water tight shell");*/
    }

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

        if (!BRep_Tool::Degenerated(TopoDS::Edge(e)))
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

    // attemps to extract patch names from STEP file
    if(fromStep)
    {
        int nb = Model->NbEntities();
        for (int i = 1; i <= nb; i++)
        {
            if (!Model->Value(i)->DynamicType()->SubType(
                    "StepRepr_RepresentationItem"))
                continue;

            Handle(StepRepr_RepresentationItem) enti =
                Handle(StepRepr_RepresentationItem)::DownCast(Model->Value(i));
            Handle(TCollection_HAsciiString) name = enti->Name();

            if (name->IsEmpty())
                continue;

            Handle(Transfer_Binder) binder = TP->Find(Model->Value(i));
            if (binder.IsNull() || !binder->HasResult())
                continue;

            TopoDS_Shape S = TransferBRep::ShapeResult(TP, binder);

            if (S.IsNull())
                continue;

            if (S.ShapeType() == TopAbs_FACE)
            {
                string s(name->ToCString());

                if (mapOfFaces.Contains(S))
                {
                    int id = mapOfFaces.FindIndex(S);

                    m_surfs[id]->SetName(s);
                }
                else
                {
                    filterModShape(modShape, S);
                    if (mapOfFaces.Contains(S))
                    {
                        int id = mapOfFaces.FindIndex(S);
                        m_surfs[id]->SetName(s);
                    }
                    else
                    {
                        ASSERTL0(false, "Name error");
                    }
                }
            }
        }
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
        for (auto &i : m_curves)
        {
            ASSERTL0(i.second->GetAdjSurf().size() == 2,
                     "topolgy error found, surface not closed");
        }
    }

    if(m_verbose)
    {
        Report();
    }

    return true;
}

void CADSystemOCE::AddVert(int i, TopoDS_Shape in)
{
    CADVertSharedPtr newVert = GetCADVertFactory().CreateInstance(key);

    std::static_pointer_cast<CADVertOCE>(newVert)->Initialise(i, in);

    m_verts[i] = newVert;
}

void CADSystemOCE::AddCurve(int i, TopoDS_Shape in)
{
    CADCurveSharedPtr newCurve = GetCADCurveFactory().CreateInstance(key);
    std::static_pointer_cast<CADCurveOCE>(newCurve)->Initialise(i, in);

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
    std::static_pointer_cast<CADSurfOCE>(newSurf)->Initialise(i, in);

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
    ParseUtils::GenerateVector(m_naca, data);
    ASSERTL0(data.size() == 5, "not a vaild domain");

    int n       = std::stoi(naca);
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

TopoDS_Shape CADSystemOCE::BuildGeo(string geo)
{
    ifstream f;
    f.open(geo.c_str());

    map<int, string> points;
    map<int, string> lines;
    map<int, string> splines;
    map<int, string> bsplines;
    map<int, string> circles;
    map<int, string> ellipses;
    map<int, string> lineloops;
    map<int, string> planeSurfs;
    map<int, string> ruledSurfs;
    map<int, string> surfLoops;
    map<int, string> volumes;

    string fline;
    string flinetmp;

    while (getline(f, fline))
    {
        boost::erase_all(fline, "\r");

        vector<string> tmp1, tmp2;
        boost::split(tmp1, fline, boost::is_any_of("//"));
        fline = tmp1[0];

        if (!boost::contains(fline, ";"))
        {
            flinetmp += fline;
            continue;
        }

        fline = flinetmp + fline;
        flinetmp.clear();

        boost::split(tmp1, fline, boost::is_any_of("="));

        boost::split(tmp2, tmp1[0], boost::is_any_of("("));

        string type = tmp2[0];
        boost::erase_all(tmp2[1], ")");
        boost::erase_all(tmp2[1], " ");
        int id = std::stoi(tmp2[1]);

        boost::erase_all(tmp1[1], " ");
        boost::erase_all(tmp1[1], "{");
        boost::erase_all(tmp1[1], "}");
        boost::erase_all(tmp1[1], ";");

        string var = tmp1[1];

        if (boost::iequals(type, "Point"))
        {
            ASSERTL0(points.find(id) == points.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            points[id] = var;
        }
        else if (boost::iequals(type, "Line"))
        {
            ASSERTL0(lines.find(id) == lines.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            lines[id] = var;
        }
        else if (boost::iequals(type, "Spline"))
        {
            ASSERTL0(splines.find(id) == splines.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            splines[id] = var;
        }
        else if (boost::iequals(type, "BSpline"))
        {
            ASSERTL0(bsplines.find(id) == bsplines.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            bsplines[id] = var;
        }
        else if (boost::iequals(type, "Circle"))
        {
            ASSERTL0(circles.find(id) == circles.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            circles[id] = var;
        }
        else if (boost::iequals(type, "Ellipse"))
        {
            ASSERTL0(ellipses.find(id) == ellipses.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            ellipses[id] = var;
        }
        else if (boost::iequals(type, "Line Loop"))
        {
            ASSERTL0(lineloops.find(id) == lineloops.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            boost::erase_all(var, "-");
            lineloops[id] = var;
        }
        else if (boost::iequals(type, "Plane Surface"))
        {
            ASSERTL0(planeSurfs.find(id) == planeSurfs.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            planeSurfs[id] = var;
        }
        else if (boost::iequals(type, "Ruled Surface"))
        {
            ASSERTL0(ruledSurfs.find(id) == ruledSurfs.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            ruledSurfs[id] = var;
        }
        else if (boost::iequals(type, "Surface Loop"))
        {
            ASSERTL0(surfLoops.find(id) == surfLoops.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            surfLoops[id] = var;
        }
        else if (boost::iequals(type, "Volume"))
        {
            ASSERTL0(volumes.find(id) == volumes.end(),
             "Duplicate "+ type + "(" 
             + std::to_string(id) + ")");
            boost::erase_all(var, "-");
            volumes[id] = var;
        }
        else
        {
            cout << "Geometry Command Unknown: " << type << endl;
        }
    }

    map<int, string>::iterator it;

    // build points
    map<int, gp_Pnt> cPoints;
    for (it = points.begin(); it != points.end(); it++)
    {
        vector<NekDouble> data;
        ParseUtils::GenerateVector(it->second, data);

        cPoints[it->first] =
            gp_Pnt(data[0] * 1000.0, data[1] * 1000.0, data[2] * 1000.0);
    }

    // build edges
    map<int, TopoDS_Edge> cEdges;
    for (it = lines.begin(); it != lines.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);
        BRepBuilderAPI_MakeEdge em(cPoints[data[0]], cPoints[data[1]]);
        cEdges[it->first] = em.Edge();
    }
    for (it = splines.begin(); it != splines.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);

        TColgp_Array1OfPnt pointArray(0, data.size() - 1);

        for (int i = 0; i < data.size(); i++)
        {
            pointArray.SetValue(i, cPoints[data[i]]);
        }
        GeomAPI_PointsToBSpline spline(pointArray);
        Handle(Geom_BSplineCurve) curve = spline.Curve();

        BRepBuilderAPI_MakeEdge em(curve);
        cEdges[it->first] = em.Edge();
    }
    for (it = bsplines.begin(); it != bsplines.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);

        TColgp_Array1OfPnt pointArray(0, data.size() - 1);

        for (int i = 0; i < data.size(); i++)
        {
            pointArray.SetValue(i, cPoints[data[i]]);
        }
        Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(pointArray);

        BRepBuilderAPI_MakeEdge em(curve);
        cEdges[it->first] = em.Edge();
    }
    for (it = circles.begin(); it != circles.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);

        ASSERTL0(data.size() == 3, "Wrong definition of circle arc");
        gp_Pnt start  = cPoints[data[0]];
        gp_Pnt centre = cPoints[data[1]];
        gp_Pnt end    = cPoints[data[2]];

        Standard_Real radius = start.Distance(centre);
        gce_MakeCirc mkArc(
            centre, gce_MakePln(start, centre, end).Value(), radius);

        const gp_Circ &circ = mkArc.Value();
        Handle(Geom_Circle) c = new Geom_Circle(circ);

        Standard_Real alpha1 = ElCLib::Parameter(circ, start);
        Standard_Real alpha2 = ElCLib::Parameter(circ, end);
        if ((alpha1 > alpha2) ^ (fabs(alpha2 - alpha1) > M_PI))
        {
            std::swap(alpha1, alpha2);
        }

        Handle(Geom_TrimmedCurve) tc = new Geom_TrimmedCurve(
            c, alpha1, alpha2, false);

        BRepBuilderAPI_MakeEdge em(tc);
        em.Build();
        cEdges[it->first] = em.Edge();
    }
    for (it = ellipses.begin(); it != ellipses.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);

        ASSERTL0(data.size() == 4, "Wrong definition of ellipse arc");
        gp_Pnt start  = cPoints[data[0]];
        gp_Pnt centre = cPoints[data[1]];
        // data[2] useless??
        gp_Pnt end = cPoints[data[3]];

        NekDouble major = start.Distance(centre);

        gp_Vec v1(centre, start);
        gp_Vec v2(centre, end);

        gp_Vec vx(1.0, 0.0, 0.0);
        NekDouble angle = v1.Angle(vx);
        // Check for negative rotation
        if (v1.Y() < 0)
        {
            angle *= -1;
        }

        v2.Rotate(gp_Ax1(), angle);
        NekDouble minor = fabs(v2.Y() / sqrt(1.0 - v2.X() * v2.X() / (major * major)));

        gp_Elips e;
        e.SetLocation(centre);
        e.SetMajorRadius(major);
        e.SetMinorRadius(minor);
        e.Rotate(e.Axis(), angle);
        Handle(Geom_Ellipse) ge = new Geom_Ellipse(e);

        ShapeAnalysis_Curve sac;
        NekDouble p1, p2;
        sac.Project(ge, start, 1e-8, start, p1);
        sac.Project(ge, end, 1e-8, end, p2);

        // Make sure the arc is always of length less than pi
        if (fabs(p2 - p1) > M_PI)
        {
            swap(p1, p2);
        }

        Handle(Geom_TrimmedCurve) tc = new Geom_TrimmedCurve(ge, p1, p2, true);

        BRepBuilderAPI_MakeEdge em(tc);
        em.Build();
        cEdges[it->first] = em.Edge();
    }

    // build wires
    map<int, TopoDS_Wire> cWires;
    for (it = lineloops.begin(); it != lineloops.end(); it++)
    {
        vector<unsigned int> data;
        ParseUtils::GenerateVector(it->second, data);
        BRepBuilderAPI_MakeWire wm;
        for (int i = 0; i < data.size(); i++)
        {
            wm.Add(cEdges[data[i]]);
        }
        cWires[it->first] = wm.Wire();
    }
    //Check if 2D
    if (ruledSurfs.size() == 0)
    {
        if (planeSurfs.size() == 1)
        {
            //Build Plane Surface
            it = planeSurfs.begin();
            vector<unsigned int> data;
            ParseUtils::GenerateVector(it->second, data);
            BRepBuilderAPI_MakeFace face(cWires[data[0]], true);
            for (int i = 1; i < data.size(); i++)
            {
                face.Add(cWires[data[i]]);
            }

            ASSERTL0(face.Error() == BRepBuilderAPI_FaceDone, "build geo failed");

            ShapeFix_Face sf(face.Face());
            sf.FixOrientation();

            return sf.Face();    
        }
                
    }
    //Run 3D
    else
    {
        // Build Surfaces
        map<int, TopoDS_Face> cFaces;
        // Build Place Surfaces
        for (auto &surf : planeSurfs)
        {
            vector<unsigned int> data;
            ParseUtils::GenerateVector(surf.second, data);
            ASSERTL0(data.size() == 1,
                    "Ruled surface should have only one curve loop");

            BRepBuilderAPI_MakeFace faceBuilder(cWires[data[0]]);

            cFaces[surf.first] = faceBuilder.Face();
        }
        // Build Ruled Surfaces
        for (auto &surf : ruledSurfs)
        {
            vector<unsigned int> data;
            ParseUtils::GenerateVector(surf.second, data);
            ASSERTL0(data.size() == 1,
                    "Ruled surface should have only one curve loop");

            BRepFill_Filling fill(2, 30);
            TopExp_Explorer ex;

            for (ex.Init(cWires[data[0]], TopAbs_EDGE); ex.More(); ex.Next())
            {
                const TopoDS_Edge &shape = TopoDS::Edge(ex.Current());
                fill.Add(shape, GeomAbs_C0);
            }

            fill.Build();
            TopoDS_Face fixedFace = fill.Face();
            cFaces[surf.first] = fixedFace;

            /*
            TopExp_Explorer ex;
            std::vector<Handle(GeomFill_SimpleBound)> bounds;

            for (ex.Init(cWires[data[0]], TopAbs_EDGE); ex.More(); ex.Next())
            {
                double s0, s1;
                Handle(Geom_Curve) c = BRep_Tool::Curve(
                    TopoDS::Edge(ex.Current()), s0, s1);
                Handle(GeomAdaptor_HCurve) Curve = new GeomAdaptor_HCurve(c, s0, s1);
                bounds.push_back(new GeomFill_SimpleBound(Curve, 0.001, 0.001));
            }

            GeomFill_ConstrainedFilling fill(2, 8);
            fill.Init(bounds[0], bounds[1], bounds[2]);

            Handle(Geom_BSplineSurface) fillSurf = fill.Surface();
            BRepBuilderAPI_MakeFace mf(fillSurf, cWires[data[0]]);
            mf.Build();
            cFaces[surf.first] = mf.Face();

            StlAPI_Writer writer;
            BRepMesh_IncrementalMesh Mesh(cFaces[surf.first], 10.0);
            Mesh.Perform();
            std::string outfname = "out-" + boost::lexical_cast<std::string>(
                surf.first) + ".stl";
            writer.Write(cFaces[surf.first], outfname.c_str());
            */
        }
        // Build Shells
        map<int, TopoDS_Shell> cShells;
        for (auto &sloop : surfLoops)
        {
            vector<unsigned int> data;
            ParseUtils::GenerateVector(sloop.second, data);

            BRepBuilderAPI_Sewing shellMaker;

            for (int i = 0; i < data.size(); ++i)
            {
                shellMaker.Add(cFaces[data[i]]);
            }

            shellMaker.Perform();

            cShells[sloop.first] = TopoDS::Shell(shellMaker.SewedShape());

            // BRepMesh_IncrementalMesh Mesh( shell, 200.0 );
            // Mesh.Perform();
            // StlAPI_Writer asd;
            // std::string outfname = "out.stl";
            // asd.Write(shell, outfname.c_str());
        }
        // Build Volumes
        map<int, TopoDS_Shape> cVolumes;
        for (auto &vol : volumes)
        {
            vector<unsigned int> data;
            ParseUtils::GenerateVector(vol.second, data);

            BRepBuilderAPI_MakeSolid solidMaker;

            for (int i = 0; i < data.size(); ++i)
            {
                solidMaker.Add(cShells[data[i]]);

                if (i == 0)
                {
                    continue;
                }

                // For each shell that's being removed from the solid, add centroid
                // for tetrahedralisation purposes. In general this won't work so
                // well because we could have a non-convex shell, but let's assume
                // it is a decent guess for now.
                GProp_GProps props;
                BRepGProp::VolumeProperties(cShells[data[i]], props);
                gp_Pnt centroid = props.CentreOfMass();

                Array<OneD, NekDouble> voidPt(3);
                voidPt[0] = centroid.X();
                voidPt[1] = centroid.Y();
                voidPt[2] = centroid.Z();
                m_voidPoints.push_back(voidPt);
            }

            TopoDS_Solid s = solidMaker.Solid();

            // Perform fix on solid
            ShapeFix_Solid solidFix(s);

            solidFix.Perform();

            return solidFix.Solid();
        }
    }
}
}
}
