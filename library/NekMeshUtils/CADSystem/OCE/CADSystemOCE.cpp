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

#include <NekMeshUtils/CADSystem/OCE/GeoParser.hpp>
#include <NekMeshUtils/CADSystem/OCE/TransfiniteSurface.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <ElCLib.hxx>
#include <gce_MakeCirc.hxx>
#include <gce_MakePln.hxx>

#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Solid.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <StlAPI_Writer.hxx>
#include <TopoDS_Solid.hxx>
#include <BRepClass3d_SolidClassifier.hxx>

#include <Geom2d_Line.hxx>

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

/**
 * @brief For a given shell @p shell, attempt to find a point that is strictly
 * on the interior of the shape.
 */
gp_Pnt FindInteriorPoint(TopoDS_Shell &shell)
{
    // Grab a face.
    TopExp_Explorer ex;
    ex.Init(shell, TopAbs_FACE);

    BRepBuilderAPI_MakeSolid makeSolid(shell);
    TopoDS_Solid solid = makeSolid.Solid();

    // Compute solid volume. if it's negative, then the shell was constructed
    // inside-out, so reverse the solid.
    GProp_GProps solidProps;
    BRepGProp::VolumeProperties(solid, solidProps);
    if (solidProps.Mass() < 0.0)
    {
        solid.Reverse();
    }

    BRepClass3d_SolidClassifier cls(solid);

    for (ex.Init(shell, TopAbs_FACE); ex.More(); ex.Next())
    {
        // Get bounds for the face and grab a handle to a Geom_Surface.
        TopoDS_Face face = TopoDS::Face(ex.Current());
        Handle(Geom_Surface) s = BRep_Tool::Surface(face);

        NekDouble umin, umax, vmin, vmax;
        BRepTools::UVBounds(face, umin, umax, vmin, vmax);

        // Take central point on parametrisation and compute normal, which
        // should be outwards facing assuming the shell has been fixed.
        GeomLProp_SLProps d(s, 2, Precision::Confusion());
        d.SetParameters(0.5 * (umin + umax), 0.5 * (vmin + vmax));

        gp_XYZ pnt = s->Value(0.5 * (umin + umax), 0.5 * (vmin + vmax)).XYZ();

        // If for some reason we don't have a normal for this face, continue to
        // the next one.
        if (!d.IsNormalDefined())
        {
            continue;
        }

        // Compute inwards facing normal and multiply by a scale factor which is
        // hopefully somewhat proportional to the direction we care about.
        GProp_GProps faceProp;
        BRepGProp::SurfaceProperties(face, faceProp);

        // Follow the inwards normal direction and backtrack. Use the classifier
        // to figure out if we're inside the volume (or not).
        gp_XYZ n = d.Normal().XYZ() * sqrt(faceProp.Mass());
        NekDouble alpha = 1.0;

        for (int i = 0; i < 2; ++i)
        {
            // Try inwards and outwards facing directions.
            n *= -1.0;

            while (alpha > 1e-4)
            {
                gp_Pnt testPnt(pnt + n * alpha);
                cls.Perform(testPnt, Precision::Confusion());
                TopAbs_State state = cls.State();

                if (state == TopAbs_IN)
                {
                    return testPnt;
                }

                alpha /= 2.0;
            }
        }

        // Failed on this face, try the next one.
    }

    ASSERTL0(false, "Failed to find point internal to geometry.");
    return gp_Pnt();
}

bool CADSystemOCE::LoadCAD()
{
    Handle(XSControl_WorkSession) WS;
    Handle(Interface_InterfaceModel) Model;
    Handle(XSControl_TransferReader) TR;
    Handle(Transfer_TransientProcess) TP;
    Handle(XCAFDoc_ShapeTool) STool;

    bool fromStep = false;

    if (m_config.count("UseNACA") == 0)
    {
        // not a naca profile behave normally
        // could be a geo
        string ext = boost::filesystem::extension(m_name);

        if (boost::iequals(ext, ".geo"))
        {
            m_shape = BuildGeo(m_name);
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
            m_shape = reader.OneShape();
            if (m_shape.IsNull())
            {
                return false;
            }
            fromStep = true;
        }
    }
    else
    {
        m_shape = BuildNACA(m_name);
    }

    TopExp_Explorer explr;

    TopTools_DataMapOfShapeShape modShape;

    if(!m_2d)
    {
        const NekDouble sewTolerance = 0.1;
        BRepBuilderAPI_Sewing sew(sewTolerance);

        for (explr.Init(m_shape, TopAbs_FACE); explr.More(); explr.Next())
        {
            sew.Add(explr.Current());
        }

        sew.Perform();

        for (explr.Init(m_shape, TopAbs_FACE); explr.More(); explr.Next())
        {
            if (sew.IsModified(explr.Current()))
            {
                modShape.Bind(explr.Current(), sew.Modified(explr.Current()));
            }
        }

        m_shape = sew.SewedShape();
    }

    // build map of verticies
    for (explr.Init(m_shape, TopAbs_VERTEX); explr.More(); explr.Next())
    {
        TopoDS_Shape v = explr.Current();
        if (m_mapOfVerts.Contains(v))
        {
            continue;
        }
        int i = m_mapOfVerts.Add(v);
        AddVert(i, v);
    }

    // For each face of the geometry, get the local edges which bound it. If
    // they are valid (their type != 7), then add them to an edge map. This
    // filters out the dummy edges which OCC uses.
    for (explr.Init(m_shape, TopAbs_EDGE); explr.More(); explr.Next())
    {
        TopoDS_Shape e = explr.Current().Oriented(TopAbs_FORWARD);
        if (m_mapOfEdges.Contains(e))
        {
            continue;
        }

        if (!BRep_Tool::Degenerated(TopoDS::Edge(e)))
        {
            int i = m_mapOfEdges.Add(e);
            AddCurve(i, e);
        }
    }

    for (explr.Init(m_shape, TopAbs_FACE); explr.More(); explr.Next())
    {
        TopoDS_Shape f = explr.Current();
        ASSERTL0(!m_mapOfFaces.Contains(f), "duplicated faces");
        int i = m_mapOfFaces.Add(f);

        AddSurf(i, f);
    }

    // Attempts to extract patch names from STEP file
    if (fromStep)
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

                if (m_mapOfFaces.Contains(S))
                {
                    int id = m_mapOfFaces.FindIndex(S);

                    m_surfs[id]->SetName(s);
                }
                else
                {
                    filterModShape(modShape, S);
                    if (m_mapOfFaces.Contains(S))
                    {
                        int id = m_mapOfFaces.FindIndex(S);
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
    for (int i = 1; i <= m_mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = m_mapOfFaces.FindKey(i).Oriented(TopAbs_FORWARD);

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

                m_verts[m_mapOfVerts.FindIndex(TopExp::FirstVertex(
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
    vs.push_back(m_verts[m_mapOfVerts.FindIndex(fv)]);
    vs.push_back(m_verts[m_mapOfVerts.FindIndex(lv)]);
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

            if (m_mapOfEdges.Contains(edge))
            {
                int e = m_mapOfEdges.FindIndex(edge);
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
    ASSERTL0(naca.length() == 4, "not a 4 digit code: " + naca);
    vector<NekDouble> data;
    ParseUtils::GenerateVector(m_config["UseNACA"], data);
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

/**
 * @brief Helper function for Gmsh file construction.
 */
template<typename T>
inline bool ContainsIDs(std::string                name,
                        unsigned int               id,
                        std::string                contName,
                        std::vector<int>          &facetids,
                        std::map<unsigned int, T> &toSearch)
{
    bool valid = true;
    for (auto &fid : facetids)
    {
        if (toSearch.find(fid) == toSearch.end())
        {
            NEKERROR(ErrorUtil::ewarning,
                     name + " ID " + std::to_string(id) + " refers to "
                     "unknown " + contName + " ID " + std::to_string(fid));
        }
    }

    return valid;
}

/**
 * @brief Helper function for Gmsh file construction.
 */
template<typename T>
inline void CheckWarning(std::string                name,
                         unsigned int               id,
                         std::map<unsigned int, T> &toSearch)
{
    if (toSearch.find(id) != toSearch.end())
    {
        NEKERROR(ErrorUtil::ewarning,
                 "Duplicate .geo " + name + " " + std::to_string(id) +
                 " found, will be overwritten.");
    }
}

/**
 * @brief Create a OpenCASCADE object from a Gmsh file @p geo.
 *
 * This routine implements a reader for simple Gmsh .geo files. Currently, it
 * supports:
 *
 * - 0D: points;
 * - 1D: lines, splines, bsplines, circles, ellipses and line loops;
 * - 2D: surfaces and ruled surfaces;
 * - 3D: volumes
 *
 * Most other information, including physical lines, surfaces and volumes is all
 * ignored at present.
 *
 * @param geo   Name of the .geo file to read.
 * @return      An OpenCASCADE TopoDS_Shape that defines the geometry.
 */
TopoDS_Shape CADSystemOCE::BuildGeo(string geo)
{
    using GeoAst::Geom;

    ifstream f(geo.c_str());
    const std::string str((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());

    // Construct a parser for the geo file. Ensure that we use the correct
    // skipper so that comments are ignored.
    LibUtilities::Interpreter interp;
    GeoParser<std::string::const_iterator> geoParser(interp);
    CommentSkipper<std::string::const_iterator> skip;
    GeoAst::GeoFile geoFile;

    // Run parser.
    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bool r = qi::phrase_parse(iter, end, geoParser, skip, geoFile);

    // If parsing failed, output an error.
    if (!r)
    {
        std::string::const_iterator some = iter+30;
        std::string context(iter, (some>end)?end:some);
        std::cout << "Parsing of geo file failed\n";
        std::cout << "stopped at: \": " << context << "...\"\n";
    }

    // Build points.
    map<unsigned int, gp_Pnt> cPoints;
    for (auto &point : geoFile.points)
    {
        CheckWarning("point", point.id, cPoints);
        cPoints[point.id] =
            gp_Pnt(point.x * 1000.0, point.y * 1000.0, point.z * 1000.0);
    }

    //
    // Build edges.
    //
    map<unsigned int, TopoDS_Edge> cEdges;
    map<unsigned int, std::string> cEdgeTypes;
    map<unsigned int, Geom> edgeMap;

    // Get straight lines.
    for (auto &line : geoFile.lines)
    {
        if (line.ids.size() != 2)
        {
            NEKERROR(ErrorUtil::ewarning, "Line " + std::to_string(line.id) +
                     " contains more than two points, ignoring.");
            continue;
        }

        CheckWarning("line", line.id, cEdges);
        if (!ContainsIDs("Line", line.id, "point", line.ids, cPoints))
        {
            continue;
        }

        BRepBuilderAPI_MakeEdge em(cPoints[line.ids[0]], cPoints[line.ids[1]]);
        cEdges[line.id] = em.Edge();
        edgeMap[line.id] = line;
        cEdgeTypes[line.id] = "line";
    }

    // Get splines.
    for (auto &spline : geoFile.splines)
    {
        if (spline.ids.size() < 2)
        {
            NEKERROR(ErrorUtil::ewarning, "Spline " + std::to_string(spline.id)
                     + " does not contain enough points, ignoring.");
            continue;
        }

        CheckWarning("spline", spline.id, cEdges);
        if (!ContainsIDs("Spline", spline.id, "point", spline.ids, cPoints))
        {
            continue;
        }

        TColgp_Array1OfPnt pointArray(0, spline.ids.size() - 1);

        for (int i = 0; i < spline.ids.size(); i++)
        {
            pointArray.SetValue(i, cPoints[spline.ids[i]]);
        }
        GeomAPI_PointsToBSpline oceSpline(pointArray);
        Handle(Geom_BSplineCurve) curve = oceSpline.Curve();

        BRepBuilderAPI_MakeEdge em(curve);
        cEdges[spline.id] = em.Edge();
        edgeMap[spline.id] = spline;
        cEdgeTypes[spline.id] = "spline";
    }

    // Get B-Splines.
    for (auto &bspline : geoFile.bsplines)
    {
        if (bspline.ids.size() < 2)
        {
            NEKERROR(ErrorUtil::ewarning,
                     "BSpline " + std::to_string(bspline.id) + " does not "
                     "contain enough points, ignoring.");
            continue;
        }

        CheckWarning("bspline", bspline.id, cEdges);
        if (!ContainsIDs("BSpline", bspline.id, "point", bspline.ids, cPoints))
        {
            continue;
        }

        TColgp_Array1OfPnt pointArray(0, bspline.ids.size() - 1);

        for (int i = 0; i < bspline.ids.size(); i++)
        {
            pointArray.SetValue(i, cPoints[bspline.ids[i]]);
        }
        Handle(Geom_BezierCurve) curve = new Geom_BezierCurve(pointArray);

        BRepBuilderAPI_MakeEdge em(curve);
        cEdges[bspline.id] = em.Edge();
        edgeMap[bspline.id] = bspline;
        cEdgeTypes[bspline.id] = "bspline";
    }

    // Get circles.
    for (auto &circle : geoFile.circles)
    {
        if (circle.ids.size() != 3)
        {
            NEKERROR(ErrorUtil::ewarning, "Circle " + std::to_string(circle.id)
                     + " should contain only three points.");
            continue;
        }

        CheckWarning("circle", circle.id, cEdges);
        if (!ContainsIDs("Circle", circle.id, "point", circle.ids, cPoints))
        {
            continue;
        }

        gp_Pnt start  = cPoints[circle.ids[0]];
        gp_Pnt centre = cPoints[circle.ids[1]];
        gp_Pnt end    = cPoints[circle.ids[2]];

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
        cEdges[circle.id] = em.Edge();
        edgeMap[circle.id] = circle;
        cEdgeTypes[circle.id] = "circle";
    }

    // Get ellipses.
    for (auto &ellipse : geoFile.ellipses)
    {
        if (ellipse.ids.size() != 4)
        {
            NEKERROR(ErrorUtil::ewarning,
                     "Ellipse " + std::to_string(ellipse.id) + " should contain"
                     " only four points.");
            continue;
        }

        CheckWarning("ellipse", ellipse.id, cEdges);
        if (!ContainsIDs("Ellipse", ellipse.id, "point", ellipse.ids, cPoints))
        {
            continue;
        }

        gp_Pnt start  = cPoints[ellipse.ids[0]];
        gp_Pnt centre = cPoints[ellipse.ids[1]];
        // data[2] useless??
        gp_Pnt end = cPoints[ellipse.ids[3]];

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
        NekDouble minor = fabs(
            v2.Y() / sqrt(1.0 - v2.X() * v2.X() / (major * major)));

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
            std::swap(p1, p2);
        }

        Handle(Geom_TrimmedCurve) tc = new Geom_TrimmedCurve(ge, p1, p2, true);

        BRepBuilderAPI_MakeEdge em(tc);
        em.Build();
        cEdges[ellipse.id] = em.Edge();
        edgeMap[ellipse.id] = ellipse;
        cEdgeTypes[ellipse.id] = "ellipse";
    }

    //
    // Build wires.
    //
    map<unsigned int, TopoDS_Wire> cWires;
    map<unsigned int, Geom> loopMap;
    for (auto &loop : geoFile.lineLoops)
    {
        // Make IDs all positive since we don't care about the orientations and
        // will get OCC to construct this for us.
        for (auto &id : loop.ids)
        {
            id = std::abs(id);
        }

        CheckWarning("line loop", loop.id, cWires);
        if (!ContainsIDs("Line Loop", loop.id, "edge", loop.ids, cEdges))
        {
            continue;
        }

        BRepBuilderAPI_MakeWire wm;
        for (auto &edgeId : loop.ids)
        {
            wm.Add(cEdges[edgeId]);
        }
        cWires[loop.id] = wm.Wire();
        loopMap[loop.id] = loop;
    }

    //
    // Build faces.
    //
    map<unsigned int, TopoDS_Face> cFaces;
    map<unsigned int, Geom> faceMap;

    // Construct plane surfaces.
    for (auto &planeSurf : geoFile.planeSurfs)
    {
        BRepBuilderAPI_MakeFace face(cWires[planeSurf.ids[0]], true);
        for (int i = 1; i < planeSurf.ids.size(); i++)
        {
            face.Add(cWires[planeSurf.ids[i]]);
        }

        ASSERTL0(face.Error() == BRepBuilderAPI_FaceDone, "build geo failed");

        ShapeFix_Face sf(face.Face());
        sf.FixOrientation();

        cFaces[planeSurf.id] = sf.Face();
        faceMap[planeSurf.id] = planeSurf;
    }

    // At this point, 2D simulations should be good because we have (at most)
    // one surface.
    if (m_2d)
    {
        ASSERTL0(cFaces.size() == 1,
                 "2D simulations should define at most one plane surface in "
                 "the .geo file.");

        return cFaces.begin()->second;
    }

    // Build ruled surfaces. This is probably amongst the most complex parts of
    // this code, since we need to construct a filling for the wire which may be
    // based off a couple of different approaches:
    //
    // 1) If Gmsh is using the OpenCASCADE kernel, then we should always use the
    //    BRepFill_Filling class, which uses some kind of magic to blend
    //    curvature of the edges to construct a filling.
    // 2) If Gmsh is using its built-in kernel, then ruled surfaces are must
    //    have either three of four curves, and then the filling is defined
    //    using a transfinite interpolation.
    // 3) ...unless those three or four curves are circle arcs, in which case we
    //    test to see if we're on a sphere first.
    //
    // Most of this is straightforward, asides from (2) which requires a custom
    // OpenCASCADE class to handle the transfinite interpolation.
    for (auto &surf : geoFile.ruledSurfs)
    {
        CheckWarning("surface", surf.id, cFaces);
        if (!ContainsIDs("Ruled Surface", surf.id, "line loop", surf.ids, cWires))
        {
            continue;
        }
        if (surf.ids.size() != 1)
        {
            NEKERROR(ErrorUtil::ewarning,
                     "Surface " + std::to_string(surf.id) + " should only "
                     "contain a single line loop, ignoring.");
            continue;
        }

        Geom &loop = loopMap[surf.ids[0]];

        // Get edges.
        const int nEdges = loop.ids.size();

        // Special case (3) above: 3 or 4 circular edges should define a
        // spherical surface. This is what the built-in geo kernel supports, but
        // OpenCASCADE engine gives a regular filling.
        bool isSphere = true;
        if (nEdges == 3 || nEdges == 4)
        {
            for (int i = 0; i < nEdges; ++i)
            {
                if (cEdgeTypes[loop.ids[i]] != "circle")
                {
                    isSphere = false;
                    break;
                }
            }

            // Now check circle origins are defined at the same point.
            if (isSphere)
            {
                int originId = edgeMap[loop.ids[0]].ids[1];
                for (int i = 1; i < nEdges; ++i)
                {
                    if (originId != edgeMap[loop.ids[i]].ids[1])
                    {
                        isSphere = false;
                        break;
                    }
                }
            }
        }

        // If we're a sphere, all is good and construct a spherical surface.
        if (isSphere)
        {
            gp_Pnt origin = cPoints[edgeMap[loop.ids[0]].ids[1]];
            gp_Sphere sph;
            sph.SetLocation(origin);
            sph.SetRadius(origin.Distance(cPoints[edgeMap[loop.ids[0]].ids[0]]));
            BRepBuilderAPI_MakeFace makeFace(sph, cWires[surf.ids[0]]);

            ShapeFix_Face sf(makeFace.Face());
            sf.FixOrientation();
            cFaces[surf.id] = sf.Face();
            continue;
        }

        // Attempt to emulate the built-in CAD kernel, which performs a
        // transfinite interpolation.
        bool builtIn = true;
        if (nEdges == 4 && builtIn)
        {
            // Attempt to reconstruct shape from standard Gmsh kernel by
            // creating a custom transfinite surface. So far only 4 edges are
            // supported. This requires a few things:
            //
            // - edges are handles to our OCC curves
            // - clims define the start and end parameter for each curve
            // - vertIds define the vertex IDs for the patch
            // - fwd[i] is true if the CCW-orientated edge i on the patch is the
            //   same direction of curve[i].
            // - verts[i] are the coordinate points of the vertex given by
            //   vertIds[i].
            std::vector<Handle(Geom_Curve)>   edges    (nEdges);
            std::vector<pair<double, double>> clims    (nEdges);
            std::vector<int>                  vertIds  (nEdges);
            std::vector<bool>                 fwd      (nEdges);
            std::vector<gp_Pnt>               verts    (nEdges);

            for (int i = 0; i < nEdges; ++i)
            {
                auto edge = edgeMap[loop.ids[i]];
                auto nextEdge = edgeMap[loop.ids[(i + 1) % nEdges]];

                std::pair<double, double> clim;
                edges[i] = BRep_Tool::Curve(
                    cEdges[loop.ids[i]], clim.first, clim.second);
                clims[i] = clim;

                // Determine orientation.
                fwd[i] = true;
                if ((edge.ids[0] == nextEdge.ids[0]) ||
                    (edge.ids[0] == nextEdge.ids.back()))
                {
                    fwd[i] = false;
                }
                else if ((edge.ids.back() != nextEdge.ids[0]) &&
                         (edge.ids.back() != nextEdge.ids.back()))
                {
                    ASSERTL0(false, "connectivity issue");
                }

                verts[i] = fwd[i] ? cPoints[edge.ids[0]] :
                    cPoints[edge.ids.back()];
            }

            // Create new transfinite surface.
            Handle(Geom_TransfiniteSurface) tf = new Geom_TransfiniteSurface(
                edges, fwd, clims, verts);

            BRepBuilderAPI_MakeFace mkFace(tf, cWires[surf.ids[0]]);
            TopoDS_Face face = mkFace.Face();

            // This is an attempt to figure out the pcurves for each face. For
            // some reason, it isn't working and seems to get nuked by the
            // sewing operation above. Not quite sure why, but it means some
            // OpenCASCADE operations won't work properly. For now the most
            // serious implication is that we need the user to supply a list of
            // void points so that holes in 3D geometries are reproduced.
            if (nEdges == 4 && false)
            {
                BRep_Builder B;
                TopLoc_Location L;
                BRep_Tool::Surface(face, L);
                Handle(Geom2d_Line) e0, e1, e2, e3;
                e0 = new Geom2d_Line(gp_Pnt2d(0, 0), gp_Dir2d( 1,  0));
                e1 = new Geom2d_Line(gp_Pnt2d(1, 0), gp_Dir2d( 0,  1));
                e2 = new Geom2d_Line(gp_Pnt2d(0, 1), gp_Dir2d( 1,  0));
                e3 = new Geom2d_Line(gp_Pnt2d(0, 0), gp_Dir2d( 0,  1));

                B.UpdateEdge(cEdges[loop.ids[0]], e0, tf, L, 0.);
                B.UpdateEdge(cEdges[loop.ids[1]], e1, tf, L, 0.);
                B.UpdateEdge(cEdges[loop.ids[2]], e2, tf, L, 0.);
                B.UpdateEdge(cEdges[loop.ids[3]], e3, tf, L, 0.);
                B.Range(cEdges[loop.ids[0]], face, 0, 1);
                B.Range(cEdges[loop.ids[1]], face, 0, 1);
                B.Range(cEdges[loop.ids[2]], face, 0, 1);
                B.Range(cEdges[loop.ids[3]], face, 0, 1);
            }

            cFaces[surf.id] = face;
            continue;
        }

        // Otherwise we will just go ahead and construct a filling from
        // OpenCASCADE.
        BRepOffsetAPI_MakeFilling fill;
        TopExp_Explorer ex;

        // Add all edges from our line loop.
        for (ex.Init(cWires[surf.ids[0]], TopAbs_EDGE); ex.More(); ex.Next())
        {
            const TopoDS_Edge &shape = TopoDS::Edge(ex.Current());
            fill.Add(shape, GeomAbs_C0);
        }

        fill.Build();
        cFaces[surf.id] = TopoDS::Face(fill.Shape());
    }

    //
    // Construct shells.
    //
    map<unsigned int, TopoDS_Shell> cShells;
    for (auto &sloop : geoFile.surfLoops)
    {
        CheckWarning("surface loop", sloop.id, cShells);
        if (!ContainsIDs("Surface Loop", sloop.id, "surface", sloop.ids, cFaces))
        {
            continue;
        }

        BRep_Builder builder;
        BRepPrim_Builder b(builder);
        TopoDS_Shell shell;
        b.MakeShell(shell);
        for (auto &id : sloop.ids)
        {
            b.AddShellFace(shell, cFaces[id]);
        }
        cShells[sloop.id] = shell;
    }

    map<unsigned int, TopoDS_Shape> cVolumes;
    for (auto &vol : geoFile.volumes)
    {
        for (auto &id : vol.ids)
        {
            id = std::abs(id);
        }

        CheckWarning("volume", vol.id, cVolumes);
        if (!ContainsIDs("Volume", vol.id, "surface loop", vol.ids, cShells))
        {
            continue;
        }

        BRepBuilderAPI_MakeSolid solidMaker;

        for (int i = 0; i < vol.ids.size(); ++i)
        {
            solidMaker.Add(cShells[vol.ids[i]]);
        }

        cVolumes[vol.id] = solidMaker.Solid();
    }

    ASSERTL0(cVolumes.size() == 1,
             "3D .geo file should define exactly one volume.");

    return cVolumes.begin()->second;
}

}
}
