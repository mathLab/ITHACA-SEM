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

#include <NekMeshUtils/CADSystem/OCE/CADSystemOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADVertOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADCurveOCE.h>
#include <NekMeshUtils/CADSystem/OCE/CADSurfOCE.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

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
        //not a naca profile behave normally
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

    // faces and verts can be extracted straight from shape
    TopTools_IndexedMapOfShape mapOfVerts, mapOfFaces;
    TopExp::MapShapes(shape, TopAbs_VERTEX, mapOfVerts);
    TopExp::MapShapes(shape, TopAbs_FACE, mapOfFaces);

    // edges need to be built from loops around faces to elimiate degen and
    // hanging edges
    TopTools_IndexedMapOfShape mapOfEdges;

    // build map of verticies
    for (int i = 1; i <= mapOfVerts.Extent(); i++)
    {
        TopoDS_Shape v = mapOfVerts.FindKey(i);
        AddVert(i, v);
    }

    // For each face of the geometry, get the local edges which bound it. If
    // they are valid (their type != 7), then add them to an edge map. This
    // filters out the dummy edges which OCC uses.
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape localEdges;
        TopExp::MapShapes(face, TopAbs_EDGE, localEdges);

        for (int j = 1; j <= localEdges.Extent(); j++)
        {
            TopoDS_Shape edge       = localEdges.FindKey(j);
            BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
            if (curve.GetType() != 7)
            {
                if (!(mapOfEdges.Contains(edge)))
                {
                    mapOfEdges.Add(edge);
                }
            }
        }
    }

    // from id of curve to list of ids of surfs
    map<int, vector<int> > adjsurfmap;

    // Adds edges to our type and map
    for (int i = 1; i <= mapOfEdges.Extent(); i++)
    {
        TopoDS_Shape edge = mapOfEdges.FindKey(i);
        TopoDS_Vertex fv =
            TopExp::FirstVertex(TopoDS::Edge(edge), Standard_True);
        TopoDS_Vertex lv =
            TopExp::LastVertex(TopoDS::Edge(edge), Standard_True);

        if (edge.Orientation() == 0)
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(fv),
                     mapOfVerts.FindIndex(lv));
        }
        else
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(lv),
                     mapOfVerts.FindIndex(fv));
        }
    }

    // For each face, examine all the wires (i.e. bounding loops) and
    // investigates the loop. Using this information, connectivity is determined
    // and edges are associated with surfaces.
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape mapOfWires;
        TopExp::MapShapes(face, TopAbs_WIRE, mapOfWires);

        // this pice of code does an idiot check on the loops to make sure
        // they dont cross or touch
        if (mapOfWires.Extent() > 1)
        {
            TopoDS_Wire ow = BRepTools::OuterWire(TopoDS::Face(face));

            vector<TopoDS_Shape> wirefacecuts;
            vector<gp_Pnt> centersofcutfaces;

            for (int j = 1; j <= mapOfWires.Extent(); j++)
            {
                TopoDS_Shape wire = mapOfWires.FindKey(j);

                if (wire != ow)
                {
                    BRepBuilderAPI_MakeFace build(
                        BRep_Tool::Surface(TopoDS::Face(face)), 1e-7);
                    build.Add(TopoDS::Wire(wire));
                    TopoDS_Shape newface = build.Shape();
                    wirefacecuts.push_back(newface);
                    BRepAdaptor_Surface b =
                        BRepAdaptor_Surface(TopoDS::Face(newface));
                    NekDouble u, v;
                    u = (b.LastUParameter() - b.FirstUParameter()) / 2.0;
                    v = (b.LastVParameter() - b.FirstVParameter()) / 2.0;
                    centersofcutfaces.push_back(b.Value(u, v));
                }
            }
            for (int j = 0; j < wirefacecuts.size(); j++)
            {
                for (int k = 0; k < wirefacecuts.size(); k++)
                {
                    if (j == k)
                        continue;

                    /// TODO fix this test
                    BRepClass_FaceClassifier fc(TopoDS::Face(wirefacecuts[j]),
                                                centersofcutfaces[k], 1e-7);
                    // ASSERTL0(fc.State() == 1, "Internal face loops make this
                    // cad impossible to mesh");
                }
            }
        }

        vector<EdgeLoop> edgeloops;

        // now we acutally analyse the loops for cad building
        for (int j = 1; j <= mapOfWires.Extent(); j++)
        {
            EdgeLoop edgeloop;

            TopoDS_Shape wire = mapOfWires.FindKey(j);

            ShapeAnalysis_Wire wiretest(TopoDS::Wire(wire), TopoDS::Face(face),
                                        1E-6);

            // calculate the center of the wire
            GProp_GProps massProps;
            BRepGProp::SurfaceProperties(wire, massProps);
            gp_Pnt gPt = massProps.CentreOfMass();
            // alternative locuv methods
            ShapeAnalysis_Surface sas(BRep_Tool::Surface(TopoDS::Face(face)));
            sas.SetDomain(
                BRepAdaptor_Surface(TopoDS::Face(face)).FirstUParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).LastUParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).FirstVParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).LastVParameter());
            gp_Pnt2d p2 = sas.ValueOfUV(gPt, 1e-7);
            Array<OneD, NekDouble> cen(2);
            cen[0]          = p2.X();
            cen[1]          = p2.Y();
            edgeloop.center = cen;

            BRepTools_WireExplorer exp;

            exp.Init(TopoDS::Wire(wire));

            while (exp.More())
            {
                TopoDS_Shape edge = exp.Current();

                if (mapOfEdges.Contains(edge))
                {
                    int e = mapOfEdges.FindIndex(edge);
                    edgeloop.edges.push_back(m_curves[e]);
                    edgeloop.edgeo.push_back(exp.Orientation());
                    adjsurfmap[e].push_back(i);
                }

                exp.Next();
            }

            edgeloops.push_back(edgeloop);
        }

        AddSurf(i, face, edgeloops);
    }

    // attempts to identify properties of the vertex on the degen edge
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

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
    for (map<int, vector<int> >::iterator it = adjsurfmap.begin();
         it != adjsurfmap.end(); it++)
    {
        if(!m_2d)
        {
            ASSERTL0(it->second.size() == 2, "no three curve surfaces");
        }

        vector<CADSurfSharedPtr> sfs;
        for (int i = 0; i < it->second.size(); i++)
        {
            sfs.push_back(m_surfs[it->second[i]]);
        }
        m_curves[it->first]->SetAdjSurf(sfs);
    }
    return true;
}

void CADSystemOCE::AddVert(int i, TopoDS_Shape in)
{
    CADVertSharedPtr newVert = GetCADVertFactory().CreateInstance(key);

    boost::static_pointer_cast<CADVertOCE>(newVert)->Initialise(i, in);

    m_verts[i] = newVert;
}

void CADSystemOCE::AddCurve(int i, TopoDS_Shape in, int fv, int lv)
{
    CADCurveSharedPtr newCurve = GetCADCurveFactory().CreateInstance(key);
    boost::static_pointer_cast<CADCurveOCE>(newCurve)->Initialise(i, in);

    vector<CADVertSharedPtr> vs;
    vs.push_back(m_verts[fv]);
    vs.push_back(m_verts[lv]);
    m_curves[i] = newCurve;
    m_curves[i]->SetVert(vs);
}

void CADSystemOCE::AddSurf(int i, TopoDS_Shape in, vector<EdgeLoop> ein)
{
    CADSurfSharedPtr newSurf = GetCADSurfFactory().CreateInstance(key);
    boost::static_pointer_cast<CADSurfOCE>(newSurf)->Initialise(i, in, ein);
    m_surfs[i] = newSurf;

    if (in.Orientation() == 0)
    {
        m_surfs[i]->SetReverseNomral();
    }

    int tote = 0;
    for (int i = 0; i < ein.size(); i++)
    {
        tote += ein[i].edges.size();
    }

    ASSERTL0(tote != 1, "cannot handle periodic curves");
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
        CADCurveSharedPtr c         = GetCurve(i);
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

    int n = boost::lexical_cast<int>(naca);
    NekDouble T = (n%100) / 100.0;
    n/=100;
    NekDouble P = (n%10)/10.0;
    n/=10;
    NekDouble M = (n%10)/100.0;

    int np = 25;

    Array<OneD, NekDouble> xc(np);
    NekDouble dtheta = M_PI/(np-1);
    for(int i = 0; i < np; i++)
    {
        xc[i] = (1.0 - cos(i*dtheta)) / 2.0;
    }

    Array<OneD, NekDouble> yc(np), dyc(np);
    for(int i = 0; i < np; i++)
    {
        if(xc[i] < P)
        {
            yc[i] = M / P / P * (2.0 * P * xc[i] - xc[i] * xc[i]);
            dyc[i] = 2.0 * M / P / P * (P - xc[i]);
        }
        else
        {
            yc[i]  = M / (1.0 - P) / (1.0 - P) * (
                         1.0 - 2.0 * P + 2.0 * P * xc[i] - xc[i] * xc[i]);
            dyc[i] = 2.0 * M / (1.0 - P) / (1.0 - P) * (P - xc[i]);
        }
    }

    Array<OneD, NekDouble> yt(np);
    for(int i = 0; i < np; i++)
    {
        yt[i] = T / 0.2 * ( 0.2969 * sqrt(xc[i])
                           -0.1260 * xc[i]
                           -0.3516 * xc[i] * xc[i]
                           +0.2843 * xc[i] * xc[i] * xc[i]
                           -0.1015 * xc[i] * xc[i] * xc[i] * xc[i]);
    }

    Array<OneD, NekDouble> x(2*np-1), y(2*np-1);
    int l = 0;
    for(int i = np-1; i >= 0; i--, l++)
    {
        NekDouble theta = atan(dyc[i]);

        x[l] = xc[i] - yt[i] * sin(theta);
        y[l] = yc[i] + yt[i] * cos(theta);

    }
    for(int i = 1; i < np; i++)
    {
        NekDouble theta = atan(dyc[i]);

        x[i+np-1] = xc[i] + yt[i] * sin(theta);
        y[i+np-1] = yc[i] - yt[i] * cos(theta);
    }

    TColgp_Array1OfPnt pointArray(0,2*np-2);

    for(int i = 0; i < 2*np-1; i++)
    {
        pointArray.SetValue(i,gp_Pnt(x[i]*1000.0,y[i]*1000.0,0.0));
    }

    GeomAPI_PointsToBSpline spline(pointArray);
    Handle(Geom_BSplineCurve) curve = spline.Curve();

    BRepBuilderAPI_MakeEdge areoEdgeBuilder(curve);
    TopoDS_Edge aeroEdge = areoEdgeBuilder.Edge();
    BRepBuilderAPI_MakeEdge aeroTEBuilder(
        gp_Pnt(x[0]*1000.0, y[0]*1000.0, 0.0),
        gp_Pnt(x[2*np-2]*1000.0, y[2*np-2]*1000.0, 0.0));
    TopoDS_Edge TeEdge = aeroTEBuilder.Edge();

    BRepBuilderAPI_MakeWire aeroWireBuilder(aeroEdge, TeEdge);
    TopoDS_Wire aeroWire = aeroWireBuilder.Wire();

    gp_Trsf transform;
    gp_Ax1 rotAx(gp_Pnt(500.0,0.0,0.0),gp_Dir(gp_Vec(0.0,0.0,-1.0)));
    transform.SetRotation(rotAx, data[4]/180.0*M_PI);
    TopLoc_Location mv(transform);
    aeroWire.Move(mv);

    BRepBuilderAPI_MakeEdge domInlBuilder(gp_Pnt(data[0]*1000.0,data[1]*1000.0,0.0),
                                          gp_Pnt(data[0]*1000.0,data[3]*1000.0,0.0));
    TopoDS_Edge inlEdge = domInlBuilder.Edge();

    BRepBuilderAPI_MakeEdge domTopBuilder(gp_Pnt(data[0]*1000.0,data[3]*1000.0,0.0),
                                          gp_Pnt(data[2]*1000.0,data[3]*1000.0,0.0));
    TopoDS_Edge topEdge = domTopBuilder.Edge();

    BRepBuilderAPI_MakeEdge domOutBuilder(gp_Pnt(data[2]*1000.0,data[3]*1000.0,0.0),
                                          gp_Pnt(data[2]*1000.0,data[1]*1000.0,0.0));
    TopoDS_Edge outEdge = domOutBuilder.Edge();

    BRepBuilderAPI_MakeEdge domBotBuilder(gp_Pnt(data[2]*1000.0,data[1]*1000.0,0.0),
                                          gp_Pnt(data[0]*1000.0,data[1]*1000.0,0.0));
    TopoDS_Edge botEdge = domBotBuilder.Edge();

    BRepBuilderAPI_MakeWire domWireBuilder(inlEdge, topEdge, outEdge, botEdge);
    TopoDS_Wire domWire = domWireBuilder.Wire();

    BRepBuilderAPI_MakeFace face(domWire, true);
    face.Add(aeroWire);

    return face.Face();
}

}
}
