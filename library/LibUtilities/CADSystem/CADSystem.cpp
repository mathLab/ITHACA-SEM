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

/**
 * @brief Return the name of the CAD system.
 */
string CADSystem::GetName()
{
    return m_name;
}

/**
 * @brief Reports basic properties to screen.
 */
void CADSystem::Report()
{
    cout << endl << "CAD report:" << endl;
    cout << "\tCAD has: " << m_curves.size() << " curves." << endl;
    cout << "\tCAD has: " << m_surfs.size() << " surfaces." << endl;
    cout << "\tCAD Euler-Poincaré characteristic: " << m_epc << endl;
}

/**
 * @brief Returns bounding box of the domain.
 *
 * Gets the bounding box of the domain by considering the start and end
 * points of each curve in the geometry.
 *
 * @return Array with 6 entries: xmin, xmax, ymin, ymax, zmin and zmax.
 */
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

/**
 * @brief Initialises CAD and makes surface and curve maps.
 *
 * @return true if completed successfully
 */
bool CADSystem::LoadCAD()
{
    if (!boost::filesystem::exists(m_name.c_str()))
    {
        return false;
    }

    string ext;
    size_t pos = m_name.find(".");
    ext = m_name.substr(pos);

    TopoDS_Shape shape;

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

    // From OpenCascade maps calculate Euler-Poincare number.
    TopTools_IndexedMapOfShape fc, vc, ec;
    TopExp::MapShapes(shape, TopAbs_FACE, fc);
    TopExp::MapShapes(shape, TopAbs_EDGE, ec);
    TopExp::MapShapes(shape, TopAbs_VERTEX, vc);

    m_epc = vc.Extent() - ec.Extent() + fc.Extent();

    TopTools_IndexedMapOfShape mapOfFaces;
    TopTools_IndexedMapOfShape mapOfEdges;
    TopExp::MapShapes(shape, TopAbs_FACE, mapOfFaces);

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
        AddCurve(i, edge);
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

void CADSystem::AddCurve(int i, TopoDS_Shape in)
{
    CADCurveSharedPtr newCurve = MemoryManager<CADCurve>::
                                            AllocateSharedPtr(i,in);
    m_curves[i] = newCurve;
}

void CADSystem::AddSurf(int i, TopoDS_Shape in,
                        std::vector<std::vector<std::pair<int,int> > > ein)
{
    CADSurfSharedPtr newSurf = MemoryManager<CADSurf>::
                                            AllocateSharedPtr(i,in,ein);
    m_surfs[i] = newSurf;
}

}
}
