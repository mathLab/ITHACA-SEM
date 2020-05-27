///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessSpherigon.cpp
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
//  Description: Apply Spherigon surface smoothing technique to a 3D mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>

#include <NekMeshUtils/MeshElements/Element.h>

#include "ProcessSpherigon.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;
using namespace Nektar::NekMeshUtils;

#define TOL_BLEND 1.0e-8

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessSpherigon::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "spherigon"), ProcessSpherigon::create);

/**
 * @class ProcessSpherigon
 *
 * This class implements the spherigon surface smoothing technique which
 * is documented in
 *
 *   "The SPHERIGON: A Simple Polygon Patch for Smoothing Quickly your
 *   Polygonal Meshes": P. Volino and N. Magnenat Thalmann, Computer
 *   Animation Proceedings (1998).
 *
 * This implementation works in both a 2D manifold setting (for
 * triangles and quadrilaterals embedded in 3-space) and in a full 3D
 * enviroment (prisms, tetrahedra and hexahedra).
 *
 * No additional information needs to be supplied directly to the module
 * in order for it to work. However, 3D elements do rely on the mapping
 * Mesh::spherigonSurfs to be populated by the relevant input modules.
 *
 * Additionally, since the algorithm assumes normals are supplied which
 * are perpendicular to the true surface defined at each vertex. If
 * these are specified in Mesh::vertexNormals by the input module,
 * better smoothing results can be obtained. Otherwise, normals are
 * estimated by taking the average of all edge/face normals which
 * connect to the vertex.
 */

/**
 * @brief Default constructor.
 */
ProcessSpherigon::ProcessSpherigon(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["N"] =
        ConfigOption(false, "5", "Number of points to add to face edges.");
    m_config["surf"] =
        ConfigOption(false, "-1", "Tag identifying surface to process.");
    m_config["BothTriFacesOnPrism"] = ConfigOption(
        true, "-1", "Curve both triangular faces of prism on boundary.");
    m_config["usenormalfile"] = ConfigOption(
        false, "NoFile", "Use alternative file for Spherigon definition");
    m_config["scalefile"] = ConfigOption(
        false, "1.0", "Apply scaling factor to coordinates in file ");
    m_config["normalnoise"] =
        ConfigOption(false,
                     "NotSpecified",
                     "Add randowm noise to normals of amplitude AMP "
                     "in specified region. input string is "
                     "Amp,xmin,xmax,ymin,ymax,zmin,zmax");
}

/**
 * @brief Destructor.
 */
ProcessSpherigon::~ProcessSpherigon()
{
}

/*
 * @brief Calcuate the normalised cross product \f$ \vec{c} =
 * \vec{a}\times\vec{b} / \| \vec{a}\times\vec{b} \| \f$.
 */
void ProcessSpherigon::UnitCrossProd(Node &a, Node &b, Node &c)
{
    double inv_mag;

    c.m_x = a.m_y * b.m_z - a.m_z * b.m_y;
    c.m_y = a.m_z * b.m_x - a.m_x * b.m_z;
    c.m_z = a.m_x * b.m_y - a.m_y * b.m_x;

    inv_mag = 1.0 / sqrt(c.m_x * c.m_x + c.m_y * c.m_y + c.m_z * c.m_z);

    c.m_x = c.m_x * inv_mag;
    c.m_y = c.m_y * inv_mag;
    c.m_z = c.m_z * inv_mag;
}

/**
 * @brief Calculate the magnitude of the cross product \f$
 * \vec{a}\times\vec{b} \f$.
 */
double ProcessSpherigon::CrossProdMag(Node &a, Node &b)
{
    double tmp1 = a.m_y * b.m_z - a.m_z * b.m_y;
    double tmp2 = a.m_z * b.m_x - a.m_x * b.m_z;
    double tmp3 = a.m_x * b.m_y - a.m_y * b.m_x;
    return sqrt(tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3);
}

/**
 * @brief Calculate the \f$ C^0 \f$ blending function for spherigon
 * implementation.
 *
 * See equation (5) of the paper.
 *
 * @param r  Barycentric coordinate.
 */
double ProcessSpherigon::Blend(double r)
{
    return r * r;
}

/**
 * @brief Calculate the \f$ C^1 \f$ blending function for spherigon
 * implementation.
 *
 * See equation (10) of the paper.
 *
 * @param r      Generalised barycentric coordinates of the point P.
 * @param Q      Vector of vertices denoting this triangle/quad.
 * @param P      Point in the triangle to apply blending to.
 * @param blend  The resulting blending components for each vertex.
 */
void ProcessSpherigon::SuperBlend(vector<double> &r,
                                  vector<Node> &Q,
                                  Node &P,
                                  vector<double> &blend)
{
    vector<double> tmp(r.size());
    double totBlend = 0.0;
    int i;
    int nV = r.size();

    for (i = 0; i < nV; ++i)
    {
        blend[i] = 0.0;
        tmp[i]   = (Q[i] - P).abs2();
    }

    for (i = 0; i < nV; ++i)
    {
        int ip = (i + 1) % nV, im = (i - 1 + nV) % nV;

        if (r[i] > TOL_BLEND && r[i] < (1 - TOL_BLEND))
        {
            blend[i] =
                r[i] * r[i] * (r[im] * r[im] * tmp[im] / (tmp[im] + tmp[i]) +
                               r[ip] * r[ip] * tmp[ip] / (tmp[ip] + tmp[i]));
            totBlend += blend[i];
        }
    }

    for (i = 0; i < nV; ++i)
    {
        blend[i] /= totBlend;
        if (r[i] >= (1 - TOL_BLEND))
        {
            blend[i] = 1.0;
        }
        if (r[i] <= TOL_BLEND)
        {
            blend[i] = 0.0;
        }
    }
}

void ProcessSpherigon::FindNormalFromPlyFile(MeshSharedPtr &plymesh,
                                             map<int,NodeSharedPtr> &surfverts)
{
    int      cnt = 0;
    int      j = 0;
    int      prog=0,cntmin;

    typedef bg::model::point<NekDouble, 3, bg::cs::cartesian> Point;
    typedef pair<Point, unsigned int> PointI;

    int n_neighbs = 1;

    map<int,int>  TreeidtoPlyid;

    //Fill vertex array into tree format
    vector<PointI> dataPts;
    for (auto &it : plymesh->m_vertexSet)
    {
        dataPts.push_back(make_pair(Point(it->m_x, it->m_y, it->m_z), j));
        TreeidtoPlyid[j++] = it->m_id;
    }

    //Build tree
    bgi::rtree<PointI, bgi::rstar<16> > rtree;
    rtree.insert(dataPts.begin(), dataPts.end());

    //Find neipghbours
    for (auto &vIt : surfverts)
    {
        if(m_mesh->m_verbose)
        {
            prog = LibUtilities::PrintProgressbar(cnt,surfverts.size(),
                                                  "Nearest ply verts",prog);
        }

        Point queryPt(vIt.second->m_x, vIt.second->m_y, vIt.second->m_z);
        vector<PointI> result;
        rtree.query(bgi::nearest(queryPt, n_neighbs),
                    std::back_inserter(result));

        cntmin = TreeidtoPlyid[result[0].second];

        ASSERTL1(cntmin < plymesh->m_vertexNormals.size(),
                 "cntmin is out of range");

        m_mesh->m_vertexNormals[vIt.first] =
            plymesh->m_vertexNormals[cntmin];
        ++cnt;
    }
}

/**
 * @brief Generate a set of approximate vertex normals to a surface
 * represented by line segments in 2D and a hybrid
 * triangular/quadrilateral mesh in 3D.
 *
 * This routine approximates the true vertex normals to a surface by
 * averaging the normals of all edges/faces which connect to the
 * vertex. It is better to use the exact surface normals which can be
 * set in Mesh::vertexNormals, but where they are not supplied this
 * routine calculates an approximation for the spherigon implementation.
 *
 * @param el  Vector of elements denoting the surface mesh.
 */
void ProcessSpherigon::GenerateNormals(std::vector<ElementSharedPtr> &el,
                                       MeshSharedPtr &mesh)
{
    for (int i = 0; i < el.size(); ++i)
    {
        ElementSharedPtr e = el[i];

        // Ensure that element is a line, triangle or quad.
        ASSERTL0(e->GetConf().m_e == LibUtilities::eSegment ||
                     e->GetConf().m_e == LibUtilities::eTriangle ||
                     e->GetConf().m_e == LibUtilities::eQuadrilateral,
                 "Spherigon expansions must be lines, triangles or "
                 "quadrilaterals.");

        // Calculate normal for this element.
        int nV = e->GetVertexCount();
        vector<NodeSharedPtr> node(nV);

        for (int j = 0; j < nV; ++j)
        {
            node[j] = e->GetVertex(j);
        }

        Node n;

        if (mesh->m_spaceDim == 3)
        {
            // Create two tangent vectors and take unit cross product.
            Node v1 = *(node[1]) - *(node[0]);
            Node v2 = *(node[2]) - *(node[0]);
            UnitCrossProd(v1, v2, n);
        }
        else
        {
            // Calculate gradient vector and invert.
            Node dx = *(node[1]) - *(node[0]);
            NekDouble diff = dx.abs2();
            dx /= sqrt(diff);
            n.m_x = -dx.m_y;
            n.m_y = dx.m_x;
            n.m_z = 0;
        }

        // Insert face normal into vertex normal list or add to existing
        // value.
        for (int j = 0; j < nV; ++j)
        {
            auto nIt = mesh->m_vertexNormals.find(e->GetVertex(j)->m_id);
            if (nIt == mesh->m_vertexNormals.end())
            {
                mesh->m_vertexNormals[e->GetVertex(j)->m_id] = n;
            }
            else
            {
                nIt->second += n;
            }
        }
    }

    // Normalize resulting vectors.
    for (auto &nIt : mesh->m_vertexNormals)
    {
        Node &n = nIt.second;
        n /= sqrt(n.abs2());
    }

}

/**
 * @brief Perform the spherigon smoothing technique on the mesh.
 */
void ProcessSpherigon::Process()
{
    ASSERTL0(m_mesh->m_spaceDim == 3 || m_mesh->m_spaceDim == 2,
             "Spherigon implementation only valid in 2D/3D.");

    std::unordered_set<int> visitedEdges;

    // First construct vector of elements to process.
    vector<ElementSharedPtr> el;

    if (m_mesh->m_verbose)
    {
        cout << "ProcessSpherigon: Smoothing mesh..." << endl;
    }

    if (m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        // Manifold case - copy expansion dimension.
        el = m_mesh->m_element[m_mesh->m_expDim];
    }
    else if (m_mesh->m_spaceDim == m_mesh->m_expDim)
    {
        // Full 2D or 3D case - iterate over stored edges/faces and
        // create segments/triangles representing those edges/faces.
        set<pair<int, int> >::iterator it;
        vector<int> t;
        t.push_back(0);

        // Construct list of spherigon edges/faces from a tag.
        string surfTag = m_config["surf"].as<string>();
        bool prismTag  = m_config["BothTriFacesOnPrism"].beenSet;

        if (surfTag != "")
        {
            vector<unsigned int> surfs;
            ParseUtils::GenerateSeqVector(surfTag, surfs);
            sort(surfs.begin(), surfs.end());

            m_mesh->m_spherigonSurfs.clear();
            for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
            {
                ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
                int nSurf = m_mesh->m_expDim == 3 ? el->GetFaceCount()
                                                  : el->GetEdgeCount();

                for (int j = 0; j < nSurf; ++j)
                {
                    int bl = el->GetBoundaryLink(j);
                    if (bl == -1)
                    {
                        continue;
                    }

                    ElementSharedPtr bEl =
                        m_mesh->m_element[m_mesh->m_expDim - 1][bl];
                    vector<int> tags = bEl->GetTagList();
                    vector<int> inter;

                    sort(tags.begin(), tags.end());
                    set_intersection(surfs.begin(),
                                     surfs.end(),
                                     tags.begin(),
                                     tags.end(),
                                     back_inserter(inter));

                    if (inter.size() == 1)
                    {
                        m_mesh->m_spherigonSurfs.insert(make_pair(i, j));

                        // Curve other tri face on Prism. Note could be
                        // problem on pyramid when implemented.
                        if (nSurf == 5 && prismTag)
                        {
                            // add other end of prism on boundary for
                            // smoothing
                            int triFace = j == 1 ? 3 : 1;
                            m_mesh->m_spherigonSurfs.insert(
                                make_pair(i, triFace));
                        }
                    }
                }
            }
        }

        if (m_mesh->m_spherigonSurfs.size() == 0)
        {
            cerr << "WARNING: Spherigon surfaces have not been defined "
                 << "-- ignoring smoothing." << endl;
            return;
        }

        if (m_mesh->m_expDim == 3)
        {
            for (it = m_mesh->m_spherigonSurfs.begin();
                 it != m_mesh->m_spherigonSurfs.end();
                 ++it)
            {
                FaceSharedPtr f =
                    m_mesh->m_element[m_mesh->m_expDim][it->first]->GetFace(
                        it->second);
                vector<NodeSharedPtr> nodes = f->m_vertexList;
                LibUtilities::ShapeType eType =
                    (LibUtilities::ShapeType)(nodes.size());
                ElmtConfig conf(eType, 1, false, false);

                // Create 2D element.
                ElementSharedPtr elmt =
                    GetElementFactory().CreateInstance(eType, conf, nodes, t);

                // Copy vertices/edges from face.
                for (int i = 0; i < f->m_vertexList.size(); ++i)
                {
                    elmt->SetVertex(i, f->m_vertexList[i]);
                }
                for (int i = 0; i < f->m_edgeList.size(); ++i)
                {
                    elmt->SetEdge(i, f->m_edgeList[i]);
                }

                el.push_back(elmt);
            }
        }
        else
        {
            for (it = m_mesh->m_spherigonSurfs.begin();
                 it != m_mesh->m_spherigonSurfs.end();
                 ++it)
            {
                EdgeSharedPtr edge =
                    m_mesh->m_element[m_mesh->m_expDim][it->first]->GetEdge(
                        it->second);
                vector<NodeSharedPtr> nodes;
                LibUtilities::ShapeType eType = LibUtilities::eSegment;
                ElmtConfig conf(eType, 1, false, false);

                nodes.push_back(edge->m_n1);
                nodes.push_back(edge->m_n2);

                // Create 2D element.
                ElementSharedPtr elmt =
                    GetElementFactory().CreateInstance(eType, conf, nodes, t);

                // Copy vertices/edges from original element.
                elmt->SetVertex(0, nodes[0]);
                elmt->SetVertex(1, nodes[1]);
                elmt->SetEdge(
                    0,
                    m_mesh->m_element[m_mesh->m_expDim][it->first]->GetEdge(
                        it->second));
                el.push_back(elmt);
            }
        }
    }
    else
    {
        ASSERTL0(false, "Spherigon expansions must be 2/3 dimensional");
    }

    // See if vertex normals have been generated. If they have not,
    // approximate them by summing normals of surrounding elements.
    bool normalsGenerated = false;

    // Read Normal file if one exists
    std::string normalfile = m_config["usenormalfile"].as<string>();
    if (normalfile.compare("NoFile") != 0)
    {
        NekDouble scale = m_config["scalefile"].as<NekDouble>();

        if (m_mesh->m_verbose)
        {
            cout << "Inputing normal file: " << normalfile
                 << " with scaling of " << scale << endl;
        }

        ifstream inplyTmp;
        io::filtering_istream inply;
        InputPlySharedPtr plyfile;

        inplyTmp.open(normalfile.c_str());
        ASSERTL0(inplyTmp,
                 string("Could not open input ply file: ") + normalfile);

        inply.push(inplyTmp);

        MeshSharedPtr m = std::shared_ptr<Mesh>(new Mesh());
        plyfile = std::shared_ptr<InputPly>(new InputPly(m));
        plyfile->ReadPly(inply, scale);
        plyfile->ProcessVertices();

        MeshSharedPtr plymesh = plyfile->GetMesh();
        if (m_mesh->m_verbose)
        {
            cout << "\t Generating ply normals" << endl;
        }
        GenerateNormals(plymesh->m_element[plymesh->m_expDim], plymesh);

        // finally find nearest vertex and set normal to mesh surface file
        // normal.  probably should have a hex tree search ?
        Node minx(0, 0.0, 0.0, 0.0), tmp, tmpsav;
        NodeSet::iterator it;
        map<int, NodeSharedPtr>::iterator vIt;
        map<int, NodeSharedPtr> surfverts;

        // make a map of normal vertices to visit based on elements el
        for (int i = 0; i < el.size(); ++i)
        {
            ElementSharedPtr e = el[i];
            int nV = e->GetVertexCount();
            for (int j = 0; j < nV; ++j)
            {
                int id        = e->GetVertex(j)->m_id;
                surfverts[id] = e->GetVertex(j);
            }
        }


        if (m_mesh->m_verbose)
        {
            cout << "\t Processing surface normals " << endl;
        }

        // loop over all element in ply mesh and determine
        // and set normal to nearest point in ply mesh
        FindNormalFromPlyFile(plymesh,surfverts);

        normalsGenerated = true;
    }
    else if (m_mesh->m_vertexNormals.size() == 0)
    {
        GenerateNormals(el, m_mesh);
        normalsGenerated = true;
    }

    // See if we should add noise to normals
    std::string normalnoise = m_config["normalnoise"].as<string>();
    if (normalnoise.compare("NotSpecified") != 0)
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateVector(normalnoise, values),
                 "Failed to interpret normal noise string");

        int nvalues   = values.size() / 2;
        NekDouble amp = values[0];

        if (m_mesh->m_verbose)
        {
            cout << "\t adding noise to normals of amplitude " << amp
                 << " in range: ";
            for (int i = 0; i < nvalues; ++i)
            {
                cout << values[2 * i + 1] << "," << values[2 * i + 2] << " ";
            }
            cout << endl;
        }

        map<int, NodeSharedPtr>::iterator vIt;
        map<int, NodeSharedPtr> surfverts;

        // make a map of normal vertices to visit based on elements el
        for (int i = 0; i < el.size(); ++i)
        {
            ElementSharedPtr e = el[i];
            int nV = e->GetVertexCount();
            for (int j = 0; j < nV; ++j)
            {
                int id        = e->GetVertex(j)->m_id;
                surfverts[id] = e->GetVertex(j);
            }
        }

        for (vIt = surfverts.begin(); vIt != surfverts.end(); ++vIt)
        {
            bool AddNoise = false;

            for (int i = 0; i < nvalues; ++i)
            {
                // check to see if point is in range
                switch (nvalues)
                {
                    case 1:
                    {
                        if (((vIt->second)->m_x > values[2 * i + 1]) &&
                            ((vIt->second)->m_x < values[2 * i + 2]))
                        {
                            AddNoise = true;
                        }
                    }
                    break;
                    case 2:
                    {
                        if (((vIt->second)->m_x > values[2 * i + 1]) &&
                            ((vIt->second)->m_x < values[2 * i + 2]) &&
                            ((vIt->second)->m_y > values[2 * i + 3]) &&
                            ((vIt->second)->m_y < values[2 * i + 4]))
                        {
                            AddNoise = true;
                        }
                    }
                    break;
                    case 3:
                    {
                        if (((vIt->second)->m_x > values[2 * i + 1]) &&
                            ((vIt->second)->m_x < values[2 * i + 2]) &&
                            ((vIt->second)->m_y > values[2 * i + 3]) &&
                            ((vIt->second)->m_y < values[2 * i + 4]) &&
                            ((vIt->second)->m_z > values[2 * i + 5]) &&
                            ((vIt->second)->m_z < values[2 * i + 6]))

                        {
                            AddNoise = true;
                        }
                    }
                    break;
                }

                if (AddNoise)
                {
                    // generate random unit vector;
                    Node rvec(0, rand(), rand(), rand());
                    rvec *= values[0] / sqrt(rvec.abs2());

                    Node normal = m_mesh->m_vertexNormals[vIt->first];

                    normal += rvec;
                    normal /= sqrt(normal.abs2());

                    m_mesh->m_vertexNormals[vIt->first] = normal;
                }
            }
        }
    }

    // Allocate storage for interior points.
    int nq    = m_config["N"].as<int>();
    int nquad = m_mesh->m_spaceDim == 3 ? nq * nq : nq;
    Array<OneD, NekDouble> x(nq * nq);
    Array<OneD, NekDouble> y(nq * nq);
    Array<OneD, NekDouble> z(nq * nq);

    Array<OneD, NekDouble> xc(nq * nq);
    Array<OneD, NekDouble> yc(nq * nq);
    Array<OneD, NekDouble> zc(nq * nq);

    ASSERTL0(nq > 2, "Number of points must be greater than 2.");

    LibUtilities::BasisKey B0(
        LibUtilities::eOrtho_A,
        nq,
        LibUtilities::PointsKey(nq, LibUtilities::eGaussLobattoLegendre));
    LibUtilities::BasisKey B1(
        LibUtilities::eOrtho_B,
        nq,
        LibUtilities::PointsKey(nq, LibUtilities::eGaussRadauMAlpha1Beta0));
    StdRegions::StdNodalTriExpSharedPtr stdtri =
        MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
            B0, B1, LibUtilities::eNodalTriElec);

    Array<OneD, NekDouble> xnodal(nq * (nq + 1) / 2), ynodal(nq * (nq + 1) / 2);
    stdtri->GetNodalPoints(xnodal, ynodal);

    int edgeMap[3][4][2] = {
        {{0, 1}, {-1, -1}, {-1, -1}, {-1, -1}},                         // seg
        {{0, 1}, {nq - 1, nq}, {nq * (nq - 1), -nq}, {-1, -1}},         // tri
        {{0, 1}, {nq - 1, nq}, {nq * nq - 1, -1}, {nq * (nq - 1), -nq}} // quad
    };

    int vertMap[3][4][2] = {
        {{0, 1}, {0, 0}, {0, 0}, {0, 0}}, // seg
        {{0, 1}, {1, 2}, {2, 3}, {0, 0}}, // tri
        {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, // quad
    };

    for (int i = 0; i < el.size(); ++i)
    {
        // Construct a Nektar++ element to obtain coordinate points
        // inside the element. TODO: Add options for various
        // nodal/tensor point distributions + number of points to add.
        ElementSharedPtr e = el[i];

        LibUtilities::BasisKey B2(
            LibUtilities::eModified_A,
            nq,
            LibUtilities::PointsKey(nq, LibUtilities::eGaussLobattoLegendre));

        if (e->GetConf().m_e == LibUtilities::eSegment)
        {
            SpatialDomains::SegGeomSharedPtr geom =
                std::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                    e->GetGeom(m_mesh->m_spaceDim));
            LocalRegions::SegExpSharedPtr seg =
                MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(B2,
                                                                       geom);
            seg->GetCoords(x, y, z);
            nquad = nq;
        }
        else if (e->GetConf().m_e == LibUtilities::eTriangle)
        {
            SpatialDomains::TriGeomSharedPtr geom =
                std::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                    e->GetGeom(3));
            LocalRegions::NodalTriExpSharedPtr tri =
                MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(
                    B0, B1, LibUtilities::eNodalTriElec, geom);

            Array<OneD, NekDouble> coord(2);
            tri->GetCoords(xc, yc, zc);
            nquad = nq * (nq + 1) / 2;

            for (int j = 0; j < nquad; ++j)
            {
                coord[0] = xnodal[j];
                coord[1] = ynodal[j];
                x[j]     = stdtri->PhysEvaluate(coord, xc);
            }

            for (int j = 0; j < nquad; ++j)
            {
                coord[0] = xnodal[j];
                coord[1] = ynodal[j];
                y[j]     = stdtri->PhysEvaluate(coord, yc);
            }

            for (int j = 0; j < nquad; ++j)
            {
                coord[0] = xnodal[j];
                coord[1] = ynodal[j];
                z[j]     = stdtri->PhysEvaluate(coord, zc);
            }
        }
        else if (e->GetConf().m_e == LibUtilities::eQuadrilateral)
        {
            SpatialDomains::QuadGeomSharedPtr geom =
                std::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                    e->GetGeom(3));
            LocalRegions::QuadExpSharedPtr quad =
                MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                    B2, B2, geom);
            quad->GetCoords(x, y, z);
            nquad = nq * nq;
        }
        else
        {
            ASSERTL0(false, "Unknown expansion type.");
        }

        // Zero z-coordinate in 2D.
        if (m_mesh->m_spaceDim == 2)
        {
            Vmath::Zero(nquad, z, 1);
        }

        // Find vertex normals.
        int nV = e->GetVertexCount();
        vector<Node> v, vN;
        for (int j = 0; j < nV; ++j)
        {
            v.push_back(*(e->GetVertex(j)));
            ASSERTL1(m_mesh->m_vertexNormals.count(v[j].m_id) != 0,
                     "Normal has not been defined");
            vN.push_back(m_mesh->m_vertexNormals[v[j].m_id]);
        }

        vector<Node> tmp(nV);
        vector<NekDouble> r(nV);
        vector<Node> K(nV);
        vector<Node> Q(nV);
        vector<Node> Qp(nV);
        vector<NekDouble> blend(nV);
        vector<Node> out(nquad);

        // Calculate segment length for 2D spherigon routine.
        NekDouble segLength = sqrt((v[0] - v[1]).abs2());

        // Perform Spherigon method to smooth manifold.
        for (int j = 0; j < nquad; ++j)
        {
            Node P(0, x[j], y[j], z[j]);
            Node N(0, 0, 0, 0);

            // Calculate generalised barycentric coordinates r[] and the
            // Phong normal N = vN . r for this point of the element.
            if (m_mesh->m_spaceDim == 2)
            {
                // In 2D the coordinates are given by a ratio of the
                // segment length to the distance from one of the
                // endpoints.
                r[0] = sqrt((P - v[0]).abs2()) / segLength;
                r[0] = max(min(1.0, r[0]), 0.0);
                r[1] = 1.0 - r[0];

                // Calculate Phong normal.
                N = vN[0] * r[0] + vN[1] * r[1];
            }
            else if (m_mesh->m_spaceDim == 3)
            {
                for (int k = 0; k < nV; ++k)
                {
                    tmp[k] = P - v[k];
                }

                // Calculate generalized barycentric coordinate system
                // (see equation 6 of paper).
                NekDouble weight = 0.0;
                for (int k = 0; k < nV; ++k)
                {
                    r[k] = 1.0;
                    for (int l = 0; l < nV - 2; ++l)
                    {
                        r[k] *= CrossProdMag(tmp[(k + l + 1) % nV],
                                             tmp[(k + l + 2) % nV]);
                    }
                    weight += r[k];
                }

                // Calculate Phong normal (equation 1).
                for (int k = 0; k < nV; ++k)
                {
                    r[k] /= weight;
                    N += vN[k] * r[k];
                }
            }

            // Normalise Phong normal.
            N /= sqrt(N.abs2());

            for (int k = 0; k < nV; ++k)
            {
                // Perform steps denoted in equations 2, 3, 8 for C1
                // smoothing.
                NekDouble tmp1;
                K[k]  = P + N * ((v[k] - P).dot(N));
                tmp1  = (v[k] - K[k]).dot(vN[k]) / (1.0 + N.dot(vN[k]));
                Q[k]  = K[k] + N * tmp1;
                Qp[k] = v[k] - N * ((v[k] - P).dot(N));
            }

            // Apply C1 blending function to the surface. TODO: Add
            // option to do (more efficient) C0 blending function.
            SuperBlend(r, Qp, P, blend);
            P.m_x = P.m_y = P.m_z = 0.0;

            // Apply blending (equation 4).
            for (int k = 0; k < nV; ++k)
            {
                P += Q[k] * blend[k];
            }

            if ((boost::math::isnan)(P.m_x) || (boost::math::isnan)(P.m_y) ||
                (boost::math::isnan)(P.m_z))
            {
                ASSERTL0(false,
                         "spherigon point is a nan. Check to see if "
                         "ply file is correct if using input normal file");
            }
            else
            {
                out[j] = P;
            }
        }

        // Push nodes into lines - TODO: face interior nodes.
        // offset = 0 (seg), 1 (tri) or 2 (quad)
        int offset = (int)e->GetConf().m_e - 2;

        for (int edge = 0; edge < e->GetEdgeCount(); ++edge)
        {
            auto eIt = visitedEdges.find(e->GetEdge(edge)->m_id);
            if (eIt == visitedEdges.end())
            {
                bool reverseEdge =
                    !(v[vertMap[offset][edge][0]] == *(e->GetEdge(edge)->m_n1));

                // Clear existing curvature.
                e->GetEdge(edge)->m_edgeNodes.clear();

                if (e->GetConf().m_e != LibUtilities::eTriangle)
                {
                    for (int j = 1; j < nq - 1; ++j)
                    {
                        int v = edgeMap[offset][edge][0] +
                                j * edgeMap[offset][edge][1];
                        e->GetEdge(edge)->m_edgeNodes.push_back(
                            NodeSharedPtr(new Node(out[v])));
                    }
                }
                else
                {
                    for (int j = 0; j < nq - 2; ++j)
                    {
                        int v = 3 + edge * (nq - 2) + j;
                        e->GetEdge(edge)->m_edgeNodes.push_back(
                            NodeSharedPtr(new Node(out[v])));
                    }
                }

                if (reverseEdge)
                {
                    reverse(e->GetEdge(edge)->m_edgeNodes.begin(),
                            e->GetEdge(edge)->m_edgeNodes.end());
                }

                e->GetEdge(edge)->m_curveType =
                    LibUtilities::eGaussLobattoLegendre;

                visitedEdges.insert(e->GetEdge(edge)->m_id);
            }
        }

        // Add face nodes in manifold and full 3D case, but not for 2D.
        if (m_mesh->m_spaceDim == 3)
        {
            vector<NodeSharedPtr> volNodes;

            if (e->GetConf().m_e == LibUtilities::eQuadrilateral)
            {
                volNodes.resize((nq - 2) * (nq - 2));
                for (int j = 1; j < nq - 1; ++j)
                {
                    for (int k = 1; k < nq - 1; ++k)
                    {
                        int v = j * nq + k;
                        volNodes[(j - 1) * (nq - 2) + (k - 1)] =
                            NodeSharedPtr(new Node(out[v]));
                    }
                }
            }
            else
            {
                for (int j = 3 + 3 * (nq - 2); j < nquad; ++j)
                {
                    volNodes.push_back(NodeSharedPtr(new Node(out[j])));
                }
            }

            e->SetVolumeNodes(volNodes);
        }
    }

    // Copy face nodes back into 3D element faces
    if (m_mesh->m_expDim == 3)
    {
        int elmt = 0;
        for (auto &it : m_mesh->m_spherigonSurfs)
        {
            FaceSharedPtr f =
                m_mesh->m_element[m_mesh->m_expDim][it.first]->GetFace(
                    it.second);

            f->m_faceNodes = el[elmt++]->GetVolumeNodes();
            f->m_curveType = f->m_vertexList.size() == 3
                                 ? LibUtilities::eNodalTriElec
                                 : LibUtilities::eGaussLobattoLegendre;
        }
    }

    if (normalsGenerated)
    {
        m_mesh->m_vertexNormals.clear();
    }
}
}
}
