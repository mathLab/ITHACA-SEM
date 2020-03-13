#include <cstdlib>
#include <string>
#include <cassert>

#include <iostream>
#include <fstream>
#include <vector>

#include <boost/lexical_cast.hpp>

using namespace std;


void print_usage_info(char* binary_name)
{
        //                     argv:      1        2        3   4   5   6   7         8
        cerr << "Usage "<<binary_name<<" meshtype  splits   nx  ny  sx  sy  nummodes  output_file.xml\n";
        cerr << "where 'meshtype' is\n";
        cerr << "        = 0 for regular diamond mesh\n";
        cerr << "        = 1 for regular diamond mesh with even vertical stripes of diamonds not being split\n";
        cerr << "        = 2 for regular diamond mesh with even horizontal stripes of diamonds shifted one quad right\n";
        cerr << "        = 3 for a single quad split into triangles symmetrically relative to its central point\n";
        cerr << "      'splits' is:\n";
        cerr << "        - for 'meshtype' = 0 - a number of points that splits every edge on each even horizontal line,\n";
        cerr << "        - for 'meshtype' = 1 - a number of points that splits every odd edge on each even horizontal line,\n";
        cerr << "        - for 'meshtype' = 2 - a number of points that splits every interior edge on each even horizontal line,\n";
        cerr << "        - for 'meshtype' = 3 - a number of points splitting every symmetrical part of boundary,\n";
        cerr << "      nx   is the number of points in x direction that forms quadrilateral grid skeleton,\n";
        cerr << "      ny   is the number of points in y direction that forms quadrilateral grid skeleton,\n";
        cerr << "      sx   is the coordinate scaling in x direction,\n";
        cerr << "      sy   is the coordinate scaling in y direction,\n";
        cerr << "      nummodes is the number of boundary modes.\n";
        cerr << "For meshtype in {0,1,2} it generates regular mesh with triangles filling quadrilateral grid (aka skeleton)\n";
        cerr << "in the way that forms vertices of different valence.\n";
        cerr << "All vertex coordinates of quadrilateral grid skeleton are evenly distributed within\n";
        cerr << "unit square but then scaled by means of XSCALE and YSCALE attributes of VERTEX section.\n";
}




void  PrintConditions(ofstream& output);


struct Vertex
{
    Vertex() {}
    Vertex(int id, double x, double y):
        m_id(id),
        m_x(x),
        m_y(y)
    {}
    void init(int id, double x, double y)
    {
        m_id = id;
        m_x  = x;
        m_y  = y;
    }

    int m_id;
    double m_x;
    double m_y;
};

struct Segment
{
    Segment() {}
    Segment(int id, Vertex v1, Vertex v2):
        m_id(id),
        m_v1(v1),
        m_v2(v2)
    {}
    void init(int id, Vertex v1, Vertex v2)
    {
        m_id = id;
        m_v1 = v1;
        m_v2 = v2;
    }

    int     m_id;
    Vertex  m_v1, m_v2;
};

struct Triangle
{
    Triangle() {}
    Triangle(Segment s1, Segment s2, Segment s3):
        m_s1(s1),
        m_s2(s2),
        m_s3(s3)
    {}

    Segment m_s1, m_s2, m_s3;
};

struct Mesh
{
    vector<Vertex>   vert;
    vector<Segment>  seg;
    vector<Triangle> tri;

    // boundary
    vector<Segment>  south, north, east, west;
};

enum MeshType
{
    eRegularGridOfSimilarDiamonds,
    eRegularGridOfDiamondsDifferentlySplit,
    eRegularGridofDiamondsWithHorizontalShifts,
    eStarcutSingleQuadrilateral,
    eDummy
};

// Calculates geometrical info for four adjacent quadrilaterals
// subdivided by diagonal segments to form triangular submesh.
Mesh generateSimilarDiamondsMesh(vector<double>& xc, vector<double>& yc, int splits)
{
    // at least 3 points
    assert(xc.size() > 2);
    assert(yc.size() > 2);
    // number of points should be even
    assert(xc.size() % 2 == 1);
    assert(yc.size() % 2 == 1);

    double x_split_inc = (xc[1]-xc[0])/(splits+1);

    Mesh mesh;

    int vertex_id  = 0;
    int segment_id = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    // prepare regular grid of vertices
    vector<vector<Vertex> > verts(yc.size());
    for (j = 0; j < yc.size(); j++)
    {
        verts[j].resize(xc.size());
        for (i = 0; i < xc.size(); i++)
        {
            verts[j][i].init (vertex_id++, xc[i], yc[j]);
            mesh.vert.push_back( verts[j][i] );
        }
    }

    // prepare carcass of edges
    vector<vector<Segment> > hseg(yc.size());
    vector<vector<Segment> > vseg(yc.size()-1);

    for (j = 0; j < verts.size(); j++)
    {
        hseg[j].resize(xc.size()-1);
    }
    // prepare only odd horizontal lines.
    // even lines are split below
    for (j = 0; j < verts.size(); j+=2)
    {
        for (i = 0; i < verts[j].size()-1; i++)
        {
            hseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+0][i+1]);
            mesh.seg.push_back( hseg[j][i] );

            if (j == 0)
            {
                mesh.south.push_back( hseg[j][i] );
            }
            if (j == verts[j].size()-1)
            {
                mesh.north.push_back( hseg[j][i] );
            }
        }
    }

    for (j = 0; j < verts.size()-1; j++)
    {
        vseg[j].resize(xc.size());
        for (i = 0; i < verts[j].size(); i++)
        {
            vseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+1][i+0]);
            mesh.seg.push_back( vseg[j][i] );

            if (i == 0)
            {
                mesh.west.push_back( vseg[j][i] );
            }
            if (i == verts[j].size()-1)
            {
                mesh.east.push_back( vseg[j][i] );
            }
        }
    }

    // loop through groups of 4 adjacent quadrilaterals
    // and define internal triangulation of these groups
    for (j = 0; j < yc.size()-1; j+=2)
    {
        for (i = 0; i < xc.size()-1; i+=2)
        {
            // vertices splitting even horizontal lines, including border vertices
            vector<Vertex> v_cli, v_cri;
            v_cli.push_back( verts[j+1][i+0] );
            v_cri.push_back( verts[j+1][i+2] );
            for (k = 0; k < splits; k++)
            {
                // these number from outside of the group towards the center
                v_cli.push_back( Vertex (vertex_id++, xc[i+0]+(k+1)*x_split_inc, yc[j+1]) );
                v_cri.push_back( Vertex (vertex_id++, xc[i+2]-(k+1)*x_split_inc, yc[j+1]) );
            }
            v_cli.push_back( verts[j+1][i+1] );
            v_cri.push_back( verts[j+1][i+1] );

            // saving vertices
            mesh.vert.insert( mesh.vert.end(), v_cli.begin()+1, v_cli.end()-1 );
            mesh.vert.insert( mesh.vert.end(), v_cri.begin()+1, v_cri.end()-1 );

            // define diagonal segments: upper left diagonal, ..
            vector<Segment> s_uld, s_urd, s_lld, s_lrd;
            // as well as horizontal segments splitting even lines
            vector<Segment> s_clh, s_crh;

            for (k = 0; k < splits+1; k++)
            {
                s_uld.push_back( Segment(segment_id++, verts[j+2][i+1], v_cli[k]) );
                s_lld.push_back( Segment(segment_id++, verts[j+0][i+1], v_cli[k]) );
                s_urd.push_back( Segment(segment_id++, verts[j+2][i+1], v_cri[k]) );
                s_lrd.push_back( Segment(segment_id++, verts[j+0][i+1], v_cri[k]) );
            }
            //s_clh.push_back( hseg[j+1][i+0] );
            //s_crh.push_back( hseg[j+1][i+1] );
            for (k = 0; k < splits+1; k++)
            {
                s_clh.push_back( Segment(segment_id++, v_cli[k], v_cli[k+1]) );
                s_crh.push_back( Segment(segment_id++, v_cri[k], v_cri[k+1]) );
            }

            // saving segments
            mesh.seg.insert( mesh.seg.end(), s_uld.begin(), s_uld.end() );
            mesh.seg.insert( mesh.seg.end(), s_lld.begin(), s_lld.end() );
            mesh.seg.insert( mesh.seg.end(), s_urd.begin(), s_urd.end() );
            mesh.seg.insert( mesh.seg.end(), s_lrd.begin(), s_lrd.end() );
            mesh.seg.insert( mesh.seg.end(), s_clh.begin(), s_clh.end() );
            mesh.seg.insert( mesh.seg.end(), s_crh.begin(), s_crh.end() );

            // corner triangles first
            mesh.tri.push_back( Triangle( hseg[j+2][i+0], vseg[j+1][i+0], s_uld[0]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+0], s_lld[0], vseg[j+0][i+0]) );
            mesh.tri.push_back( Triangle( hseg[j+2][i+1], s_urd[0], vseg[j+1][i+2]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+1], vseg[j+0][i+2], s_lrd[0]) );

            // internal triangles
            for (k = 0; k < splits; k++)
            {
                mesh.tri.push_back( Triangle( s_uld[k+0], s_clh[k+0], s_uld[k+1]) );
                mesh.tri.push_back( Triangle( s_lld[k+0], s_lld[k+1], s_clh[k+0]) );
                mesh.tri.push_back( Triangle( s_urd[k+0], s_urd[k+1], s_crh[k+0]) );
                mesh.tri.push_back( Triangle( s_lrd[k+1], s_lrd[k+0], s_crh[k+0]) );
            }

            // inner central triangles
            mesh.tri.push_back( Triangle( s_clh[splits], vseg[j+1][i+1], s_uld[splits]) );
            mesh.tri.push_back( Triangle( s_clh[splits], s_lld[splits], vseg[j+0][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[splits], s_urd[splits], vseg[j+1][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[splits], vseg[j+0][i+1], s_lrd[splits]) );
        }
    }
    return mesh;
}


// Calculates geometrical info for four adjacent quadrilaterals
// subdivided by diagonal segments to form triangular submesh.
//
// Generate regular grid of diamonds where within each group of
// 4 quadrilaterals their left quads have 'splits' number of
// horizontal splits (same as for regular grid above) while each
// right hand side pair of quadrilaterals is not split at all.
Mesh generateDiamondMeshDifferentlySplit(vector<double>& xc, vector<double>& yc, int splits)
{
    // at least 3 points
    assert(xc.size() > 2);
    assert(yc.size() > 2);
    // number of points should be even
    assert(xc.size() % 2 == 1);
    assert(yc.size() % 2 == 1);

    double x_split_inc = (xc[1]-xc[0])/(splits+1);

    Mesh mesh;

    int vertex_id  = 0;
    int segment_id = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    // prepare regular grid of vertices
    vector<vector<Vertex> > verts(yc.size());
    for (j = 0; j < yc.size(); j++)
    {
        verts[j].resize(xc.size());
        for (i = 0; i < xc.size(); i++)
        {
            verts[j][i].init (vertex_id++, xc[i], yc[j]);
            mesh.vert.push_back( verts[j][i] );
        }
    }

    // prepare carcass of edges
    vector<vector<Segment> > hseg(yc.size());
    vector<vector<Segment> > vseg(yc.size()-1);

    for (j = 0; j < verts.size(); j++)
    {
        hseg[j].resize(xc.size()-1);
    }
    // prepare only odd horizontal lines.
    // even lines are split below
    for (j = 0; j < verts.size(); j+=2)
    {
        for (i = 0; i < verts[j].size()-1; i++)
        {
            hseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+0][i+1]);
            mesh.seg.push_back( hseg[j][i] );

            if (j == 0)
            {
                mesh.south.push_back( hseg[j][i] );
            }
            if (j == verts[j].size()-1)
            {
                mesh.north.push_back( hseg[j][i] );
            }
        }
    }

    for (j = 0; j < verts.size()-1; j++)
    {
        vseg[j].resize(xc.size());
        for (i = 0; i < verts[j].size(); i++)
        {
            vseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+1][i+0]);
            mesh.seg.push_back( vseg[j][i] );

            if (i == 0)
            {
                mesh.west.push_back( vseg[j][i] );
            }
            if (i == verts[j].size()-1)
            {
                mesh.east.push_back( vseg[j][i] );
            }
        }
    }

    // loop through groups of 4 adjacent quadrilaterals
    // and define internal triangulation of these groups
    for (j = 0; j < yc.size()-1; j+=2)
    {
        for (i = 0; i < xc.size()-1; i+=2)
        {
            // make every even diamond in a horizontal direction
            // be simple diamond with only one diagonal edge splitting
            // each participating quad
            int this_diamond_splits = splits;
            if ((i % 4) >= 2)
            {
                this_diamond_splits = 0;
            }

            // vertices splitting even horizontal lines, including border vertices
            vector<Vertex> v_cli, v_cri;
            v_cli.push_back( verts[j+1][i+0] );
            v_cri.push_back( verts[j+1][i+2] );
            for (k = 0; k < this_diamond_splits; k++)
            {
                // these number from outside of the group towards the center
                v_cli.push_back( Vertex (vertex_id++, xc[i+0]+(k+1)*x_split_inc, yc[j+1]) );
                v_cri.push_back( Vertex (vertex_id++, xc[i+2]-(k+1)*x_split_inc, yc[j+1]) );
            }
            v_cli.push_back( verts[j+1][i+1] );
            v_cri.push_back( verts[j+1][i+1] );

            // saving vertices
            mesh.vert.insert( mesh.vert.end(), v_cli.begin()+1, v_cli.end()-1 );
            mesh.vert.insert( mesh.vert.end(), v_cri.begin()+1, v_cri.end()-1 );

            // define diagonal segments: upper left diagonal, ..
            vector<Segment> s_uld, s_urd, s_lld, s_lrd;
            // as well as horizontal segments splitting even lines
            vector<Segment> s_clh, s_crh;

            for (k = 0; k < this_diamond_splits+1; k++)
            {
                s_uld.push_back( Segment(segment_id++, verts[j+2][i+1], v_cli[k]) );
                s_lld.push_back( Segment(segment_id++, verts[j+0][i+1], v_cli[k]) );
                s_urd.push_back( Segment(segment_id++, verts[j+2][i+1], v_cri[k]) );
                s_lrd.push_back( Segment(segment_id++, verts[j+0][i+1], v_cri[k]) );
            }
            //s_clh.push_back( hseg[j+1][i+0] );
            //s_crh.push_back( hseg[j+1][i+1] );
            for (k = 0; k < this_diamond_splits+1; k++)
            {
                s_clh.push_back( Segment(segment_id++, v_cli[k], v_cli[k+1]) );
                s_crh.push_back( Segment(segment_id++, v_cri[k], v_cri[k+1]) );
            }

            // saving segments
            mesh.seg.insert( mesh.seg.end(), s_uld.begin(), s_uld.end() );
            mesh.seg.insert( mesh.seg.end(), s_lld.begin(), s_lld.end() );
            mesh.seg.insert( mesh.seg.end(), s_urd.begin(), s_urd.end() );
            mesh.seg.insert( mesh.seg.end(), s_lrd.begin(), s_lrd.end() );
            mesh.seg.insert( mesh.seg.end(), s_clh.begin(), s_clh.end() );
            mesh.seg.insert( mesh.seg.end(), s_crh.begin(), s_crh.end() );

            // corner triangles first
            mesh.tri.push_back( Triangle( hseg[j+2][i+0], vseg[j+1][i+0], s_uld[0]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+0], s_lld[0], vseg[j+0][i+0]) );
            mesh.tri.push_back( Triangle( hseg[j+2][i+1], s_urd[0], vseg[j+1][i+2]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+1], vseg[j+0][i+2], s_lrd[0]) );

            // internal triangles
            for (k = 0; k < this_diamond_splits; k++)
            {
                mesh.tri.push_back( Triangle( s_uld[k+0], s_clh[k+0], s_uld[k+1]) );
                mesh.tri.push_back( Triangle( s_lld[k+0], s_lld[k+1], s_clh[k+0]) );
                mesh.tri.push_back( Triangle( s_urd[k+0], s_urd[k+1], s_crh[k+0]) );
                mesh.tri.push_back( Triangle( s_lrd[k+1], s_lrd[k+0], s_crh[k+0]) );
            }

            // inner central triangles
            mesh.tri.push_back( Triangle( s_clh[this_diamond_splits], vseg[j+1][i+1], s_uld[this_diamond_splits]) );
            mesh.tri.push_back( Triangle( s_clh[this_diamond_splits], s_lld[this_diamond_splits], vseg[j+0][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[this_diamond_splits], s_urd[this_diamond_splits], vseg[j+1][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[this_diamond_splits], vseg[j+0][i+1], s_lrd[this_diamond_splits]) );
        }
    }
    return mesh;
}


// Calculates geometrical info for four adjacent quadrilaterals
// subdivided by diagonal segments to form triangular submesh.
Mesh generateDiamondMeshWithHorizontalShifts(vector<double>& xc, vector<double>& yc, int splits)
{
    // at least 3 points
    assert(xc.size() > 2);
    assert(yc.size() > 2);
    // number of points should be even
    assert(xc.size() % 2 == 1);
    assert(yc.size() % 2 == 1);

    double x_split_inc = (xc[1]-xc[0])/(splits+1);

    Mesh mesh;

    int vertex_id  = 0;
    int segment_id = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    // prepare regular grid of vertices
    vector<vector<Vertex> > verts(yc.size());
    for (j = 0; j < yc.size(); j++)
    {
        verts[j].resize(xc.size());
        for (i = 0; i < xc.size(); i++)
        {
            verts[j][i].init (vertex_id++, xc[i], yc[j]);
            mesh.vert.push_back( verts[j][i] );
        }
    }

    // prepare carcass of edges
    vector<vector<Segment> > hseg(yc.size());
    vector<vector<Segment> > vseg(yc.size()-1);

    for (j = 0; j < verts.size(); j++)
    {
        hseg[j].resize(xc.size()-1);
    }
    // prepare only odd horizontal lines.
    // even lines are split below
    for (j = 0; j < verts.size(); j+=2)
    {
        for (i = 0; i < verts[j].size()-1; i++)
        {
            hseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+0][i+1]);
            mesh.seg.push_back( hseg[j][i] );

            if (j == 0)
            {
                mesh.south.push_back( hseg[j][i] );
            }
            if (j == verts.size()-1)
            {
                mesh.north.push_back( hseg[j][i] );
            }
        }
    }

    for (j = 0; j < verts.size()-1; j++)
    {
        vseg[j].resize(xc.size());
        for (i = 0; i < verts[j].size(); i++)
        {
            vseg[j][i].init (segment_id++, verts[j+0][i+0], verts[j+1][i+0]);
            mesh.seg.push_back( vseg[j][i] );

            if (i == 0)
            {
                mesh.west.push_back( vseg[j][i] );
            }
            if (i == verts[j].size()-1)
            {
                mesh.east.push_back( vseg[j][i] );
            }
        }
    }

    // loop through groups of 4 adjacent quadrilaterals
    // and define internal triangulation of these groups
    for (j = 0; j < yc.size()-1; j+=2)
    {
        for (i = ((j % 4) >= 2); i < xc.size()-2; i+=2)
        {
            // vertices splitting even horizontal lines, including border vertices
            vector<Vertex> v_cli, v_cri;
            v_cli.push_back( verts[j+1][i+0] );
            v_cri.push_back( verts[j+1][i+2] );
            for (k = 0; k < splits; k++)
            {
                // these number from outside of the group towards the center
                v_cli.push_back( Vertex (vertex_id++, xc[i+0]+(k+1)*x_split_inc, yc[j+1]) );
                v_cri.push_back( Vertex (vertex_id++, xc[i+2]-(k+1)*x_split_inc, yc[j+1]) );
            }
            v_cli.push_back( verts[j+1][i+1] );
            v_cri.push_back( verts[j+1][i+1] );

            // saving vertices
            mesh.vert.insert( mesh.vert.end(), v_cli.begin()+1, v_cli.end()-1 );
            mesh.vert.insert( mesh.vert.end(), v_cri.begin()+1, v_cri.end()-1 );

            // define diagonal segments: upper left diagonal, ..
            vector<Segment> s_uld, s_urd, s_lld, s_lrd;
            // as well as horizontal segments splitting even lines
            vector<Segment> s_clh, s_crh;

            for (k = 0; k < splits+1; k++)
            {
                s_uld.push_back( Segment(segment_id++, verts[j+2][i+1], v_cli[k]) );
                s_lld.push_back( Segment(segment_id++, verts[j+0][i+1], v_cli[k]) );
                s_urd.push_back( Segment(segment_id++, verts[j+2][i+1], v_cri[k]) );
                s_lrd.push_back( Segment(segment_id++, verts[j+0][i+1], v_cri[k]) );
            }
            //s_clh.push_back( hseg[j+1][i+0] );
            //s_crh.push_back( hseg[j+1][i+1] );
            for (k = 0; k < splits+1; k++)
            {
                s_clh.push_back( Segment(segment_id++, v_cli[k], v_cli[k+1]) );
                s_crh.push_back( Segment(segment_id++, v_cri[k], v_cri[k+1]) );
            }

            // saving segments
            mesh.seg.insert( mesh.seg.end(), s_uld.begin(), s_uld.end() );
            mesh.seg.insert( mesh.seg.end(), s_lld.begin(), s_lld.end() );
            mesh.seg.insert( mesh.seg.end(), s_urd.begin(), s_urd.end() );
            mesh.seg.insert( mesh.seg.end(), s_lrd.begin(), s_lrd.end() );
            mesh.seg.insert( mesh.seg.end(), s_clh.begin(), s_clh.end() );
            mesh.seg.insert( mesh.seg.end(), s_crh.begin(), s_crh.end() );

            // corner triangles first
            mesh.tri.push_back( Triangle( hseg[j+2][i+0], vseg[j+1][i+0], s_uld[0]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+0], s_lld[0], vseg[j+0][i+0]) );
            mesh.tri.push_back( Triangle( hseg[j+2][i+1], s_urd[0], vseg[j+1][i+2]) );
            mesh.tri.push_back( Triangle( hseg[j+0][i+1], vseg[j+0][i+2], s_lrd[0]) );

            // internal triangles
            for (k = 0; k < splits; k++)
            {
                mesh.tri.push_back( Triangle( s_uld[k+0], s_clh[k+0], s_uld[k+1]) );
                mesh.tri.push_back( Triangle( s_lld[k+0], s_lld[k+1], s_clh[k+0]) );
                mesh.tri.push_back( Triangle( s_urd[k+0], s_urd[k+1], s_crh[k+0]) );
                mesh.tri.push_back( Triangle( s_lrd[k+1], s_lrd[k+0], s_crh[k+0]) );
            }

            // inner central triangles
            mesh.tri.push_back( Triangle( s_clh[splits], vseg[j+1][i+1], s_uld[splits]) );
            mesh.tri.push_back( Triangle( s_clh[splits], s_lld[splits], vseg[j+0][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[splits], s_urd[splits], vseg[j+1][i+1]) );
            mesh.tri.push_back( Triangle( s_crh[splits], vseg[j+0][i+1], s_lrd[splits]) );
        }

        // simplified triangulation for boundary groups of horizontally
        // shifted stripes
        if ((j % 4) >= 2)
        {

            int imin = 0;
            int imax = xc.size()-3;

            // horizontal segments
            Segment s_clh (segment_id++, verts[j+1][imin+0], verts[j+1][imin+1]);
            Segment s_crh (segment_id++, verts[j+1][imax+2], verts[j+1][imax+1]);

            // diagonal segments
            Segment s_uld(segment_id++, verts[j+0][imin+0], verts[j+1][imin+1] );
            Segment s_urd(segment_id++, verts[j+2][imin+0], verts[j+1][imin+1] );
            Segment s_lld(segment_id++, verts[j+0][imax+2], verts[j+1][imax+1] );
            Segment s_lrd(segment_id++, verts[j+2][imax+2], verts[j+1][imax+1] );


            // saving segments
            mesh.seg.push_back( s_clh );
            mesh.seg.push_back( s_crh );
            mesh.seg.push_back( s_uld );
            mesh.seg.push_back( s_urd );
            mesh.seg.push_back( s_lld );
            mesh.seg.push_back( s_lrd );

            // triangles
            mesh.tri.push_back( Triangle( hseg[j+0][imin+0], vseg[j+0][imin+1], s_uld) );
            mesh.tri.push_back( Triangle( s_clh, vseg[j+0][imin+0], s_uld) );

            mesh.tri.push_back( Triangle( s_clh, s_urd, vseg[j+1][imin+0]) );
            mesh.tri.push_back( Triangle( hseg[j+2][imin+0], s_urd, vseg[j+1][imin+1]) );

            mesh.tri.push_back( Triangle( hseg[j+0][imax+1], s_lld, vseg[j+0][imax+1]) );
            mesh.tri.push_back( Triangle( s_crh, s_lld, vseg[j+0][imax+2]) );

            mesh.tri.push_back( Triangle( s_crh, vseg[j+1][imax+2], s_lrd) );
            mesh.tri.push_back( Triangle( hseg[j+2][imax+1], vseg[j+1][imax+1], s_lrd) );
        }
    }
    return mesh;
}


// Calculates geometrical info for a single quadrilateral
// split into triangles in star-like manner
// to form triangular mesh.
Mesh generateStarcutSingleQuadMesh(int split_param)
{
    // at least 4 triangles
    assert(split_param > 0);

    // number of internal splits for each boundary edge of quad
    int splits = 0;
    if (split_param > 1)
    {
        splits = 2*(split_param-2)+1;
    }

    Mesh mesh;

    int vertex_id  = 0;
    int segment_id = 0;
    int j = 0;


    // prepare boundary vertices
    vector<Vertex> north_verts(splits+2);
    vector<Vertex> south_verts(splits+2);
    vector<Vertex> west_verts(splits+2);
    vector<Vertex> east_verts(splits+2);

    Vertex ll (vertex_id++, 0.0, 0.0);
    Vertex ul (vertex_id++, 0.0, 1.0);
    Vertex lr (vertex_id++, 1.0, 0.0);
    Vertex ur (vertex_id++, 1.0, 1.0);
    Vertex cc( vertex_id++, 0.5, 0.5);

    mesh.vert.push_back( ll );
    mesh.vert.push_back( ul );
    mesh.vert.push_back( ur );
    mesh.vert.push_back( lr );
    mesh.vert.push_back( cc );

    // prepare vertices inner to every boundary edge
    north_verts[0] = ul;
    south_verts[0] = ll;
    west_verts[0]  = ll;
    east_verts[0]  = lr;
    north_verts[splits+1] = ur;
    south_verts[splits+1] = lr;
    west_verts[splits+1]  = ul;
    east_verts[splits+1]  = ur;

    for (j = 1; j <= splits; j++)
    {
        north_verts[j].init (vertex_id++, (double)j * (1.0 / (splits+1)), 1.0);
        south_verts[j].init (vertex_id++, (double)j * (1.0 / (splits+1)), 0.0);
        west_verts[j].init  (vertex_id++, 0.0, (double)j * (1.0 / (splits+1)));
        east_verts[j].init  (vertex_id++, 1.0, (double)j * (1.0 / (splits+1)));

        mesh.vert.push_back( north_verts[j] );
        mesh.vert.push_back( south_verts[j] );
        mesh.vert.push_back(  west_verts[j] );
        mesh.vert.push_back(  east_verts[j] );
    }

    // prepare carcass of edges
    vector<Segment> north_seg(splits+1);
    vector<Segment> south_seg(splits+1);
    vector<Segment> west_seg(splits+1);
    vector<Segment> east_seg(splits+1);

    vector<Segment> north_diag_seg(splits+2);
    vector<Segment> south_diag_seg(splits+2);
    vector<Segment> west_diag_seg(splits+2);
    vector<Segment> east_diag_seg(splits+2);

    Segment s_ul (segment_id++, ul, cc);
    Segment s_ll (segment_id++, ll, cc);
    Segment s_lr (segment_id++, lr, cc);
    Segment s_ur (segment_id++, ur, cc);

    mesh.seg.push_back( s_ul );
    mesh.seg.push_back( s_ll );
    mesh.seg.push_back( s_lr );
    mesh.seg.push_back( s_ur );

    north_diag_seg[0] = s_ul;
    south_diag_seg[0] = s_ll;
     west_diag_seg[0] = s_ll;
     east_diag_seg[0] = s_lr;
    north_diag_seg[splits+1] = s_ur;
    south_diag_seg[splits+1] = s_lr;
     west_diag_seg[splits+1] = s_ul;
     east_diag_seg[splits+1] = s_ur;

    north_seg[0].init (segment_id++, north_verts[0], north_verts[1]);
    south_seg[0].init (segment_id++, south_verts[0], south_verts[1]);
     west_seg[0].init (segment_id++,  west_verts[0],  west_verts[1]);
     east_seg[0].init (segment_id++,  east_verts[0],  east_verts[1]);

    mesh.seg.push_back( north_seg[0] );
    mesh.seg.push_back( south_seg[0] );
    mesh.seg.push_back(  west_seg[0] );
    mesh.seg.push_back(  east_seg[0] );

    for (j = 1; j <= splits; j++)
    {
        north_seg[j].init (segment_id++, north_verts[j], north_verts[j+1]);
        south_seg[j].init (segment_id++, south_verts[j], south_verts[j+1]);
         west_seg[j].init (segment_id++,  west_verts[j],  west_verts[j+1]);
         east_seg[j].init (segment_id++,  east_verts[j],  east_verts[j+1]);

        north_diag_seg[j].init (segment_id++, north_verts[j], cc);
        south_diag_seg[j].init (segment_id++, south_verts[j], cc);
         west_diag_seg[j].init (segment_id++,  west_verts[j], cc);
         east_diag_seg[j].init (segment_id++,  east_verts[j], cc);

        mesh.seg.push_back( north_seg[j] );
        mesh.seg.push_back( south_seg[j] );
        mesh.seg.push_back(  west_seg[j] );
        mesh.seg.push_back(  east_seg[j] );

        mesh.seg.push_back( north_diag_seg[j] );
        mesh.seg.push_back( south_diag_seg[j] );
        mesh.seg.push_back(  west_diag_seg[j] );
        mesh.seg.push_back(  east_diag_seg[j] );
    }

    mesh.north.insert( mesh.north.end(), north_seg.begin(), north_seg.end() );
    mesh.south.insert( mesh.south.end(), south_seg.begin(), south_seg.end() );
    mesh.west.insert(  mesh.west.end(),   west_seg.begin(),  west_seg.end() );
    mesh.east.insert(  mesh.east.end(),   east_seg.begin(),  east_seg.end() );

    // define internal triangulation of a quad
    for (j = 0; j < splits+1; j++)
    {
        mesh.tri.push_back( Triangle( north_diag_seg[j], north_diag_seg[j+1], north_seg[j] ) );
        mesh.tri.push_back( Triangle( south_diag_seg[j], south_seg[j], south_diag_seg[j+1] ) );
        mesh.tri.push_back( Triangle(  west_diag_seg[j],  west_diag_seg[j+1],  west_seg[j] ) );
        mesh.tri.push_back( Triangle(  east_diag_seg[j],  east_seg[j],  east_diag_seg[j+1] ) );
    }
    return mesh;
}


int main(int argc, char *argv[])
{
    vector<double> xc, yc;
    int            nx = 0;
    int            ny = 0;
    double         sx = 1.0;
    double         sy = 1.0;
    int            nummodes = 7;
    int            splits   = 0;
    int            i;
    int            type;
    string         output_file;
    ofstream       output;

    if(argc != 9)
    {
        print_usage_info(argv[0]);
        exit(1);
    }

    try{
        type        = boost::lexical_cast<int>(argv[1]);
        splits      = boost::lexical_cast<int>(argv[2]);
        nx          = boost::lexical_cast<int>(argv[3]);
        ny          = boost::lexical_cast<int>(argv[4]);
        sx          = boost::lexical_cast<double>(argv[5]);
        sy          = boost::lexical_cast<double>(argv[6]);
        nummodes    = boost::lexical_cast<int>(argv[7]);
        output_file = boost::lexical_cast<string>(argv[8]);
        output.open(output_file.c_str());

        for (i = 0; i < nx; i++)
        {
            xc.push_back( (double)i * (1.0 / (nx-1)) );
        }
        for (i = 0; i < ny; i++)
        {
            yc.push_back( (double)i * (1.0 / (ny-1)) );
        }

        Mesh mesh;
        switch (type)
        {
        case eRegularGridOfSimilarDiamonds:
            mesh = generateSimilarDiamondsMesh(xc, yc, splits);
            break;
        case eRegularGridOfDiamondsDifferentlySplit:
            mesh = generateDiamondMeshDifferentlySplit(xc, yc, splits);
            break;
        case eRegularGridofDiamondsWithHorizontalShifts:
            mesh = generateDiamondMeshWithHorizontalShifts(xc, yc, splits);
            break;
        case eStarcutSingleQuadrilateral:
            mesh = generateStarcutSingleQuadMesh(splits);
            break;
        default:
            cout << "strange mesh type, aborting" << endl;
            throw 1;
        }

        output<< "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
        output<< "<NEKTAR>" << endl;
        output << "<EXPANSIONS>" << endl;
        output << "<E COMPOSITE=\"C[0]\" NUMMODES=\"" << nummodes << "\" FIELDS=\"u\" TYPE=\"MODIFIED\" />" <<endl;
        output << "</EXPANSIONS>\n" << endl;

        PrintConditions(output);

        output << "<GEOMETRY DIM=\"2\" SPACE=\"2\">" << endl;

        // -----------------------------------
        //  Vertices
        // -----------------------------------

        output << "<VERTEX XSCALE=\"" << sx << "\" YSCALE=\"" << sy << "\">" << endl;

        for (i = 0; i < mesh.vert.size(); i++)
        {
            output << "    <V ID=\"";
            output << mesh.vert[i].m_id << "\">\t";
            output << mesh.vert[i].m_x << " ";
            output << mesh.vert[i].m_y << " 0.0 </V>" << endl;
        }

        output << "  </VERTEX>\n" << endl;

        // -----------------------------------
        //  Edges
        // -----------------------------------

        output << "  <EDGE>" << endl;
        for(i = 0; i < mesh.seg.size(); i++)
        {
            output << "    <E ID=\"";
            output << mesh.seg[i].m_id << "\">\t";
            output << mesh.seg[i].m_v1.m_id << " ";
            output << mesh.seg[i].m_v2.m_id << " </E>" << endl;
        }
        output << "  </EDGE>\n" << endl;

        // -----------------------------------
        //  Elements
        // -----------------------------------

        output << "  <ELEMENT>" << endl;
        for(i = 0; i < mesh.tri.size(); i++)
        {
            output << "    <T ID=\"";
            output << i << "\">\t";
            output << mesh.tri[i].m_s1.m_id << " ";
            output << mesh.tri[i].m_s2.m_id << " ";
            output << mesh.tri[i].m_s3.m_id << " </T>" << endl;
        }
        output << "  </ELEMENT>\n" << endl;

        // -----------------------------------
        //  Composites
        // -----------------------------------

        output << "<COMPOSITE>" << endl;
        output << "<C ID=\"0\"> T[0-" << mesh.tri.size()-1 << "] </C>" << endl;

        // boundary composites coincide for both mesh element types
        output << "<C ID=\"1\"> E[";
        for(i = 0; i < mesh.south.size(); ++i)
        {
            output << mesh.south[i].m_id;
            if(i != mesh.south.size()-1) 
            {
                output << ",";
            }
        }
        output << "] </C>   // south border" << endl;

        output << "<C ID=\"2\"> E[";

        for(i = 0; i < mesh.west.size(); ++i)
        {
            output << mesh.west[i].m_id;
            if(i != mesh.west.size()-1) 
            {
                output << ",";
            }
        }
        output << "] </C>   // west border" << endl;

        output << "<C ID=\"3\"> E[";
        for(i = 0; i < mesh.north.size(); ++i)
        {
            output << mesh.north[i].m_id;
            if(i != mesh.north.size()-1) 
            {
                output << ",";
            }
        }
        output << "] </C>   // north border" << endl;

        output << "<C ID=\"4\"> E[";
        for(i = 0; i < mesh.east.size(); ++i)
        {
            output << mesh.east[i].m_id;
            if(i != mesh.east.size()-1) 
            {
                output << ",";
            }
        }
        output << "] </C>   // East border" << endl;

        output << "</COMPOSITE>\n" << endl;


        output << "<DOMAIN> C[0] </DOMAIN>\n" << endl;
        output << "</GEOMETRY>\n" << endl;
        output << "</NEKTAR>" << endl;

    }
    catch(...)
    {
        cerr << "Something went wrong. Caught an exception, stop." << endl;
        return 1;
    }

    return 0;

}



void  PrintConditions(ofstream& output)
{
    output << "<CONDITIONS>" << endl;

    output << "<SOLVERINFO>" << endl;
    output << "<I PROPERTY=\"GlobalSysSoln\"        VALUE=\"DirectFull\"/>" << endl;
    output << "</SOLVERINFO>\n" << endl;

    output << "<PARAMETERS>" << endl;
    output << "<P> TimeStep      = 0.002  </P>" << endl;
    output << "<P> Lambda        = 1      </P>" << endl;
    output << "</PARAMETERS>\n" << endl;

    output << "<VARIABLES>" << endl;
    output << "  <V ID=\"0\"> u </V>" << endl; 
    output << "</VARIABLES>\n" << endl;

    output << "<BOUNDARYREGIONS>" << endl;
    output << "<B ID=\"0\"> C[1] </B>" << endl;
    output << "<B ID=\"1\"> C[2] </B>" << endl;
    output << "<B ID=\"2\"> C[3] </B>" << endl;
    output << "<B ID=\"3\"> C[4] </B>" << endl;
    output << "</BOUNDARYREGIONS>\n" << endl;

    output << "<BOUNDARYCONDITIONS>" << endl;
    output << "  <REGION REF=\"0\"> // South border " << endl;
    output << "     <N VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;

    output << "  <REGION REF=\"1\"> // West border " << endl;
    output << "     <N VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;

    output << "  <REGION REF=\"2\"> // North border " << endl;
    output << "     <N VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;

    output << "  <REGION REF=\"3\"> // East border " << endl;
    output << "     <N VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;
    output << "</BOUNDARYCONDITIONS>" << endl;

    output << "  <FUNCTION NAME=\"Forcing\">" << endl;
    output << "     <E VAR=\"u\" VALUE=\"-(Lambda + 2*PI*PI/4)*sin(PI/2*x)*sin(PI/2*y)\" />" << endl;
    output << "  </FUNCTION>" << endl;

    output << "  <FUNCTION NAME=\"ExactSolution\">" << endl;
    output << "     <E VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />" << endl;
    output << "  </FUNCTION>" << endl;

    output << "</CONDITIONS>" << endl;
}

