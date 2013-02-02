#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <list>

#include <tinyxml/tinyxml.h>
using namespace std;


struct Vertex
{
  Vertex(int id, double x, double y, double z):
id(id),
x(x),
y(y),
z(z)
  {};
  int id;
  double x;
  double y;
  double z;
};

struct Edge
{
  Edge(int id, vector<int> vert):
id(id),
vert(vert)
  {};
  int id;
  vector<int> vert;
};

struct TwoDElement
{
  TwoDElement(int id, int type, vector<int> tags, vector<int> vert, vector<int> edge):
id(id),
type(type),
tags(tags),
verts(vert),
edges(edge)
  {};
  TwoDElement(int id, int type, vector<int> vert):
      id(id),
      type(type),
      verts(vert)
  {};
  int id;
  int type;
  vector<int> tags;
  vector<int> verts;
  vector<int> edges;
};

void readEnsite(const char* filename,
        vector<Vertex>& pVertexList,
        vector<TwoDElement>& pElementList)
{
    cout << "Reading Ensite file..." << flush;
    TiXmlDocument vEnsite(filename);
    TiXmlElement * vElmt;
    stringstream vTmp;
    int n;

    if (!vEnsite.LoadFile()) {
        cout << "Unable to open file " << filename << endl;
        exit(-1);
    }

    vElmt = vEnsite.FirstChildElement("DIF");
    if (!vElmt) cout << "Unable to find DIF tag" << endl;
    vElmt = vElmt->FirstChildElement("DIFBody");
    if (!vElmt) cout << "Unable to find DIFBody tag" << endl;
    vElmt = vElmt->FirstChildElement("Volumes");
    if (!vElmt) cout << "Unable to find Volumes tag" << endl;
    vElmt = vElmt->FirstChildElement("Volume");
    if (!vElmt) cout << "Unable to find Volume tag" << endl;

    TiXmlElement * vVertices = vElmt->FirstChildElement("Vertices");
    if (!vVertices) cout << "Unable to find Vertices tag" << endl;
    int vNVertices = atoi(vVertices->Attribute("number"));
    TiXmlElement * vPolygons = vElmt->FirstChildElement("Polygons");
    if (!vPolygons) cout << "Unable to find Polygons tag" << endl;
    int vNPolygons = atoi(vPolygons->Attribute("number"));

    // Read vertices
    vTmp.write(vVertices->GetText(), strlen(vVertices->GetText()));
    double x, y, z;
    for (n = 0; n < vNVertices; ++n)
    {
        vTmp >> x >> y >> z;
        pVertexList.push_back(Vertex(n,x,y,z));
    }

    // Read polygons
    vTmp.clear();
    vTmp.write(vPolygons->GetText(), strlen(vPolygons->GetText()));
    int a, b, c;
    for (n = 0; n < vNPolygons; ++n)
    {
        vTmp >> a >> b >> c;
        vector<int> vert;
        vert.push_back(a-1);
        vert.push_back(b-1);
        vert.push_back(c-1);
        pElementList.push_back(TwoDElement(n,2,vert));
    }
    cout << "done." << endl;
}

int GetEdge(vector<int> &vert, vector<Edge>& edges, int elm_type)
{
  int i;
  int size = edges.size();
  for(i = 0; i < size; i++)
{
  if( ((vert[0] == edges[i].vert[0])&&
       (vert[1] == edges[i].vert[1]))
                ||
      ((vert[1] == edges[i].vert[0])&&
       (vert[0] == edges[i].vert[1])) )
    {
      return i;
    }
}

  // for the last edge of an element
  vector<int> edgevert;
  edgevert.push_back(vert[0]);
  edgevert.push_back(vert[1]);
  if(elm_type==9)
{
  edgevert.push_back(vert[2]);
}
  edges.push_back( Edge(size,edgevert) );
  return size;
}

void createEdgeList(vector<Vertex>& pVertexList,
                    vector<Edge>& pEdgeList,
                    vector<TwoDElement>& pElementList)
{
    int i,j,edgeid;

    cout << "Creating edge list..." << flush;
    for (i = 0; i < pElementList.size(); ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            vector<int> verts;
            verts.push_back(pElementList[i].verts[j%3]);
            verts.push_back(pElementList[i].verts[(j+1)%3]);
            edgeid = GetEdge(verts, pEdgeList, 2);
            pElementList[i].edges.push_back(edgeid);
        }
    }
    cout << "done." << endl;
}

void WriteToXMLFile(const char* outfile, const vector<Vertex> & nodes, const vector<Edge> & edges,
        const vector<TwoDElement> & twoDElements)
{
    cout << "Writing XML file..." << flush;

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
    doc.LinkEndChild( decl );

    TiXmlElement * root = new TiXmlElement( "NEKTAR" );
    doc.LinkEndChild( root );

    //---------------------------------------------
    // Write DIM and SPACE
    TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
    geomTag->SetAttribute("DIM", "2");
    geomTag->SetAttribute("SPACE", "3");
    root->LinkEndChild( geomTag );
    //---------------------------------------------

    //---------------------------------------------
    // Write VERTEX
    TiXmlElement* verTag = new TiXmlElement( "VERTEX" );

    for( int i = 0; i < nodes.size(); ++i ) {
        stringstream s;
        s << scientific << setprecision(3) <<  nodes[i].x << " "
          << nodes[i].y << " " << nodes[i].z;
        TiXmlElement * v = new TiXmlElement( "V" );
        v->SetAttribute("ID",nodes[i].id);
        v->LinkEndChild( new TiXmlText(s.str()) );
        verTag->LinkEndChild(v);
    }
    geomTag->LinkEndChild( verTag );
    //---------------------------------------------


    //---------------------------------------------
    // Write EDGE
    int edgecnt = 0;
    verTag = new TiXmlElement( "EDGE" );

    for( int i = 0; i < edges.size(); ++i ) {
        stringstream s;

        s << setw(5) << edges[i].vert[0] << "  " << edges[i].vert[1] << "   ";
        TiXmlElement * e = new TiXmlElement( "E" );
        e->SetAttribute("ID",edgecnt++);
        e->LinkEndChild( new TiXmlText(s.str()) );
        verTag->LinkEndChild(e);
    }
    geomTag->LinkEndChild( verTag );
    //---------------------------------------------


    //--------------------------------------------
    // Write ELEMENT
    verTag = new TiXmlElement( "ELEMENT" );
    int elmcnt = 0;

    for(int i=0; i<twoDElements.size(); ++i){
        stringstream st;
        for(int j=0; j<twoDElements[i].edges.size(); ++j){
            st << setw(5) << twoDElements[i].edges[j] << " ";
        }

        TiXmlElement *elm_tag;
        elm_tag = new TiXmlElement("T");
        elm_tag->SetAttribute("ID", elmcnt++);
        elm_tag->LinkEndChild( new TiXmlText(st.str()) );
        verTag->LinkEndChild(elm_tag);
    }
    geomTag->LinkEndChild( verTag );
    //--------------------------------------------------


    //--------------------------------------------------
    // Write COMPOSITE
    verTag = new TiXmlElement("COMPOSITE");

    list<int>::iterator it;
    list<int> compList;

    stringstream st;
    st << " T[0-" << twoDElements.size()-1 << "]";

    TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
    comp_tag->SetAttribute("ID", "0");
    comp_tag->LinkEndChild( new TiXmlText(st.str()) );
    verTag->LinkEndChild(comp_tag);
    geomTag->LinkEndChild( verTag );
    //--------------------------------------------------

    TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
    domain->LinkEndChild( new TiXmlText( " C[0] " ));
    geomTag->LinkEndChild( domain );

    doc.SaveFile(outfile );

    cout << "done." << endl;
} // end of function WriteToXMLFile

void WriteMshFile(const char* outfile, const vector<Vertex> & nodes, const vector<Edge> & edges,
        const vector<TwoDElement> & twoDElements)
{
    cout << "Write GMSH file..." << flush;
    int n;
    ofstream out(outfile);

    out << "$MeshFormat" << endl;
    out << "2.1 0 8" << endl;
    out << "$EndMeshFormat" << endl;
    out << "$Nodes" << endl;
    out << nodes.size() << endl;
    for (n = 0; n < nodes.size(); ++n)
    {
        out << n+1 << " " << nodes[n].x << " " << nodes[n].y << " "
                << nodes[n].z << endl;
    }
    out << "$EndNodes" << endl;
    out << "$Elements" << endl;
    out << twoDElements.size() << endl;
    for (n = 0; n < twoDElements.size(); ++n)
    {
        out << n+1 << " 2 3 1 0 0 " << twoDElements[n].verts[0] + 1 << " "
                << twoDElements[n].verts[1] + 1
                << " " << twoDElements[n].verts[2] + 1 << endl;
    }
    out << "$EndElements" << endl;
    cout << "done." << endl;
}

int main(int argc, char *argv[])
{
    if(argc != 4)
    {
        fprintf(stderr,"Usage: EnsiteToXml ensite.xml session.xml mesh.msh\n");
        exit(1);
    }

    vector<Vertex> vVertexList;
    vector<Edge> vEdgeList;
    vector<TwoDElement> vElementList;

    readEnsite(argv[1], vVertexList, vElementList);
    createEdgeList(vVertexList, vEdgeList, vElementList);
    WriteToXMLFile(argv[2], vVertexList, vEdgeList, vElementList);
    WriteMshFile(argv[3], vVertexList, vEdgeList, vElementList);
}
