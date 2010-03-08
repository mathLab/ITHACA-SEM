#include <algorithm>

#include <PreProcessing/Gmsh.h>

#define ERROR(msg)                              \
    string strmsg(msg);                         \
    cerr << strmsg << endl;                     \
    throw runtime_error(strmsg);

namespace Utilities
{
  namespace Gmsh
  {

    void ParseGmshFile(const char* inFile, const char* outfile)
    {
      string line;
      int nVertices        = 0;
      int nEntities        = 0;
      int nBoundComposites = 0;
      int expDim           = 3;
      int spaceDim   = 2;
      int elm_type = 0;

      vector<Vertex> vertices;
      vector<Edge> edges;
      vector<Face> faces;
      vector<ZeroDElement> zeroDElements;
      vector<OneDElement> oneDElements;
      vector<TwoDElement> twoDElements;
      vector<ThreeDElement> threeDElements;
      vector<Composite> composites;


      //---------------------------------------------
      // Read the *.msh file and fill the
      // nodes and geometric enteties structures

      fstream mshFile(inFile);

      if(!mshFile.is_open())
    {
      ERROR("Unable to find msh file");
    }
      else
    {
      cout <<"Start reading gmsh..." << endl;
      while(!mshFile.eof())
        {
          getline(mshFile, line);
          stringstream s(line);
          string word;
          s >> word;

          if(word == "$Nodes")
        {
          getline(mshFile, line);
          stringstream s(line);
          s >> nVertices;
          int id=0;
          for(int i=0; i<nVertices; ++i)
            {
              getline(mshFile, line);
              stringstream st(line);
              double x=0, y=0, z=0;
              st >> id >> x >> y >> z;

              if((z*z)>0.000001)
            {
              spaceDim = 3;
            }
              id -= 1; // counter starts at 0
              vertices.push_back( Vertex(id, x, y, z) );
            }
        }
          else if(word == "$Elements")
        {
          int zeroDid=0, oneDid=0, twoDid=0, threeDid=0;
          getline(mshFile, line);
          stringstream s(line);
          s >> nEntities;
          for(int i=0; i<nEntities; ++i)
            {
              getline(mshFile, line);
              stringstream st(line);
              int id=0, num_tag=0, num_nodes=0;

              st >> id >> elm_type >> num_tag;
              id -= 1; // counter starts at 0

              // curved edge
              if(elm_type==8)
            {
              num_nodes = 3;
            }

              // curved triangular or quadratic
              else if(elm_type==9)
            {
              num_nodes = 6;
            }

              else
            {
              num_nodes = GetNnodes(elm_type);
            }

              vector<int> tags;
              for(int j=0; j<num_tag; ++j)
            {
              int tag=0;
              st >> tag;
              tags.push_back(tag);

              // physical entities (used to define boundary entities)
              // are stored in tags[0]
              nBoundComposites = max(nBoundComposites,tags[0]);
            }

              vector<int> edgeList;
              vector<int> nodeList;
              vector<int> faceList;
              for(int k=0; k<num_nodes; ++k)
            {
              int node = 0;
              st >> node;
              node -= 1; // counter starts at 0
              nodeList.push_back(node);
            }

              switch(elm_type)
            {
            case 15:
              zeroDElements.push_back( ZeroDElement(zeroDid++,elm_type, tags, nodeList) );
              break;
            case 1:
            case 8:
              oneDElements.push_back( OneDElement(oneDid++,elm_type, tags, nodeList) );
              break;
            case 2:
            case 3:
            case 9:
              twoDElements.push_back( TwoDElement(twoDid++,elm_type, tags, nodeList, edgeList) );
              break;
            case 4:
            case 5:
              threeDElements.push_back( ThreeDElement(threeDid++, elm_type, tags, nodeList, faceList) );
              break;
            default:
              {
                ERROR("Gmsh element type not yet supported");
              }
            }
            }

          // find out what expansion dimension
          if (threeDElements.empty() && twoDElements.empty() && oneDElements.empty())
            {
               ERROR("Illegal expansion dimension");
             }
          else if (threeDElements.empty() && twoDElements.empty())
            expDim = 1;
          else if (threeDElements.empty())
            expDim = 2;
          else
            expDim = 3;


          cout << "Expansion dimension is " << expDim << endl;
          cout << "Read " << nVertices << " vertices" << endl;
          cout << "Read " << nEntities << " geometric entities" << endl;
          cout << "Read " << nBoundComposites << " boundary entities"  << endl;
        }

        }
      mshFile.close();
      //---------------------------------------------
    }

      switch(expDim)
    {
    case 1:
      SortZeroDElements(zeroDElements,vertices);
      SortOneDComposites(oneDElements, composites, nBoundComposites, expDim);
      break;
    case 2:
      // Fill the edges vectors from twoDElements.edge
      // When twoDElements are triangles, twoDElements.vert has three components.
      // With two vertices we find edge id and update it to edges.
      SortEdgeToVertex(twoDElements,edges);

      // Correct edges id of oneDElements from edges
      SortOneDElements(oneDElements,edges);

      // Set elements composite
      SortTwoDComposites(twoDElements, composites, nBoundComposites, expDim);
      SortOneDComposites(oneDElements, composites, nBoundComposites, expDim);
      break;
    case 3:
      {
        OrientTets(threeDElements, vertices);
        // Fill the faces vector from threeDElements.face
        SortFaceToVertex(threeDElements,faces);
        SortFaceToEdge(faces, edges);
        // Correct faces id of twoDElements from faces
        //SortTwoDElements(twoDElements,faces);
        SortThreeDComposites(threeDElements, composites, nBoundComposites, expDim);
        SortTwoDComposites(twoDElements, composites, nBoundComposites, expDim);
        SortOneDComposites(oneDElements, composites, nBoundComposites, expDim);
      }
      break;
    default:
      {
        ERROR("illegal expansion dimension");
      }
    }

      WriteToXMLFile(outfile, expDim, spaceDim, vertices, edges, faces, oneDElements, twoDElements, threeDElements, composites, elm_type);

    }


    // note: this function is only working for 2nodes edges.
    // Any of Gmsh higher order edges are not supported
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


    int GetFace(vector<int> &vert, vector<Face>& faces, int elm_type)
    {
        int i, j;
        int size = faces.size();
        
        // Check if face vertices match an existing face
        for(i = 0; i < size; i++)
        {
            bool match = true;
            for (j = 0; j < vert.size(); ++j)
            {
                vector<int>::iterator x = std::find(faces[i].vert.begin(), faces[i].vert.end(), vert[j]);
                if (x == faces[i].vert.end()) {
                    match = false;
                }
            }
            if (match)
            {
                return i;
            }
        }

        // for the last edge of an element
        vector<int> facevert;
        for (j = 0; j < vert.size(); ++j)
        {
            facevert.push_back(vert[j]);
        }
        faces.push_back( Face(size,facevert) );
        return size;
    }


    void SortEdgeToVertex(vector<TwoDElement> & elements, vector<Edge> & edges)
    {
      int i, j, elm_type;
      int size = elements.size();

      // fill the edge vector in the elements struct with unique edges
      for (i = 0; i < size; ++i)
    {
      elm_type = elements[i].type;

      if(elm_type==9)
        {
          for (j = 0; j < 2; ++j)
        {
          int edgeid;
          vector<int> vert;
          vert.push_back(elements[i].vert[j]);
          vert.push_back(elements[i].vert[j+1]);
          vert.push_back(elements[i].vert[j+3]);
          edgeid = GetEdge(vert,edges,elm_type);
          elements[i].edge.push_back(edgeid);
        }

          {
        int edgeid;
        vector<int> vert;
        vert.push_back(elements[i].vert[2]);
        vert.push_back(elements[i].vert[0]);
        vert.push_back(elements[i].vert[5]);
        edgeid = GetEdge(vert,edges,elm_type);
        elements[i].edge.push_back(edgeid);
          }

        }

      else
        {
          for (j = 0; j < elements[i].vert.size()-1; ++j)
        {
          int edgeid;
          vector<int> vert;
          vert.push_back(elements[i].vert[j]);
          vert.push_back(elements[i].vert[j+1]);
          edgeid = GetEdge(vert,edges,elm_type);
          elements[i].edge.push_back(edgeid);
        }

          {
        int edgeid;
        vector<int> vert;
        vert.push_back(elements[i].vert[elements[i].vert.size()-1]);
        vert.push_back(elements[i].vert[0]);
        edgeid = GetEdge(vert,edges,elm_type);
        elements[i].edge.push_back(edgeid);
          }
        }
    }
      cout << "...done sorting EdgeToVertex relations" << endl;
    }

    void OrientTets(vector<ThreeDElement> &elements, vector<Vertex> &vertices)
    {
        for (int i = 0; i < elements.size(); ++i)
        {
            // Don't do anything for Hex's
            if (elements[i].type != 4) continue;
            
            // Order vertices with lowest global vertex at top degenerate point
            // Place second lowest global vertex at base degenerate point
            vector<int> everts = elements[i].vert;
            sort(elements[i].vert.begin(), elements[i].vert.end());
            reverse(elements[i].vert.begin(), elements[i].vert.end());
            
            // Check orientation of tet and order remaining two points
            double ax, ay, az, vol;
            vector<Vertex> v;
            v.push_back(vertices[elements[i].vert[0]]);
            v.push_back(vertices[elements[i].vert[1]]);
            v.push_back(vertices[elements[i].vert[2]]);
            v.push_back(vertices[elements[i].vert[3]]);
            // Compute cross produc (b x c)
            ax = (v[1].y-v[3].y)*(v[2].z-v[3].z) - (v[1].z-v[3].z)*(v[2].y-v[3].y);
            ay = (v[2].x-v[3].x)*(v[1].z-v[3].z) - (v[1].x-v[3].x)*(v[2].z-v[3].z);
            az = (v[1].x-v[3].x)*(v[2].y-v[3].y) - (v[2].x-v[3].x)*(v[1].y-v[3].y);
            // Compute signed volume: 1/6 * (a . (b x c))
            vol = 1.0/6.0*(ax*(v[0].x-v[3].x) + ay*(v[0].y-v[3].y) + az*(v[0].z-v[3].z));
            // If negative volume, reverse order to correctly orientate tet.
            if (vol < 0)
            {
                swap(elements[i].vert[2], elements[i].vert[3]);
                swap(v[2], v[3]);
            }
        }
        cout << "...done orientating 3D elements" << endl;
    }
    
    void SortFaceToVertex(vector<ThreeDElement> & elements, vector<Face> & faces)
    {
        int i, j, k, elm_type;
        int size = elements.size();

        // fill the edge vector in the elements struct with unique edges
        for (i = 0; i < size; ++i)
        {
            elm_type = elements[i].type;
            // Hex
            if (elm_type == 5)
            {
                int face_ids[6][4] = {
                    {0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};
                for (j = 0; j < 6; ++j)
                {
                    int faceid;
                    vector<int> vert;
                    for (k = 0; k < 4; ++k)
                    {
                        vert.push_back(elements[i].vert[face_ids[j][k]]);
                    }
                    faceid = GetFace(vert,faces,elm_type);
                    elements[i].face.push_back(faceid);    
                }            
            }
            // Tet
            else if (elm_type == 4)
            {
                int face_ids[4][3] = {
                    {0,1,2},{0,1,3},{1,2,3},{2,0,3}};
                for (j = 0; j < 4; ++j)
                {
                    int faceid;
                    vector<int> vert;
                    for (k = 0; k < 3; ++k)
                    {
                        vert.push_back(elements[i].vert[face_ids[j][k]]);
                    }
                    faceid = GetFace(vert,faces,elm_type);
                    elements[i].face.push_back(faceid);
                }
            }
            else
            {
                ERROR("This element type is not supported yet.");
            }
        }
        cout << "...done sorting FaceToVertex relations" << endl;
    }
    
    void SortFaceToEdge(vector<Face> &faces, vector<Edge> &edges)
    {
        int i, j, elm_type;
        int size = faces.size();
        
        for (i = 0; i < size; ++i)
        {
            int fsize = faces[i].vert.size();
            if (fsize == 3) elm_type = 2;
            else if (fsize == 4) elm_type = 3;
            for (j = 0; j < fsize; ++j)
            {
                int edgeid;
                vector<int> vert;
                vert.push_back(faces[i].vert[j%fsize]);
                vert.push_back(faces[i].vert[(j+1)%fsize]);
                edgeid = GetEdge(vert,edges,elm_type);
                faces[i].edge.push_back(edgeid);
            }
        }
        cout << "...done sorting FaceToEdge relations" << endl;
    }
    
    
    // Get the correct edge id's from the edge struct
    void SortZeroDElements(vector<ZeroDElement> & points,const vector<Vertex> & vertices)
    {
      int i, j;
      int vertsize = vertices.size();
      int pointsize  = points.size();
      for (i = 0; i < pointsize; ++i)
    {
      for(j = 0; j < vertsize; j++)
        {
          if( points[i].vert[0] == vertices[j].x )
        {
          points[i].id = vertices[j].id;
          break;
        }
        }
    }
      cout << "...done sorting points ids" << endl;

    }



    // Get the correct edge id's from the edge struct
    void SortOneDElements(vector<OneDElement> &segments,const vector<Edge>& edges)
    {
      int i, j;
      int edgesize = edges.size();
      int segsize  = segments.size();
      for (i = 0; i < segsize; ++i)
    {
      for(j = 0; j < edgesize; j++)
        {
          if( ((segments[i].vert[0] == edges[j].vert[0])&&
           (segments[i].vert[1] == edges[j].vert[1]))
          ||
          ((segments[i].vert[1] == edges[j].vert[0])&&
           (segments[i].vert[0] == edges[j].vert[1])) )
        {
          segments[i].id = edges[j].id;
          break;
        }
        }
    }
      cout << "...done sorting segment ids" << endl;

    }

    void SortOneDComposites(const vector<OneDElement> & elements,
                vector<Composite> & composites, int num_composites, int dim)
    {
        int i,j;
        int nSeg  = 0;
        int oneDcomp = 0;

        // If our mesh is of dimension one, assemble the entire domain.
        if (dim == 1) {
            // set elements composite
            // for these composites we do not store
            // the individual element ids (as the composites
            // still will be defined from 0 to nElm-1), we just store the
            // total number of tris or quads
            for (int i = 0; i < elements.size(); ++i)
            {
                if ( (elements[i].type == 1))
                    nSeg++;
            }
    
            if (nSeg != 0)
            {
                list<int> eid;
                eid.push_back(nSeg-1); // -1 since start from 0
                composites.push_back(Composite(oneDcomp,1,eid));
                oneDcomp++;
            }
        }
        // Otherwise we have, for example, a boundary on a 2D mesh
        else {
            for (i = 0; i < elements.size(); ++i) {
                bool found = false;
                for (j = 0; j < composites.size(); ++j) {
                    if (composites[j].id == elements[i].tags[0]) {
                        composites[j].eid.push_back(elements[i].id);
                        found = true;
                    }
                }
                list<int> eid;
                eid.push_back(elements[i].id);
                if (!found) {
                    composites.push_back(Composite(elements[i].tags[0],1,eid));
                }
            }
        }
        cout << "...done sorting 1D composites" << endl;
    }


    void SortTwoDComposites(const vector<TwoDElement> & elements,
                vector<Composite> & composites, int num_composites, int dim)
    {
        int i,j;
        int nTri  = 0;
        int nQuad = 0;
        int twoDcomp = 0;

        if (dim == 2) {
            // set elements composite
            // for these composites we do not store
            // the individual element ids (as the composites
            // still will be defined from 0 to nElm-1), we just store the
            // total number of tris or quads
            for (int i = 0; i < elements.size(); ++i)
            {
                if ( (elements[i].type == 2) || (elements[i].type == 9) )
                    nTri++;
                if (elements[i].type == 3)
                    nQuad++;
            }
    
            if (nTri != 0)
            {
                list<int> eid;
                eid.push_back(nTri-1); // -1 since start from 0
                composites.push_back(Composite(twoDcomp,2,eid));
                twoDcomp++;
            }
            if (nQuad != 0)
            {
                list<int> eid;
                eid.push_back(nQuad-1); // -1 since start from 0
                composites.push_back(Composite(twoDcomp,3,eid));
                twoDcomp++;
            }
        }
        else {
            for (i = 0; i < elements.size(); ++i) {
                bool found = false;
                for (j = 0; j < composites.size(); ++j) {
                    if (composites[j].id == elements[i].tags[0]) {
                        composites[j].eid.push_back(elements[i].id);
                        found = true;
                    }
                }
                list<int> eid;
                eid.push_back(elements[i].id);
                if (!found) {
                    composites.push_back(Composite(elements[i].tags[0],2,eid));
                }
            }
        }
        cout << "...done sorting 2D composites" << endl;
    }


    void SortThreeDComposites(const vector<ThreeDElement> & elements,
                vector<Composite> & composites, int num_composites, int dim)
    {
        int i;

        // set elements composite
        // for these composites we do not store
        // the individual element ids (as the composites
        // still will be defined from 0 to nElm-1), we just store the
        // total number of tets or hexes
        int nTet  = 0;
        int nHex = 0;
        int threeDcomp = 0;
        for (int i = 0; i < elements.size(); ++i)
        {
            if ( (elements[i].type == 4) )
                nTet++;
            if (elements[i].type == 5)
                nHex++;
        }

        if (nTet != 0)
        {
            list<int> eid;
            eid.push_back(nTet-1); // -1 since start from 0
            composites.push_back(Composite(threeDcomp,4,eid));
            threeDcomp++;
        }
        if (nHex != 0)
        {
            list<int> eid;
            eid.push_back(nHex-1); // -1 since start from 0
            composites.push_back(Composite(threeDcomp,5,eid));
            threeDcomp++;
        }
        cout << "...done sorting 3D composites" << endl;
    }


    void WriteToXMLFile(const char* outfile, int expDim, int spaceDim, const vector<Vertex> & nodes, const vector<Edge> & edges,
            const vector<Face> & faces, const vector<OneDElement> & oneDElements,
            const vector<TwoDElement> & twoDElements, const vector<ThreeDElement> & threeDElements,
            const vector<Composite> & composites, const int elm_type)
    {

        TiXmlDocument doc;
        TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");
        doc.LinkEndChild( decl );

        TiXmlElement * root = new TiXmlElement( "NEKTAR" );
        doc.LinkEndChild( root );

        TiXmlComment * comment = new TiXmlComment();
        comment->SetValue(" Embed a 3-dimensional object in a 3-dimensional space");
        root->LinkEndChild( comment );
        comment = new TiXmlComment();
        comment->SetValue( "DIM <= SPACE  "
             "this provides a method of optimizing code for a 3-D curve embedded in 3-d space.");
        root->LinkEndChild( comment );

        //---------------------------------------------
        // Write DIM and SPACE

        TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
        geomTag->SetAttribute("DIM", expDim);
        geomTag->SetAttribute("SPACE", spaceDim);
        root->LinkEndChild( geomTag );

        comment = new TiXmlComment();
        comment->SetValue( " Definitions that can be used below in this file. ");
        geomTag->LinkEndChild( comment );
        //---------------------------------------------


        //---------------------------------------------
        // Write VERTEX

        TiXmlElement* verTag = new TiXmlElement( "VERTEX" );

        comment = new TiXmlComment();
        comment->SetValue( " Always must have four values per entry. ");
        verTag->LinkEndChild( comment );

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

        // only if expansion dimension >= 2
        if (expDim >= 2)
        {
            int edgecnt = 0;
            int small, large;
            verTag = new TiXmlElement( "EDGE" );

            comment = new TiXmlComment();
            comment->SetValue( "Edges are vertex pairs ");
            verTag->LinkEndChild( comment );

            for( int i = 0; i < edges.size(); ++i ) {
                stringstream s;

                s << setw(5) << edges[i].vert[0] << "  " << edges[i].vert[1] << "   ";
                TiXmlElement * e = new TiXmlElement( "E" );
                e->SetAttribute("ID",edgecnt++);
                e->LinkEndChild( new TiXmlText(s.str()) );
                verTag->LinkEndChild(e);
            }
            geomTag->LinkEndChild( verTag );
        }
        //---------------------------------------------


        //--------------------------------------------
        // Write FACES
        // only if expansion dimension == 3
        if (expDim == 3)
        {
            int facecnt = 0;
            int small, large;
            verTag = new TiXmlElement( "FACE" );

            comment = new TiXmlComment();
            comment->SetValue( "Faces are sets of edges ");
            verTag->LinkEndChild( comment );

            for( int i = 0; i < faces.size(); ++i ) {
                stringstream s;

                for ( int j = 0; j < faces[i].edge.size(); ++j)
                {
                    s << setw(5) << faces[i].edge[j];
                }
                TiXmlElement * f;
                if (faces[i].edge.size() == 3)
                {
                    f = new TiXmlElement( "T" );
                }
                else
                {
                    f = new TiXmlElement( "Q" );
                }
                f->SetAttribute("ID",facecnt++);
                f->LinkEndChild( new TiXmlText(s.str()) );
                verTag->LinkEndChild(f);
            }
            geomTag->LinkEndChild( verTag );
        }
        //---------------------------------------------


        //--------------------------------------------
        // Write ELEMENT

        verTag = new TiXmlElement( "ELEMENT" );

        comment = new TiXmlComment();
        comment->SetValue(  " P - points, Q -  quads, T - triangles, "
                            "S - segments, A - tetrahedron, Y - pyramid, "
                            "R - prism, H- hex.");
        verTag->LinkEndChild( comment );

        comment = new TiXmlComment();
        comment->SetValue(  "Only certain element types are appropriate for "
                            "the given dimension (dim on mesh)");
        verTag->LinkEndChild( comment );

        comment = new TiXmlComment();
        comment->SetValue(  "Can also use faces to define 3-D elements.  "
                            "Specify with F[1] for face 1, for example. ");
        verTag->LinkEndChild( comment );

        switch(expDim){
            // for 1D elements are Segments
            case 1:
            {
                for(int i=0; i<oneDElements.size(); ++i){
                    if ( (oneDElements[i].type == 1) 
                            || (oneDElements[i].type == 8) ){ // 1 denotes Seg

                        stringstream st;
                        for(int j=0; j<oneDElements[i].vert.size(); ++j){
                            st << setw(5) << oneDElements[i].vert[j] << " ";
                        }

                        TiXmlElement *elm_tag;
                        elm_tag = new TiXmlElement("S");
                        elm_tag->SetAttribute("ID", oneDElements[i].id);
                        elm_tag->LinkEndChild( new TiXmlText(st.str()) );
                        verTag->LinkEndChild(elm_tag);
                    }
                }
                break;
            }
            // for 2D elements are Tri or Quads
            case 2:
            {
                int elmcnt = 0;

                for(int i=0; i<twoDElements.size(); ++i){
                    if ( (twoDElements[i].type == 2) || (twoDElements[i].type == 3) ){ // 2 denotes Tri/3 denotes Quad

                        stringstream st;
                        for(int j=0; j<twoDElements[i].edge.size(); ++j){
                            st << setw(5) << twoDElements[i].edge[j] << " ";
                        }

                        TiXmlElement *elm_tag;
                        if (twoDElements[i].type == 2) {
                            elm_tag = new TiXmlElement("T");
                        }
                        else if (twoDElements[i].type == 3) {
                            elm_tag = new TiXmlElement("Q");
                        }
                        else{
                            cout << "No valid 2D element" << endl;
                            exit(1);
                        }

                        elm_tag->SetAttribute("ID", elmcnt++);
                        elm_tag->LinkEndChild( new TiXmlText(st.str()) );
                        verTag->LinkEndChild(elm_tag);
                    }
                    else if (twoDElements[i].type == 9){ // 9 denotes curved Tri
                        stringstream st;
                        for(int j=0; j<twoDElements[i].edge.size(); ++j){
                            st << setw(5) << twoDElements[i].edge[j] << " ";
                        }

                        TiXmlElement *elm_tag;
                        if (twoDElements[i].type == 9) {
                            elm_tag = new TiXmlElement("T");
                        }
                        else{
                            cout << "No valid 2D element" << endl;
                            exit(1);
                        }

                        elm_tag->SetAttribute("ID", elmcnt++);
                        elm_tag->LinkEndChild( new TiXmlText(st.str()) );
                        verTag->LinkEndChild(elm_tag);
                    }
                }
                break;
            }
            // for 3D elements are Tet, Hex, Prism or Pyr
            case 3:
            {
                int elmcnt = 0;

                for(int i=0; i<threeDElements.size(); ++i){
                    if (threeDElements[i].type == 4 || threeDElements[i].type == 5 ||
                            threeDElements[i].type == 6 || threeDElements[i].type == 7    ){

                        stringstream st;
                        for(int j=0; j< threeDElements[i].face.size(); ++j){
                            st << setw(5) << threeDElements[i].face[j] << " ";
                        }

                        TiXmlElement *elm_tag;
                        switch(threeDElements[i].type) {
                            case 4: elm_tag = new TiXmlElement("A"); break; // Tetrahedron
                            case 5: elm_tag = new TiXmlElement("H"); break; // Hexahedron
                            case 6: elm_tag = new TiXmlElement("R"); break; // Prism
                            case 7: elm_tag = new TiXmlElement("Y"); break; //Pyramid
                            default: cout << "No valid 3D element" << endl; exit(1);
                        }
                        elm_tag->SetAttribute("ID", elmcnt++);
                        elm_tag->LinkEndChild( new TiXmlText(st.str()) );
                        verTag->LinkEndChild(elm_tag);
                    }
                }
                break;
            }
            default:
            {
                ERROR("Dimension > 3");
            }
        }

        geomTag->LinkEndChild( verTag );
        //--------------------------------------------------

        if(elm_type==9)
        {
            // Write Curved elements
            TiXmlElement * curved = new TiXmlElement ("CURVED" );

            int edgecnt = 0;
            vector<int> order;

            order.push_back(0);
            order.push_back(2);
            order.push_back(1);

            for( int i = 0; i < edges.size(); ++i ) {
                stringstream s;

                for (int k = 0; k < edges[i].vert.size(); ++k) {
                    int nodeid = edges[i].vert[order[k]];
                    s << scientific << setprecision(3) << "     " 
                      <<  nodes[nodeid].x << "  " << nodes[nodeid].y 
                      << "  " << nodes[nodeid].z << "  ";
                }

                TiXmlElement * e = new TiXmlElement( "E" );
                e->SetAttribute("ID",edgecnt);
                e->SetAttribute("EDGEID",edgecnt++);
                e->SetAttribute("TYPE","PolyEvenlySpaced");
                e->SetAttribute("NUMPOINTS",3);
                TiXmlText * t0 = new TiXmlText(s.str());
                e->LinkEndChild(t0);
                curved->LinkEndChild(e);
            }
            geomTag->LinkEndChild( curved );
        }


        //--------------------------------------------------
        // Write COMPOSITE

        comment = new TiXmlComment();
        comment->SetValue( " V - vertex, S - segment, E - edge, F - face ");
        geomTag->LinkEndChild( comment );

        verTag = new TiXmlElement("COMPOSITE");

        list<int>::iterator it;
        list<int> compList;
        int id = 0;

        for(int i=0; i<composites.size(); ++i){
            stringstream st,st_start,st_end;
            compList = composites[i].eid;

            for (it = compList.begin(); it != compList.end(); ++it){
                st  << *it;
                if ( *it == compList.back())
                    st << "] ";
                else
                    st << ",";
            }

            if (compList.size() > 0) {
                TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
                TiXmlElement *elm_tag;

                switch(composites[i].type)
                {
                    case 1:
                    case 8:
                    {
                        if (expDim == 1)
                            st_start << " S[0-"; // Segment
                        else
                            st_start << " E[";  // Segment -> Edges
                        break;
                    }
                    case 2:
                    case 9:
                    {
                        if (expDim == 2)
                            st_start << " T[0-";  // Triangle
                        else
                            st_start << " F[";    // Triangle -> Face
                        break;
                    }
                    case 3:
                    {
                        if (expDim == 2)
                            st_start << " Q[0-";  // Quads
                        else
                            st_start << " F[";    // Quads -> Face
                        break;
                    }
                    case 4: st_start << " A[0-"; break;   // Tetrahedron
                    case 5: st_start << " H[0-"; break;   // Hexahedron
                    case 6: st_start << " R[0-"; break;   // Prism
                    case 7: st_start << " Y[0-"; break;   // Pyramid
                    case 15:st_start << " V["; break;     // Points -> Vertex
                }

                comp_tag->SetAttribute("ID", composites[i].id);
                comp_tag->LinkEndChild( new TiXmlText(st_start.str() 
                                                                + st.str()) );
                verTag->LinkEndChild(comp_tag);
            }
        }

        geomTag->LinkEndChild( verTag );
        //--------------------------------------------------

        TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
        domain->LinkEndChild( new TiXmlText( " C[0] " ));
        geomTag->LinkEndChild( domain );

        doc.SaveFile(outfile );
    } // end of function WriteToXMLFile




    int GetNnodes(int GmshEntity)
    {
      int nNodes;

      switch(GmshEntity)
    {
    case 1:  nNodes = 2;  break;
    case 2:  nNodes = 3;  break;
    case 3:  nNodes = 4;  break;
    case 4:  nNodes = 4;  break;
    case 5:  nNodes = 8;  break;
    case 6:  nNodes = 6;  break;
    case 7:  nNodes = 5;  break;
    case 8:  nNodes = 3;  break;
    case 9:  nNodes = 6;  break;
    case 10: nNodes = 9;  break;
    case 11: nNodes = 10; break;
    case 12: nNodes = 27; break;
    case 13: nNodes = 18; break;
    case 14: nNodes = 14; break;
    case 15: nNodes = 1;  break;
    case 16: nNodes = 8;  break;
    case 17: nNodes = 20; break;
    case 18: nNodes = 15; break;
    case 19: nNodes = 13; break;
    case 20: nNodes = 9;  break;
    case 21: nNodes = 10; break;
    case 22: nNodes = 12; break;
    case 23: nNodes = 15; break;
    case 24: nNodes = 15; break;
    case 25: nNodes = 21; break;
    case 26: nNodes = 4;  break;
    case 27: nNodes = 5;  break;
    case 28: nNodes = 6;  break;
    case 29: nNodes = 20; break;
    case 30: nNodes = 35; break;
    case 31: nNodes = 56; break;
    default:
      ERROR("unknown Gmsh element type");
    }

      return nNodes;
    }


  } //end of namespace Gmsh

} //end of namespace Utilities

