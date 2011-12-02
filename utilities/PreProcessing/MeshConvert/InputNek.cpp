////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek.cpp
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
//  Description: Nektar file format converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "InputNek.h"

#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputNek::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey("rea", eInputModule), InputNek::create);

        InputNek::InputNek(MeshSharedPtr m) : InputModule(m)
        {
            
        }

        InputNek::~InputNek()
        {
            
        }

        /**
         * @brief Processes Nektar file format.
         * 
         * Nektar sessions are defined by rea files, and contain sections
         * defining a DNS simulation in a specific order. The converter only
         * reads mesh information, curve information if it exists and boundary
         * information.
         * 
         * @param pFilename Filename of Nektar session file to read.
         */
        void InputNek::Process()
        {
            string      line, word;
            int         nParam, nModes, nElements, nCurves;
            int         i, j, k, nodeCounter = 0;
            ElementType elType;

            m->expDim   = 0;
            m->spaceDim = 0;
            
            cerr << "Start reading InputNek..." << endl;
            
            // -- Read in parameters.

            // Ignore first 3 lines. 4th line contains number of parameters.
            for (i = 0; i < 4; ++i)
            {
                getline(mshFile, line);
            }

            stringstream s(line);
            s >> nParam;
            
            for (i = 0; i < nParam; ++i)
            {
                string tmp1, tmp2;
                getline(mshFile, line);
                s.str(line);
                s >> tmp1 >> tmp2;
                if (tmp2 == "MODES")
                {
                    nModes = (int)boost::lexical_cast<double>(tmp1);
                }
            }
            
            // -- Read in passive scalars (ignore)
            getline(mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(mshFile, line);
            }

            // -- Read in logical switches (ignore)
            getline(mshFile, line);
            s.clear();
            s.str(line);
            s >> j;
            for (i = 0; i < j; ++i)
            {
                getline(mshFile, line);
            }

            // -- Read in mesh data.
            
            // First hunt for MESH tag
            bool foundMesh = false;
            while (!mshFile.eof())
            {
                getline(mshFile, line);
                if (line.find("MESH") != string::npos)
                {
                    foundMesh = true;
                    break;
                }
            }
            
            if (!foundMesh)
            {
                cerr << "Couldn't find MESH tag inside file." << endl;
                abort();
            }
            
            // Now read in number of elements and space dimension.
            getline(mshFile, line);
            s.clear(); s.str(line);
            s >> nElements >> m->expDim;
            m->spaceDim = m->expDim;
            
            // Loop over and create elements.
            for (i = 0; i < nElements; ++i)
            {
                getline(mshFile, line);
                
                if (m->expDim == 2)
                {
                    if (line.find("Qua") != string::npos || 
                        line.find("qua") != string::npos)
                    {
                        elType = eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is triangle.
                        elType = eTriangle;
                    }
                } 
                else
                {
                    if (line.find("Tet") != string::npos || 
                        line.find("tet") != string::npos) 
                    {
                        elType = eTetrahedron;
                    }
                    else if (line.find("Hex") != string::npos || 
                             line.find("hex") != string::npos) 
                    {
                        elType = eHexahedron;
                    }
                    else if (line.find("Prism") != string::npos || 
                             line.find("prism") != string::npos) 
                    {
                        elType = ePrism;
                    }
                    else if (line.find("Pyr") != string::npos || 
                             line.find("pyr") != string::npos) 
                    {
                        cerr << "Pyramid elements not yet supported." << endl;
                        abort();
                    }
                    else if (line.find("Qua") != string::npos || 
                             line.find("qua") != string::npos) 
                    {
                        elType = eQuadrilateral;
                    }
                    else
                    {
                        // Default element type in 2D is tetrahedron.
                        elType = eTetrahedron;
                    }
                }
                
                // Read in number of vertices for element type.
                int nNodes = GetNnodes(elType);
                double vertex[3][nNodes];
                
                for (j = 0; j < m->expDim; ++j)
                {
                    getline(mshFile,line);
                    s.clear(); s.str(line);
                    for (k = 0; k < nNodes; ++k)
                    {
                        s >> vertex[j][k];
                    }
                }
                
                for (j = m->expDim; j < 3; ++j)
                {
                    for (k = 0; k < nNodes; ++k)
                    {
                        vertex[j][k] = 0.0;
                    }
                }

                // Nektar meshes do not contain a unique list of nodes, so
                // this block constructs a unique set so that elements can be
                // created with unique nodes.
                vector<NodeSharedPtr> nodeList;
                for (k = 0; k < nNodes; ++k)
                {
                    NodeSharedPtr n = boost::shared_ptr<Node>(
                        new Node(nodeCounter++, vertex[0][k], 
                                 vertex[1][k],  vertex[2][k]));
                    pair<NodeSet::iterator, bool> testIns = m->vertexSet.insert(n);
                    
                    if (!testIns.second)
                    {
                        n = *(testIns.first);
                        nodeCounter--;
                    }
                    nodeList.push_back(n);
                }
                
                vector<int> tags;
                tags.push_back(elType);
                
                // Create linear element
                ElmtConfig conf(elType,1,false,false);
                ElementSharedPtr E = GetElementFactory().
                    CreateInstance(elType,conf,nodeList,tags);
                m->element[E->GetDim()].push_back(E);
            }

            // -- Process rest of mesh. This is done now to avoid overwriting
            // -- edges which have high-order information later.
            ProcessEdges     ();
            ProcessFaces     ();
            ProcessElements  ();
            ProcessComposites();

            // -- Read in curved data.
            getline(mshFile, line);
            if (line.find("CURVE") == string::npos)
            {
                cerr << "Cannot find curved side data." << endl;
                abort();
            }
            
            // Read number of curves.
            getline(mshFile, line);
            s.clear(); s.str(line);
            s >> nCurves;

            if (nCurves > 0)
            {
                string curveTag;
                
                for (i = 0; i < nCurves; ++i)
                {
                    getline(mshFile, line);
                    s.clear(); s.str(line);
                    s >> word;
                    
                    if (word == "File")
                    {
                        // Next line contains filename and curve tag.
                        getline(mshFile, line);
                        s.clear(); s.str(line);
                        s >> word >> curveTag;
                        curveTags[curveTag] = word;
                    }
                    else
                    {
                        cerr << "Unsupported curve type " << word << endl;
                        abort();
                    }
                }
                
                // Load high order surface information.
                LoadHOSurfaces();
                
                // Read in connectivity information. First line should contain
                // number of curved sides.
                getline(mshFile,line);
                
                if (line.find("side") == string::npos)
                {
                    cerr << "Unable to read number of curved sides" << endl;
                    abort();
                }
                
                int nCurvedSides;
                int faceId, elId, vid1, vid2, vid3;
                map<string,string>::iterator it;
                HOSurfSet::iterator hoIt;

                s.clear(); s.str(line);
                s >> nCurvedSides;

                // Iterate over curved sides, and look up high-order surface
                // information in the HOSurfSet, then map this onto faces.
                for (i = 0; i < nCurvedSides; ++i)
                {
                    getline(mshFile, line);
                    s.clear(); s.str(line);
                    s >> faceId >> elId >> word;
                    faceId--;
                    elId--;
                    
                    it = curveTags.find(word);
                    if (it == curveTags.end())
                    {
                        cerr << "Unrecognised curve tag " << word << " in curved lines" << endl;
                        abort();
                    }
                    
                    vector<unsigned int> vertId(3);
                    s >> vertId[0] >> vertId[1] >> vertId[2];

                    // Find vertex combination in hoData.
                    hoIt = hoData[word].find(HOSurfSharedPtr(new HOSurf(vertId)));
                    
                    if (hoIt == hoData[word].end())
                    {
                        cerr << "Unable to find high-order surface data for element id " 
                             << elId+1 << endl;
                        abort();
                    }

                    // Depending on order of vertices in rea file, surface
                    // information may need to be rotated or reflected. These
                    // procedures are taken from nektar/Hlib/src/HOSurf.C
                    HOSurfSharedPtr surf = *hoIt;

                    if (vertId[0] == surf->vertId[0]) 
                    {
                        if (vertId[1] == surf->vertId[1] || 
                            vertId[1] == surf->vertId[2])
                        {
                            if (vertId[1] == surf->vertId[2])
                            {
                                surf->Rotate(1);
                                surf->Reflect();
                            }
                        }
                    }
                    else if (vertId[0] == surf->vertId[1])
                    {
                        if (vertId[1] == surf->vertId[0] ||
                            vertId[1] == surf->vertId[2]) 
                        {
                            if (vertId[1] == surf->vertId[0])
                                surf->Reflect();
                            else
                                surf->Rotate(2);
                        }
                    }
                    else if (vertId[0] == surf->vertId[2])
                    {
                        if (vertId[1] == surf->vertId[0] ||
                            vertId[1] == surf->vertId[1])
                        {
                            if (vertId[1] == surf->vertId[1])
                            {
                                surf->Rotate(2);
                                surf->Reflect();
                            }
                            else
                            {
                                surf->Rotate(1);
                            }
                        }
                    }

                    // If the element is a prism, check to see if orientation
                    // has changed and update order of surface vertices.
                    ElementSharedPtr e = m->element[m->expDim][elId];
                    int reverseSide = 2;
                    
                    // Prisms may have been rotated by OrientPrism routine and
                    // break curved faces. This block rotates faces
                    // accordingly.
                    if (e->GetConf().e == ePrism)
                    {
                        boost::shared_ptr<Prism> pr = boost::static_pointer_cast<Prism>(e);
                        if (pr->orientation == 1)
                        {
                            // Prism has been rotated clockwise; rotate face,
                            // reverse what was the last edge (now located at
                            // edge 0).
                            (*hoIt)->Rotate(1);
                            reverseSide = 0;
                        }
                        else if (pr->orientation == 2)
                        {
                            // Prism has been rotated counter-clockwise;
                            // rotate face, reverse what was the last edge
                            // (now located at edge 1).
                            (*hoIt)->Rotate(2);
                            reverseSide = 1;
                        }
                    }
                    
                    // Finally, add high order data to appropriate face. NOTE:
                    // this is a bit of a hack since the elements are
                    // technically linear, but should work just fine.
                    FaceSharedPtr f    = e->GetFace(faceId);
                    EdgeSharedPtr edge;
                    int           Ntot = (*hoIt)->surfVerts.size();
                    int           N    = ((int)sqrt(8*Ntot+1)-1)/2;
                    
                    // Apply high-order map to convert face data to Nektar++
                    // ordering (vertices->edges->internal).
                    vector<NodeSharedPtr> tmpVerts = (*hoIt)->surfVerts;
                    for (j = 0; j < tmpVerts.size(); ++j)
                    {
                        (*hoIt)->surfVerts[hoMap[j]] = tmpVerts[j];
                    }
                    
                    for (j = 0; j < f->edgeList.size(); ++j)
                    {
                        edge = f->edgeList[j];
                        
                        // Skip over edges which have already been populated,
                        // apart from those which need to be reoriented.
                        if (edge->edgeNodes.size() > 0 && reverseSide == 2)
                            continue;
                        
                        edge->edgeNodes.clear();
                        edge->curveType = LibUtilities::eGaussLobattoLegendre;
                        
                        for (int k = 0; k < N-2; ++k)
                        {
                            edge->edgeNodes.push_back((*hoIt)->surfVerts[3+j*(N-2)+k]);
                        }
                        
                        // Reverse order of modes along correct side.
                        if (j == reverseSide)
                        {
                            reverse(edge->edgeNodes.begin(), edge->edgeNodes.end());
                        }
                    }
                }
            }

            mshFile.close();
        }

        /**
         * Load high order surface information from hsf file.
         */
        void InputNek::LoadHOSurfaces()
        {
            map<string,string>::iterator it;
            int nodeId = m->GetNumEntities();
            
            for (it = curveTags.begin(); it != curveTags.end(); ++it)
            {
                ifstream hsf;
                string   line, fileName = it->second;
                size_t   pos;
                int      N, Nface, dot;

                // Replace fro extension with hsf.
                dot = fileName.find_last_of('.');
                fileName = fileName.substr(0,dot);
                fileName += ".hsf";
                
                // Open hsf file.
                hsf.open(fileName.c_str());
                if (!hsf.is_open())
                {
                    cerr << "Could not open surface file " << fileName << endl;
                    abort();
                }
                
                // Read in header line; determine element order, number of faces
                // from this line.
                getline(hsf, line);
                pos = line.find("=");
                if (pos == string::npos)
                {
                    cerr << "hsf header error: cannot read number of nodal points." << endl;
                    abort();
                }
                line = line.substr(pos+1);
                stringstream ss(line);
                ss >> N;
                
                pos = line.find("=");
                if (pos == string::npos)
                {
                    cerr << "hsf header error: cannot read number of faces." << endl;
                    abort();
                }
                line = line.substr(pos+1);
                ss.clear(); ss.str(line);
                ss >> Nface;
                
                int Ntot = N*(N+1)/2;

                // Skip a line, then read in r,s positions inside the next
                // comments.
                Array<OneD, NekDouble> r(Ntot), s(Ntot);
                getline(hsf, line);
                
                for (int i = 0; i < 2; ++i)
                {
                    string word;
                    
                    getline(hsf, line);
                    ss.clear(); ss.str(line);
                    ss >> word;
                    
                    if (word != "#")
                    {
                        cerr << "hsf header error: cannot read in r/s points" << endl;
                        abort();
                    }
                    
                    for (int j = 0; j < Ntot; ++j)
                    {
                        ss >> (i == 0 ? r[j] : s[j]);
                    }
                }

                // Generate electrostatic points so that re-mapping array can
                // be constructed.
                Array<OneD, NekDouble> rp(Ntot), sp(Ntot);
                LibUtilities::PointsKey elec(N, LibUtilities::eNodalTriElec);
                LibUtilities::PointsManager()[elec]->GetPoints(rp,sp);

                // Expensively construct remapping array nodemap. This will
                // map nodal ordering from hsf order to Nektar++ ordering
                // (i.e. vertices followed by edges followed by interior
                // points.)
                for (int i = 0; i < Ntot; ++i)
                {
                    for (int j = 0; j < Ntot; ++j)
                    {
                        if (fabs(r[i]-rp[j]) < 1e-5 && fabs(s[i]-sp[j]) < 1e-5)
                        {
                            hoMap[i] = j;
                            break;
                        }
                    }
                }
                
                // Skip variables line
                getline(hsf, line);

                // Read in nodal points for each face.
                map<unsigned int, vector<NodeSharedPtr> > faceMap;
                for (int i = 0; i < Nface; ++i)
                {
                    getline(hsf, line);
                    vector<NodeSharedPtr> faceNodes(Ntot);
                    for (int j = 0; j < Ntot; ++j, ++nodeId)
                    {
                        double x, y, z;
                        getline(hsf, line);
                        ss.clear(); ss.str(line);
                        ss >> x >> y >> z;
                        faceNodes[j] = NodeSharedPtr(new Node(nodeId, x, y, z));
                    }
                    // Skip over tecplot connectivity information.
                    for (int j = 0; j < (N-1)*(N-1); ++j)
                    {
                        getline(hsf, line);
                    }
                    faceMap[i] = faceNodes;
                }
                
                // Finally, read in connectivity information to set up after
                // reading rea file.
                getline(hsf, line);
                for (int i = 0; i < Nface; ++i)
                {
                    string               tmp;
                    int                  fid;
                    vector<unsigned int> nodeIds(3);
                    
                    getline(hsf, line);
                    ss.clear(); ss.str(line);
                    ss >> tmp >> fid >> nodeIds[0] >> nodeIds[1] >> nodeIds[2];

                    if (tmp != "#")
                    {
                        cerr << "Unable to read hsf connectivity information." << endl;
                        abort();
                    }
                    
                    hoData[it->first].insert(HOSurfSharedPtr(new HOSurf(nodeIds, faceMap[i])));
                }
                
                hsf.close();
            }
        }
        
        /**
         * This routine aids the reading of Nektar files only; it returns the
         * number of nodes for a given entity typw.
         */
        int InputNek::GetNnodes(ElementType InputNekEntity)
        {
            int nNodes;

            switch(InputNekEntity)
            {
            case ePoint:         nNodes = 1;  break;
            case eLine:          nNodes = 2;  break;
            case eTriangle:      nNodes = 3;  break;
            case eQuadrilateral: nNodes = 4;  break;
            case eTetrahedron:   nNodes = 4;  break;
            case ePrism:         nNodes = 6;  break;
            case eHexahedron:    nNodes = 8;  break;
            default:
                cerr << "unknown Nektar element type" << endl;
            }

            return nNodes;
        }

        /** 
         * @brief Compares two %HOSurf objects referred to as shared pointers.
         *
         * Two %HOSurf objects are defined to be equal if they contain
         * identical vertex ids contained in HOSurf::vertId.
         */
        bool operator==(HOSurfSharedPtr const &p1, HOSurfSharedPtr const &p2)
        {
            if (p1->vertId.size() != p2->vertId.size())
            {
                return false;
            }
            
            vector<unsigned int> ids1 = p1->vertId;
            vector<unsigned int> ids2 = p2->vertId;
            sort(ids1.begin(), ids1.end());
            sort(ids2.begin(), ids2.end());
            
            for (int i = 0; i < ids1.size(); ++i)
            {
                if (ids1[i] != ids2[i])
                    return false;
            }
            
            return true;
        }
        
        /** 
         * Rotates the triangle of data points inside #surfVerts
         * counter-clockwise nrot times.
         *
         * @param nrot Number of times to rotate triangle.
         */
        void HOSurf::Rotate(int nrot)
        {
            int n, i, j, cnt;
            int np = ((int)sqrt(8*surfVerts.size()+1)-1)/2;
            NodeSharedPtr tmp[np][np];
            
            for (n = 0; n < nrot; ++n) 
            {
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j, cnt++)
                    {
                        tmp[i][j] = surfVerts[cnt];
                    }
                }
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j,cnt++)
                    {
                        surfVerts[cnt] = tmp[np-1-i-j][i];
                    }
                }
            }
        }
        
        void HOSurf::Reflect()
        {
            int n, i, j, cnt;
            int np = ((int)sqrt(8*surfVerts.size()+1)-1)/2;
            NodeSharedPtr tmp[np][np];
            
            for (cnt = i = 0; i < np; ++i)
            {
                for (j = 0; j < np-i; ++j,cnt++)
                {
                    tmp[i][np-i-1-j] = surfVerts[cnt];
                }
            }
            
            for(cnt = i = 0; i < np; ++i)
            {
                for(j = 0; j < np-i; ++j,cnt++)
                {
                    surfVerts[cnt] = tmp[i][j];
                }
            }
        }
    }
}
