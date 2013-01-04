///////////////////////////////////////////////////////////////////////////////
//
// File MoveMeshToCriticalLayer.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: MoveMesh to critical layer
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <tinyxml/tinyxml.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2D.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

#include "ExtractCriticalLayerFunctions.h"

enum SolveType
{
    eSolveX,
    eSolveY,
    eSolveXY,
    eNoSolve
};


struct MoveVerts
{
    vector<NekDouble> kspring;
    vector<int> springVid;
    SolveType solve; 
};

void GetInterfaceVerts(const int compositeID, 
                       SpatialDomains::MeshGraphSharedPtr &mesh,
                       vector<int> &InterfaceVerts);

void GetStreakLocation(LibUtilities::SessionReaderSharedPtr &vSession, 
                       SpatialDomains::MeshGraphSharedPtr &mesh, string &fieldfile,
                       Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc);

void GetNewVertexLocation(TiXmlElement *doc,
                          SpatialDomains::MeshGraphSharedPtr &mesh,
                          vector<int> &InterfaceVerts,
                          Array<OneD, NekDouble> &xstreak,
                          Array<OneD, NekDouble> &ystreak,
                          Array<OneD,NekDouble> &vertx,
                          Array<OneD,NekDouble> &verty,
                          int maxiter);

void  TurnOffEdges(TiXmlElement *doc, 
                   SpatialDomains::SegGeomMap &meshedges, 
                   Array<OneD,MoveVerts> &verts);

void   RedefineVertices(TiXmlElement *doc, 
                        Array<OneD,NekDouble> &dvertx, Array<OneD,NekDouble> &dverty);

void EnforceRotationalSymmetry(SpatialDomains::MeshGraphSharedPtr &mesh,
                               Array<OneD, NekDouble> &dvertx,
                               Array<OneD, NekDouble> &dverty);

int main(int argc, char *argv[])
{
    int i,j;
    NekDouble cr = 0;
    
    if(argc !=4)
    {
        fprintf(stderr,"Usage: ./MoveMeshToCriticalLayer  meshfile streakfile  outfile\n");
        exit(1);
    }
    
    //------------------------------------------------------------
    // Create Session file by reading only first meshfile
    //-----------------------------------------------------------
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(2, argv);
    
    //-------------------------------------------------------------
    // Read in mesh from input file
    //------------------------------------------------------------
    string meshfile(argv[argc-3]);
    SpatialDomains::MeshGraphSharedPtr mesh = SpatialDomains::MeshGraph::Read(vSession);
    TiXmlDocument& meshdoc = vSession->GetDocument();
    TiXmlHandle docHandle(&meshdoc);
    TiXmlElement* doc = docHandle.FirstChildElement("NEKTAR").FirstChildElement("GEOMETRY").Element();
    
    //------------------------------------------------------------
    // Get list of vertices along interface
    //------------------------------------------------------------
    vector<int> InterfaceVerts;
    GetInterfaceVerts(300,mesh,InterfaceVerts);

    //------------------------------------------------------------
    // Determine location of new interface vertices
    //------------------------------------------------------------
    Array<OneD,NekDouble> xstreak(50),ystreak(50);
    string fieldfile(argv[argc-2]);
    GetStreakLocation(vSession,mesh,fieldfile,xstreak,ystreak);

    //------------------------------------------------------------
    // Move internal mesh using critical layer info and under string analogy 
    //------------------------------------------------------------
    Array<OneD,NekDouble>  dvertx(mesh->GetNvertices(),0.0), dverty(mesh->GetNvertices(),0.0); 
    int maxiter;
    vSession->LoadParameter("MoveMeshMaxIterations",maxiter,100);

    GetNewVertexLocation(doc, mesh,InterfaceVerts,xstreak,ystreak,dvertx,dverty,maxiter);

    //------------------------------------------------------------
    // Enforce rotational symmetry on mesh 
    //------------------------------------------------------------
    if(vSession->DefinesSolverInfo("EnforceRotationalSymmetry"))
    {
        EnforceRotationalSymmetry(mesh,dvertx,dverty);
    }

    //------------------------------------------------------------
    // Redfine vertices in doc 
    //------------------------------------------------------------
    RedefineVertices(doc,dvertx,dverty);

    //------------------------------------------------------------
    // Write out moved mesh file 
    //------------------------------------------------------------
    std::string outfile(argv[argc-1]);
    meshdoc.SaveFile(outfile);
}


void GetInterfaceVerts(const int compositeID, SpatialDomains::MeshGraphSharedPtr &mesh, vector<int> &InterfaceVerts)
{
    SpatialDomains::Composite composite;
    composite = mesh->GetComposite(compositeID);
    int compsize = composite->size();
    map<int,int> vertmap;

    for(int i = 0; i < compsize; ++i)
    {
        if(vertmap.count((*composite)[i]->GetVid(0)) == 0)
        {
            InterfaceVerts.push_back((*composite)[i]->GetVid(0)); 
            vertmap[(*composite)[i]->GetVid(0)]  = 1;
        }

        if(vertmap.count((*composite)[i]->GetVid(1)) == 0)
        {
            InterfaceVerts.push_back((*composite)[i]->GetVid(1)); 
            vertmap[(*composite)[i]->GetVid(1)]  = 1;
        }
    }
}

void GetStreakLocation(LibUtilities::SessionReaderSharedPtr &vSession, 
                       SpatialDomains::MeshGraphSharedPtr &mesh, string &fieldfile,
                       Array<OneD, NekDouble> &xc, Array<OneD, NekDouble> &yc)
{
    //-------------------------------------------------------------
    // Define Streak Expansion   
    MultiRegions::ExpListSharedPtr streak;   
    
    streak = MemoryManager<MultiRegions::ExpList2D>
        ::AllocateSharedPtr(vSession,mesh);
    //---------------------------------------------------------------

    //----------------------------------------------
    // Import field file.
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    mesh->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    string  streak_field("w");
    for(unsigned int i = 0; i < fielddata.size(); ++i)
    {
        streak->ExtractDataToCoeffs(fielddef [i],
                                    fielddata[i],
                                    streak_field,
                                    streak->UpdateCoeffs());
    }
    //----------------------------------------------
    
    NekDouble cr = 0.0;

    cerr << "Extracting Critical Layer ";
    Computestreakpositions(streak,xc, yc,cr);    
}


void GetNewVertexLocation(TiXmlElement *doc, 
                          SpatialDomains::MeshGraphSharedPtr &mesh,
                          vector<int> &InterfaceVerts,
                          Array<OneD, NekDouble> &xstreak,
                          Array<OneD, NekDouble> &ystreak,
                          Array<OneD,NekDouble> &dvertx,
                          Array<OneD,NekDouble> &dverty,
                          int maxiter)
{
    int i,j,k;
    int nverts = mesh->GetNvertices(); 

    SpatialDomains::SegGeomMap meshedges = mesh->GetAllSegGeoms();

    Array<OneD,MoveVerts> Verts(nverts);
    
    // loop mesh edges and fill in verts info
    std::map<int,SpatialDomains::SegGeomSharedPtr>::iterator segIter;
    
    SpatialDomains::VertexComponentSharedPtr v0,v1;
    SpatialDomains::VertexComponent dist;

    int vid0,vid1;
    NekDouble kspring;
    NekDouble x,y,x1,y1,z1,x2,y2,z2;

    // Setup intiial spring and verts
    for(segIter = meshedges.begin(); segIter != meshedges.end(); ++segIter)
    {
        vid0 = (segIter->second)->GetVid(0);
        vid1 = (segIter->second)->GetVid(1);

        v0 = (segIter->second)->GetVertex(0);
        v1 = (segIter->second)->GetVertex(1);
        
        kspring = 1.0/v0->dist(*v1); 
        
        Verts[vid0].kspring.push_back(kspring);
        Verts[vid0].springVid.push_back(vid1);
        Verts[vid1].kspring.push_back(kspring);
        Verts[vid1].springVid.push_back(vid0);

    }
    
    // Scale spring by total around vertex and turn on all vertices. 
    for(i = 0; i < nverts; ++i)
    {
        NekDouble invktot = 0.0;
        for(j = 0; j < Verts[i].kspring.size(); ++j)
        {
            invktot += Verts[i].kspring[j];
        }
        invktot = 1.0/invktot;
        for(j = 0; j < Verts[i].kspring.size(); ++j)
        {
            Verts[i].kspring[j] *= invktot;
        }

        Verts[i].solve = eSolveXY;
    }        


    // Turn off all edges defined by composite lists of correct dimension
    TurnOffEdges(doc,meshedges,Verts);
    
    NekDouble z,h0,h1,h2;
    // Set interface vertices to lie on critical layer
    for(i = 0; i < InterfaceVerts.size(); ++i)
    {
        Verts[InterfaceVerts[i]].solve = eNoSolve; 
        mesh->GetVertex(InterfaceVerts[i])->GetCoords(x,y,z);
        
        for(j = 0; j < xstreak.num_elements()-1; ++j)
        {
            if((x >= xstreak[j])&&(x <= xstreak[j+1]))
            {
                break;
            }
        }
        
        ASSERTL0(j != xstreak.num_elements(),"Did not find x location along critical layer");
        
        k = (j==0)?1: j; // offset index at beginning
        
        // quadraticalling interpolate points
        h0 = (x-xstreak[k])*(x-xstreak[k+1])/
            ((xstreak[k-1]-xstreak[k])*(xstreak[k-1]-xstreak[k+1]));
        h1 = (x-xstreak[k-1])*(x-xstreak[k+1])/
            ((xstreak[k]-xstreak[k-1])*(xstreak[k]-xstreak[k+1]));
        h2 = (x-xstreak[k-1])*(x-xstreak[k])/
            ((xstreak[k+1]-xstreak[k-1])*(xstreak[k+1]-xstreak[k]));
        
        dvertx[InterfaceVerts[i]] =  (xstreak[k-1]*h0 + xstreak[k]*h1 + xstreak[k+1]*h2) - x; 
        dverty[InterfaceVerts[i]] =  (ystreak[k-1]*h0 + ystreak[k]*h1 + ystreak[k+1]*h2) - y; 
    }

    // shift quads in critical layer to move more or less rigidly
    SpatialDomains::QuadGeomMap quadgeom = mesh->GetAllQuadGeoms();
    std::map<int,SpatialDomains::QuadGeomSharedPtr>::iterator quadIter;
    for(quadIter = quadgeom.begin(); quadIter != quadgeom.end(); ++quadIter)
    {
        for(i = 0; i < 4; ++i)
        {
            vid0 = (quadIter->second)->GetVid(i);
            
            switch(Verts[vid0].solve)
            {
            case eSolveXY:
                {
                    mesh->GetVertex(vid0)->GetCoords(x,y,z);

                    // find nearest interface vert
                    mesh->GetVertex(InterfaceVerts[0])->GetCoords(x1,y1,z1);
                    for(j = 0; j < InterfaceVerts.size()-1; ++j)
                    {
                        mesh->GetVertex(InterfaceVerts[j+1])->GetCoords(x2,y2,z2);
                        if((x >= x1)&&(x < x2))
                        {
                            break;
                        }
                        x1 = x2;
                        y1 = y2;
                    }
                    
                    // currently just shift vert as average of two sides
                    dvertx[vid0] = (x2-x)/(x2-x1)*dvertx[InterfaceVerts[j]]+
                        (x-x1)/(x2-x1)*dvertx[InterfaceVerts[j+1]];
                    dverty[vid0] = (x2-x)/(x2-x1)*dverty[InterfaceVerts[j]]+
                        (x-x1)/(x2-x1)*dverty[InterfaceVerts[j+1]];
                }
                break;
            case eSolveY:
                {
                    mesh->GetVertex(vid0)->GetCoords(x,y,z);
                    mesh->GetVertex(InterfaceVerts[0])->GetCoords(x1,y1,z1);
                    
                    if(fabs(x-x1) < 1e-6)
                    {
                        dverty[vid0] = dverty[InterfaceVerts[0]];
                    }
                    else
                    {
                        dverty[vid0] = dverty[InterfaceVerts[InterfaceVerts.size()-1]];
                    }
                }
                break;
            default:
                break;
            }
            Verts[vid0].solve = eNoSolve;
        }
    }

    


            
    // Iterate internal vertices 
    bool ContinueToIterate = true;
    int cnt  = 0;
    int nsum = 0;
    NekDouble dsum,dx,dy,sum,prev_sum = 0.0;
    NekDouble tol = 1e-3,fac;
    int blend = 50;

    while (ContinueToIterate)
    {
        
        sum = 0.0;
        nsum = 0;

        // use a ramping function to help move interior slowly 
        fac = (cnt < blend)? 1.0/(blend+1.0)*(cnt+1): 1.0;
        
        for(i = 0; i < nverts; ++i)
        {
            if(Verts[i].solve != eNoSolve)
            {
                dx = dy = 0.0;

                if((Verts[i].solve == eSolveX)||(Verts[i].solve == eSolveXY))
                {
                    for(j = 0; j < Verts[i].kspring.size(); ++j)
                    {
                        dx += fac*Verts[i].kspring[j]*dvertx[Verts[i].springVid[j]];
                    }
                }

                if((Verts[i].solve == eSolveY)||(Verts[i].solve == eSolveXY))
                {
                    for(j = 0; j < Verts[i].kspring.size(); ++j)
                    {
                        dy += fac*Verts[i].kspring[j]*dverty[Verts[i].springVid[j]];
                    }
                }
                
                dsum = (dx*dx + dy*dy);
                
                dvertx[i] = dx;
                dverty[i] = dy;
                
                if(dsum > 1e-16)
                {
                    //sum  += dsum/(dvertx[i]*dvertx[i] + dverty[i]*dverty[i]);
                    sum  += dsum;
                    nsum += 1;
                }
            }
        }

        if(nsum)
        {
            sum = sqrt(sum/(NekDouble)nsum);
            
            NekDouble chg = sum-prev_sum;
            prev_sum = sum;

            cerr << "Iteration " << cnt << " : " << chg << endl;

            if((chg < tol)&&(cnt > blend))
            {
                ContinueToIterate = false;
            }
            
        }
        else if(cnt > blend)
        {
            ContinueToIterate = false;

        }            

        if(cnt++ > maxiter)
        {
            ContinueToIterate = false;
        }
    }
}


// Read Composites from xml document and turn off verts that are along edge composites. 
void  TurnOffEdges(TiXmlElement *doc, 
                   SpatialDomains::SegGeomMap &meshedges, 
                   Array<OneD,MoveVerts> &Verts)
{
    TiXmlElement* field = doc->FirstChildElement("COMPOSITE");
    ASSERTL0(field, "Unable to find COMPOSITE tag in file.");

    int nextCompositeNumber = -1;

    /// All elements are of the form: "<C ID = "N"> ... </C>".
    /// Read the ID field first.
    TiXmlElement *composite = field->FirstChildElement("C");

    while (composite)
    {
        nextCompositeNumber++;
        
        int indx;
        int err = composite->QueryIntAttribute("ID", &indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
        
        TiXmlNode* compositeChild = composite->FirstChild();
        // This is primarily to skip comments that may be present.
        // Comments appear as nodes just like elements.
        // We are specifically looking for text in the body
        // of the definition.
        while(compositeChild && compositeChild->Type() != TiXmlNode::TEXT)
        {
            compositeChild = compositeChild->NextSibling();
        }
        
        ASSERTL0(compositeChild, "Unable to read composite definition body.");
        std::string compositeStr = compositeChild->ToText()->ValueStr();
        
        /// Parse out the element components corresponding to type of element.
        std::istringstream compositeDataStrm(compositeStr.c_str());

        try
        {
            std::string compositeElementStr;
            compositeDataStrm >> compositeElementStr;
            
            std::istringstream tokenStream(compositeElementStr);
            char type;
            
            tokenStream >> type;
            
            // in what follows we are assuming there is only one block of data 
            std::string::size_type indxBeg = compositeElementStr.find_first_of('[') + 1;
            std::string::size_type indxEnd = compositeElementStr.find_last_of(']') - 1;
        
            ASSERTL0(indxBeg <= indxEnd, (std::string("Error reading index definition:") +  compositeElementStr).c_str());
            
            std::string indxStr = compositeElementStr.substr(indxBeg, indxEnd - indxBeg + 1);
            std::vector<unsigned int> seqVector;
            std::vector<unsigned int>::iterator seqIter;
            
            bool err = ParseUtils::GenerateSeqVector(indxStr.c_str(), seqVector);
            
            ASSERTL0(err, (std::string("Error reading composite elements: ") + indxStr).c_str());
            
            switch(type)
            {
            case 'E':  // Turn off vertices along composite edges 
                {
                    int seqlen = seqVector.size();
                    NekDouble x0,y0,z0,x1,y1,z1;
                    int vid0,vid1;

                    for(int i = 0; i < seqlen; ++i)
                    {
                        meshedges[seqVector[i]]->GetVertex(0)->GetCoords(x0,y0,z0);
                        meshedges[seqVector[i]]->GetVertex(1)->GetCoords(x1,y1,z1);
                        vid0 = meshedges[seqVector[i]]->GetVid(0);
                        vid1 = meshedges[seqVector[i]]->GetVid(1);

                        if(fabs(x0-x1) < 1e-8)
                        {
                            //identify corners by double visit
                            if(Verts[vid0].solve == eSolveX)
                            {
                                Verts[vid0].solve = eNoSolve;
                            }
                            else
                            {
                                Verts[vid0].solve = eSolveY;
                            }
                            
                            if(Verts[vid1].solve == eSolveX)
                            {
                                Verts[vid1].solve = eNoSolve;
                            }
                            else
                            {
                                Verts[vid1].solve = eSolveY;
                            }
                        }

                        if(fabs(y0-y1) < 1e-8)
                        {
                            //identify corners by double visit
                            if(Verts[vid0].solve == eSolveY)
                            {
                                Verts[vid0].solve = eNoSolve;
                            }
                            else
                            {
                                Verts[vid0].solve = eSolveX;
                            }
                            
                            if(Verts[vid1].solve == eSolveY)
                            {
                                Verts[vid1].solve = eNoSolve;
                            }
                            else
                            {
                                Verts[vid1].solve = eSolveX;
                            }
                        }
                    }

                }
                break;
                
                case 'T':  case 'Q':  // do nothing
                {
                    break;
                }
            default:
                NEKERROR(ErrorUtil::efatal, (std::string("Unrecognized composite token: ") + compositeElementStr).c_str());
            }
            
        }
        catch(...)
        {
            NEKERROR(ErrorUtil::efatal,
                     (std::string("Unable to read COMPOSITE data for composite: ") + compositeStr).c_str());
        }
        
        /// Keep looking
        composite = composite->NextSiblingElement("C");
    }
}

void   RedefineVertices(TiXmlElement *doc, 
                        Array<OneD,NekDouble> &dvertx, 
                        Array<OneD,NekDouble> &dverty)
{


    TiXmlElement* element = doc->FirstChildElement("VERTEX");
    ASSERTL0(element, "Unable to find mesh VERTEX tag in file.");
    
    TiXmlElement *vertex = element->FirstChildElement("V");
    
    int indx;
    int nextVertexNumber = -1;
    int err;    /// Error value returned by TinyXML.
    
    vector<NekDouble> xpts,ypts,zpts;
    NekDouble xval, yval, zval;
        
    NekDouble yoffset = 0.0;
    while (vertex)
    {
        nextVertexNumber++;
        
        TiXmlAttribute *vertexAttr = vertex->FirstAttribute();
        std::string attrName(vertexAttr->Name());

        ASSERTL0(attrName == "ID", (std::string("Unknown attribute name: ") + attrName).c_str());
        
        err = vertexAttr->QueryIntValue(&indx);
        ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");
        
        // Now read body of vertex
        std::string vertexBodyStr;
        
        TiXmlNode *vertexBody = vertex->FirstChild();

        while (vertexBody)
        {
            // Accumulate all non-comment body data.
            if (vertexBody->Type() == TiXmlNode::TEXT)
            {
                vertexBodyStr += vertexBody->ToText()->Value();
                vertexBodyStr += " ";
            }
            
            vertexBody = vertexBody->NextSibling();
        }

        ASSERTL0(!vertexBodyStr.empty(), "Vertex definitions must contain vertex data.");

        // Get vertex data from the data string.
        std::istringstream vertexDataStrm(vertexBodyStr.c_str());

        try
        {
            while(!vertexDataStrm.fail())
            {
                vertexDataStrm >> xval >> yval >> zval;                
            }

            xval += dvertx[indx];
            yval += dverty[indx];
            
            stringstream s;
            s << scientific << setprecision(8) 
              << xval << " " << yval << " " << zval;

            vertex->ReplaceChild(vertex->FirstChild(), TiXmlText(s.str()));
        }
        catch(...)
        {
            ASSERTL0(false, "Unable to read VERTEX data.");
        }   
        vertex = vertex->NextSiblingElement("V");
    }

}

void EnforceRotationalSymmetry(SpatialDomains::MeshGraphSharedPtr &mesh,
                               Array<OneD, NekDouble> &dvertx,
                               Array<OneD, NekDouble> &dverty)
{
    int i,j;
    int nverts = mesh->GetNvertices();
    Array<OneD, NekDouble> x(nverts),y(nverts);
    NekDouble xval,yval,zval;

    for(i = 0; i < nverts; ++i)
    {
        mesh->GetVertex(i)->GetCoords(xval,yval,zval);
        x[i] = xval + dvertx[i];
        y[i] = yval + dverty[i];
    }
    
    NekDouble xmax = Vmath::Vmax(nverts,x,1);
    NekDouble tol = 1e-5, dist2,xrot,yrot;
    Array<OneD,int> index(nverts);
    // find nearest 
    for(i = 0; i < nverts; ++i)
    {
        xrot = -x[i] + xmax;
        yrot = -y[i];
        tol  = 1.0;

        for(j = 0; j < nverts; ++j)
        {
            dist2 = (x[j]-xrot)*(x[j]-xrot) + (y[j]-yrot)*(y[j]-yrot);
            if(dist2 < tol)
            {
                index[i]  = j;
                tol = dist2;
            }
        }
    }

    //average points and recalcualte dvertx, dverty
    for(i = 0; i < nverts; ++i)
    {
        mesh->GetVertex(i)->GetCoords(xval,yval,zval);

        xrot = 0.5*(-x[index[i]] + xmax + x[i]);
        yrot = 0.5*(-y[index[i]] + y[i]);
        
        dvertx[i] = xrot - xval;
        dverty[i] = yrot - yval;
    }
}
