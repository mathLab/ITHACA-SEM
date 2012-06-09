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

#include <tinyxml/tinyxml.h>
#include <SpatialDomains/MeshGraph.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList2D.h>


int main(int argc, char *argv[])
{
    int i,j;
    NekDouble cr = 0;
    
    if(argc !=3)
    {
        fprintf(stderr,"Usage: ./MoveMeshToCriticalLayer  meshfile streakfile  \n");
        exit(1);
    }
    
    //------------------------------------------------------------
    // Create Session file. 
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);
    //-----------------------------------------------------------
    
    //-------------------------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    TiXmlDocument doc(meshfile);
    bool loadOkay = doc.LoadFile();
    SpatialDomains::MeshGraph mesh;

    std::stringstream errstr;
    errstr << "Unable to load file: " << meshfile << " (";
    errstr << doc.ErrorDesc() << ", line " << doc.ErrorRow()
           << ", column " << doc.ErrorCol() << ")";
    ASSERTL0(loadOkay, errstr.str());
    
    mesh.ReadGeometry(doc);
    //------------------------------------------------------------

    std::string outfile = meshfile + ".new";
    doc.SaveFile(outfile);
}

