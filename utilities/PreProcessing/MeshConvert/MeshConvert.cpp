////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshConvert.cpp
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
//  Description: Mesh conversion utility.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "Module.h"
using namespace Nektar::Utilities;

int main(int argc, char* argv[]) {
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " [msh-file] [xml-output-file]" << endl;
        cout << "  [msh-file] - mesh output from Gmsh." << endl;
        cout << "  [xml-output-file] - output file." << endl;
        return 1;
    }
    
    string inFilename  = argv[1];
    string outFilename = argv[2];
    int    inDot       = inFilename .find_last_of('.');
    int    outDot      = outFilename.find_last_of('.');
    string inExt       = inFilename .substr(++inDot,  inFilename .length() - inDot );
    string outExt      = outFilename.substr(++outDot, outFilename.length() - outDot);
    MeshSharedPtr m    = boost::shared_ptr<Mesh>(new Mesh(inFilename, outFilename));

    ModuleSharedPtr in  = GetModuleFactory().CreateInstance(ModuleKey(inExt, eInputModule ), m);
    ModuleSharedPtr out = GetModuleFactory().CreateInstance(ModuleKey(outExt,eOutputModule), m);
    ModuleSharedPtr pr  = GetModuleFactory().CreateInstance(ModuleKey("bl",  eProcessModule),m);
    //ModuleSharedPtr pr2 = GetModuleFactory().CreateInstance(ModuleKey("jac", eProcessModule),m);
    in -> Process();
    pr -> Process();
    //pr2-> Process();
    out-> Process();
    
    return 0;
}
