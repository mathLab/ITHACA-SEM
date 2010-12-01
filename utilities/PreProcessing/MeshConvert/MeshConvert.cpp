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

#include "Convert.h"
using namespace Nektar::Utilities;

int main(int argc, char* argv[]) {
    if (argc != 3)
    {
        cout << "Usage: " << argv[0] << " [msh-file] [xml-output-file]" << endl;
        cout << "  [msh-file] - mesh output from Gmsh." << endl;
        cout << "  [xml-output-file] - output file." << endl;
        return 1;
    }

    string vFilename = argv[1];
    int vDot = vFilename.find_last_of('.');
    string ext = vFilename.substr(++vDot, vFilename.length() - vDot);

    boost::shared_ptr<Convert> G = ConvertFactory::CreateInstance(ext);

    G->ReadFile(vFilename);
    G->Process();
    G->WriteXmlFile(string(argv[2]));

    return 0;
}
