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

#include <iostream>
#include <string>

std::string errMsg =
"MeshConvert\n"
"-----------\n"
"If you're seeing this message, you tried to run the MeshConvert executable.\n\n"
"From version 4.3 onwards, MeshConvert has been renamed NekMesh. This has\n"
"been done because MeshConvert now incorporates mesh generation features,\n"
"which takes it beyond the realm of a converter.\n\n"
"Note that the usage of NekMesh is functionally identical to MeshConvert. You\n"
"simply need to change the program name to NekMesh.\n\n"
"This placeholder program will be removed in v4.4.\n";

int main()
{
    std::cout << errMsg;
    return 1;
}
