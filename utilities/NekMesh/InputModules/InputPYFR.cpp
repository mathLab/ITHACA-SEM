////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNekpp.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////


#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "InputPYFR.h"

using namespace std;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputPYFR::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pyfrm"), InputPYFR::create,
        "PYFR mesh files.");

/**
 * @brief Set up InputNekpp object.
 */
InputPYFR::InputPYFR(MeshSharedPtr m) : InputModule(m)
{

}

InputPYFR::~InputPYFR()
{
    string filename = m_config["infile"].as<string>();

    H5std_string FILE_NAME( filename );
    H5File file( FILE_NAME, H5F_ACC_RDONLY );

    Group g = file.openGroup("/");
    for(int i = 0; i < g.getNumObjs(); i++)
    {
        cout << g.getObjnameByIdx(i) << endl;
    }

}

void InputPYFR::Process()
{

}

}
}
