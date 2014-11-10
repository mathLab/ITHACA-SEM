///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_XML_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_XML_H

#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    class XmlUtil
    {
        public:
            LIB_UTILITIES_EXPORT static const char* NEKTAR_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* GEOMETRY_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* DIM_GEOMETRY_ATTRIBUTE;
            LIB_UTILITIES_EXPORT static const char* SPACE_GEOMETRY_ATTRIBUTE;
            LIB_UTILITIES_EXPORT static const char* VERTEX_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* V_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* ID_VERTEX_ATTRIBUTE_NAME;
            LIB_UTILITIES_EXPORT static const char* EDGE_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* E_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* CURVED_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* FACE_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* COMPOSITE_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* C_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* DOMAIN_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* ELEMENT_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* S_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* T_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* Q_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* A_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* P_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* R_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* H_ELEMENT_NAME;            
            LIB_UTILITIES_EXPORT static const char* PARAMETERS_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* SOLVERINFO_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* VARIABLES_ELEMENT_NAME;
            LIB_UTILITIES_EXPORT static const char* BOUNDARYCONDITIONS_ELEMENT_NAME;

    };
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_XML_H
