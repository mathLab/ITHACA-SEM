///////////////////////////////////////////////////////////////////////////////
//
// File PtsIO.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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
// Description: Pts IO prototype definitions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_PTSIO_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_PTSIO_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <tinyxml.h>

#include <LibUtilities/Communication/Comm.h>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/FieldIO.h>

namespace Nektar
{
namespace LibUtilities
{
using namespace std;

typedef std::map<std::string, std::string> PtsMetaDataMap;
static PtsMetaDataMap NullPtsMetaDataMap;

LIB_UTILITIES_EXPORT void Import(const string &inFile,
                                 PtsFieldSharedPtr &ptsField);

LIB_UTILITIES_EXPORT void Write(const string &outFile,
                                const PtsFieldSharedPtr &ptsField);

class PtsIO : protected FieldIO
{
public:
    LIB_UTILITIES_EXPORT PtsIO(LibUtilities::CommSharedPtr pComm,
                               bool sharedFilesystem = false);

    LIB_UTILITIES_EXPORT void Import(
        const string &inFile,
        PtsFieldSharedPtr &ptsField,
        FieldMetaDataMap &fieldmetadatamap = NullFieldMetaDataMap);

    LIB_UTILITIES_EXPORT void Write(const string &outFile,
                                    const PtsFieldSharedPtr &ptsField);

    LIB_UTILITIES_EXPORT void ImportFieldData(TiXmlDocument docInput,
                                              PtsFieldSharedPtr &ptsField);

protected:
    LIB_UTILITIES_EXPORT void SetUpFieldMetaData(const std::string outname);

    LIB_UTILITIES_EXPORT virtual std::string GetFileEnding() const
    {
        return "pts";
    };
};

typedef boost::shared_ptr<PtsIO> PtsIOSharedPtr;
}
}
#endif
