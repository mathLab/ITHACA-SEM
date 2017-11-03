///////////////////////////////////////////////////////////////////////////////
//
// File LegacyFilter.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Base class for legacy filters.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <SolverUtils/Filters/LegacyFilter.h>

namespace Nektar
{
namespace SolverUtils
{

LegacyFilter::LegacyFilter(const LibUtilities::SessionReaderSharedPtr &pSession)
    : Filter(pSession)
{
}

LegacyFilter::~LegacyFilter()
{
}

void LegacyFilter::SolvedVarsOnly(
    const LibUtilities::FieldMetaDataMap &fieldMetaDataMap,
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    Array<OneD, MultiRegions::ExpListSharedPtr> &sFields)
{
    if (fieldMetaDataMap.empty())
    {
        sFields = pFields;
    }
    else
    {
        std::vector<std::string> variables;
        LibUtilities::FieldMetaDataMap tmp = fieldMetaDataMap;
        ParseUtils::GenerateVector(tmp["Variables"],
                                   variables);

        sFields = Array<OneD, MultiRegions::ExpListSharedPtr>(variables.size());
        for (int i = 0; i < variables.size(); ++i)
        {
            sFields[i] = pFields[i];
        }
    }

}
}
}
