////////////////////////////////////////////////////////////////////////////////
//
// File: PtsIO.cpp
//
// For more information, please see: http://www.nektar.info/
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
// Description: I/O routines relating to Pts data
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/PtsIO.h>

namespace Nektar
{
namespace LibUtilities
{

void Import(const string &inFile, PtsFieldSharedPtr &ptsField)
{
    PtsIO p;
    p.Import(inFile, ptsField);
}


void Write(const string &outFile, const PtsFieldSharedPtr &ptsField)
{
    PtsIO p;
    p.Write(outFile, ptsField);
}


/**
 * @brief Import a pts field from file
 *
 * @param inFile    filename of the file to read
 * @param ptsField  the resulting pts field.
 */
void PtsIO::Import(const string &inFile, PtsFieldSharedPtr &ptsField)
{
    TiXmlDocument docInput;
    ASSERTL0(docInput.LoadFile(inFile), "Unable to open file '" + inFile + "'.");

    TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");
    TiXmlElement *points = nektar->FirstChildElement("POINTS");
    int dim;
    int err = points->QueryIntAttribute("DIM", &dim);

    ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute DIM.");

    std::string fields = points->Attribute("FIELDS");

    vector<string> fieldNames;
    if (!fields.empty())
    {
        bool valid = ParseUtils::GenerateOrderedStringVector(
                         fields.c_str(), fieldNames);
        ASSERTL0(valid,
                 "Unable to process list of field variable in  FIELDS attribute:  " + fields);
    }

    int nfields = fieldNames.size();
    int totvars = dim + nfields;

    TiXmlNode *pointsBody = points->FirstChild();

    std::istringstream pointsDataStrm(pointsBody->ToText()->Value());

    vector<NekDouble> ptsSerial;
    Array<OneD,  Array<OneD,  NekDouble> > pts(totvars);

    try
    {
        NekDouble      ptsStream;
        while (!pointsDataStrm.fail())
        {
            pointsDataStrm >> ptsStream;

            ptsSerial.push_back(ptsStream);
        }
    }
    catch (...)
    {
        ASSERTL0(false, "Unable to read Points data.");
    }

    int npts = ptsSerial.size() / totvars;

    for (int i = 0; i < totvars; ++i)
    {
        pts[i] = Array<OneD, NekDouble>(npts);
    }

    for (int i = 0; i < npts; ++i)
    {
        for (int j = 0; j < totvars; ++j)
        {
            pts[j][i] = ptsSerial[i * totvars + j];
        }
    }

    ptsField = MemoryManager<PtsField>::AllocateSharedPtr(dim, fieldNames, pts);
}


/**
 * @brief Save a pts field to a file
 *
 * @param inFile    filename of the file
 * @param ptsField  the pts field
 */
void PtsIO::Write(const string &outFile, const PtsFieldSharedPtr &ptsField)
{

    ASSERTL0(false, "Not implemented yet");
}


}
}
