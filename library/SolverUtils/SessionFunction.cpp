////////////////////////////////////////////////////////////////////////////////
//
// File: PtsField.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Session Function
//
////////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/SessionFunction.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

SessionFunction::SessionFunction(LibUtilities::SessionReaderSharedPtr session,
                                 MultiRegions::ExpListSharedPtr field,
                                 std::string functionName)
    : m_session(session), m_field(field), m_functionName(functionName)
{
    ASSERTL0(m_session->DefinesFunction(m_functionName),
             "Function '" + m_functionName + "' does not exist.");

    // TODO: this is a hack because GetFunctionType wants a fieldname although
    // it doesnt need one
    std::string pFieldName = "*";

    m_type = m_session->GetFunctionType(m_functionName, pFieldName);
}

/**
 * Evaluates a physical function at each quadrature point in the domain.
 *
 * @param   pArray          The array into which to write the values.
 * @param   pEqn            The equation to evaluate.
 */
void SessionFunction::Evaluate(Array<OneD, Array<OneD, NekDouble> > &pArray,
                               const NekDouble pTime,
                               const int domain)
{
    std::vector<std::string> vFieldNames = m_session->GetVariables();

    for (int i = 0; i < vFieldNames.size(); i++)
    {
        Evaluate(vFieldNames[i], pArray[i], pTime, domain);
    }
}

/**
* Populates a forcing function for each of the dependent variables
* using the expression provided by the BoundaryConditions object.
* @param   force           Array of fields to assign forcing.
*/
void SessionFunction::Evaluate(std::vector<std::string> pFieldNames,
                               Array<OneD, Array<OneD, NekDouble> > &pArray,
                               const NekDouble &pTime,
                               const int domain)
{
    ASSERTL1(pFieldNames.size() == pArray.num_elements(),
             "Function '" + m_functionName +
                 "' variable list size mismatch with array storage.");

    for (int i = 0; i < pFieldNames.size(); i++)
    {
        Evaluate(pFieldNames[i], pArray[i], pTime, domain);
    }
}

/**
* Populates a function for each of the dependent variables using
* the expression or filenames provided by the SessionReader object.
* @param   force           Array of fields to assign forcing.
*/
void SessionFunction::Evaluate(
    std::vector<std::string> pFieldNames,
    Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &pTime,
    const int domain)
{
    ASSERTL0(pFieldNames.size() == pFields.num_elements(),
             "Field list / name list size mismatch.");

    for (int i = 0; i < pFieldNames.size(); i++)
    {
        Evaluate(pFieldNames[i], pFields[i]->UpdatePhys(), pTime, domain);
        pFields[i]->FwdTrans_IterPerExp(pFields[i]->GetPhys(),
                                        pFields[i]->UpdateCoeffs());
    }
}

void SessionFunction::Evaluate(std::string pFieldName,
                               Array<OneD, NekDouble> &pArray,
                               const NekDouble &pTime,
                               const int domain)
{
    if (m_type == LibUtilities::eFunctionTypeExpression)
    {
        EvaluateExp(pFieldName, pArray, pTime, domain);
    }
    else if (m_type == LibUtilities::eFunctionTypeFile ||
             m_type == LibUtilities::eFunctionTypeTransientFile)
    {
        std::string filename =
            m_session->GetFunctionFilename(m_functionName, pFieldName, domain);

        if (boost::filesystem::path(filename).extension() == ".pts")
        {
            EvaluatePts(pFieldName, pArray, pTime, domain);
        }
        else
        {
            EvaluateFld(pFieldName, pArray, pTime, domain);
        }
    }
}

/**
         * @brief Provide a description of a function for a given field name.
         *
         * @param pFieldName     Field name.
         * @param pFunctionName  Function name.
         */
std::string SessionFunction::Describe(std::string pFieldName,
                                              const int domain)
{
    std::string retVal;
    if (m_type == LibUtilities::eFunctionTypeExpression)
    {
        LibUtilities::EquationSharedPtr ffunc =
            m_session->GetFunction(m_functionName, pFieldName, domain);
        retVal = ffunc->GetExpression();
    }
    else if (m_type == LibUtilities::eFunctionTypeFile ||
             LibUtilities::eFunctionTypeTransientFile)
    {
        std::string filename =
            m_session->GetFunctionFilename(m_functionName, pFieldName, domain);
        retVal = "from file " + filename;
    }

    return retVal;
}

void SessionFunction::EvaluateExp(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.num_elements() < nq)
    {
        pArray = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> x2(nq);

    // Get the coordinates (assuming all fields have the same
    // discretisation)
    m_field->GetCoords(x0, x1, x2);
    LibUtilities::EquationSharedPtr ffunc =
        m_session->GetFunction(m_functionName, pFieldName, domain);

    ffunc->Evaluate(x0, x1, x2, pTime, pArray);
}

void SessionFunction::EvaluateFld(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.num_elements() < nq)
    {
        pArray = Array<OneD, NekDouble>(nq);
    }

    std::string filename =
        m_session->GetFunctionFilename(m_functionName, pFieldName, domain);
    std::string fileVar = m_session->GetFunctionFilenameVariable(
        m_functionName, pFieldName, domain);

    if (fileVar.length() == 0)
    {
        fileVar = pFieldName;
    }

    //  In case of eFunctionTypeTransientFile, generate filename from
    //  format string
    if (m_type == LibUtilities::eFunctionTypeTransientFile)
    {
        try
        {
#if (defined _WIN32 && _MSC_VER < 1900)
            // We need this to make sure boost::format has always
            // two digits in the exponents of Scientific notation.
            unsigned int old_exponent_format;
            old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
            filename            = boost::str(boost::format(filename) % pTime);
            _set_output_format(old_exponent_format);
#else
            filename = boost::str(boost::format(filename) % pTime);
#endif
        }
        catch (...)
        {
            ASSERTL0(false,
                     "Invalid Filename in function \"" + m_functionName +
                         "\", variable \"" + fileVar + "\"")
        }
    }

    // Define list of global element ids
    int numexp = m_field->GetExpSize();
    Array<OneD, int> ElementGIDs(numexp);
    for (int i = 0; i < numexp; ++i)
    {
        ElementGIDs[i] = m_field->GetExp(i)->GetGeom()->GetGlobalID();
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> > FieldData;
    LibUtilities::FieldIOSharedPtr fldIO =
        LibUtilities::FieldIO::CreateForFile(m_session, filename);
    fldIO->Import(filename,
                  FieldDef,
                  FieldData,
                  LibUtilities::NullFieldMetaDataMap,
                  ElementGIDs);

    int idx = -1;
    Array<OneD, NekDouble> vCoeffs(m_field->GetNcoeffs(), 0.0);
    // Loop over all the expansions
    for (int i = 0; i < FieldDef.size(); ++i)
    {
        // Find the index of the required field in the
        // expansion segment
        for (int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
        {
            if (FieldDef[i]->m_fields[j] == fileVar)
            {
                idx = j;
            }
        }

        if (idx >= 0)
        {
            m_field->ExtractDataToCoeffs(
                FieldDef[i], FieldData[i], FieldDef[i]->m_fields[idx], vCoeffs);
        }
        else
        {
            cout << "Field " + fileVar + " not found." << endl;
        }
    }

    m_field->BwdTrans_IterPerExp(vCoeffs, pArray);
}

void SessionFunction::EvaluatePts(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.num_elements() < nq)
    {
        pArray = Array<OneD, NekDouble>(nq);
    }

    std::string filename =
        m_session->GetFunctionFilename(m_functionName, pFieldName, domain);
    std::string fileVar = m_session->GetFunctionFilenameVariable(
        m_functionName, pFieldName, domain);

    if (fileVar.length() == 0)
    {
        fileVar = pFieldName;
    }

    //  In case of eFunctionTypeTransientFile, generate filename from
    //  format string
    if (m_type == LibUtilities::eFunctionTypeTransientFile)
    {
        try
        {
#if (defined _WIN32 && _MSC_VER < 1900)
            // We need this to make sure boost::format has always
            // two digits in the exponents of Scientific notation.
            unsigned int old_exponent_format;
            old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
            filename            = boost::str(boost::format(filename) % pTime);
            _set_output_format(old_exponent_format);
#else
            filename = boost::str(boost::format(filename) % pTime);
#endif
        }
        catch (...)
        {
            ASSERTL0(false,
                     "Invalid Filename in function \"" + m_functionName +
                         "\", variable \"" + fileVar + "\"")
        }
    }

    LibUtilities::PtsFieldSharedPtr outPts;
    // check if we already loaded this file. For transient files,
    // funcFilename != filename so we can make sure we only keep the
    // latest pts field per funcFilename.
    std::string funcFilename =
        m_session->GetFunctionFilename(m_functionName, pFieldName, domain);

    LibUtilities::PtsFieldSharedPtr inPts;
    LibUtilities::PtsIO ptsIO(m_session->GetComm());
    ptsIO.Import(filename, inPts);

    Array<OneD, Array<OneD, NekDouble> > pts(inPts->GetDim() +
                                             inPts->GetNFields());
    for (int i = 0; i < inPts->GetDim() + inPts->GetNFields(); ++i)
    {
        pts[i] = Array<OneD, NekDouble>(nq);
    }
    if (inPts->GetDim() == 1)
    {
        m_field->GetCoords(pts[0]);
    }
    else if (inPts->GetDim() == 2)
    {
        m_field->GetCoords(pts[0], pts[1]);
    }
    else if (inPts->GetDim() == 3)
    {
        m_field->GetCoords(pts[0], pts[1], pts[2]);
    }
    outPts = MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(
        inPts->GetDim(), inPts->GetFieldNames(), pts);

    if (!m_interpolator.HasWeights())
    {
        m_interpolator = FieldUtils::Interpolator(Nektar::FieldUtils::eShepard);
        if (m_session->GetComm()->GetRank() == 0)
        {
            m_interpolator.SetProgressCallback(
                &SessionFunction::PrintProgressbar, this);
        }
        m_interpolator.CalcWeights(inPts, outPts);
        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << endl;
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                m_interpolator.PrintStatistics();
            }
        }
    }

    // TODO: only interpolate the field we actually want
    m_interpolator.Interpolate(inPts, outPts);

    int fieldInd;
    vector<string> fieldNames = outPts->GetFieldNames();
    for (fieldInd = 0; fieldInd < fieldNames.size(); ++fieldInd)
    {
        if (outPts->GetFieldName(fieldInd) == pFieldName)
        {
            break;
        }
    }
    ASSERTL0(fieldInd != fieldNames.size(), "field not found");

    pArray = outPts->GetPts(fieldInd + outPts->GetDim());
}

// end of namespaces
}
}
