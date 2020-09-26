////////////////////////////////////////////////////////////////////////////////
//
// File: SessionFunction.cpp
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

#include <FieldUtils/Interpolator.h>
#include <SolverUtils/Core/SessionFunction.h>

#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <boost/format.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

/**
 * Representation of a FUNCTION defined in the session xml file.
 *
 * @param   session       The session where the function was defined.
 * @param   field         The field the function is defined on.
 * @param   functionName  The name of the function.
 * @param   toCache       Store the evaluated function for later use.
 */
SessionFunction::SessionFunction(const LibUtilities::SessionReaderSharedPtr &session,
                                 const MultiRegions::ExpListSharedPtr &field,
                                 std::string functionName,
                                 bool toCache)
    : m_session(session), m_field(field), m_name(functionName), m_toCache(toCache)
{
    ASSERTL0(m_session->DefinesFunction(m_name),
             "Function '" + m_name + "' does not exist.");
}

/**
 * Evaluates a function defined in the xml session file at each quadrature point.
 *
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
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
 * Evaluates a function defined in the xml session file at each quadrature point.
 *
 * @param   pFieldNames  The names of the fields to evaluate the function for.
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::Evaluate(std::vector<std::string> pFieldNames,
                               Array<OneD, Array<OneD, NekDouble> > &pArray,
                               const NekDouble &pTime,
                               const int domain)
{
    ASSERTL1(pFieldNames.size() == pArray.size(),
             "Function '" + m_name +
             "' variable list size mismatch with array storage.");

    for (int i = 0; i < pFieldNames.size(); i++)
    {
        Evaluate(pFieldNames[i], pArray[i], pTime, domain);
    }
}

/**
 * Evaluates a function defined in the xml session file at each quadrature point.
 *
 * @param   pFieldNames  The names of the fields to evaluate the function for.
 * @param   pFields      The fields into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::Evaluate(
    std::vector<std::string> pFieldNames,
    Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &pTime,
    const int domain)
{
    ASSERTL0(pFieldNames.size() == pFields.size(),
             "Field list / name list size mismatch.");

    for (int i = 0; i < pFieldNames.size(); i++)
    {
        Evaluate(pFieldNames[i], pFields[i]->UpdatePhys(), pTime, domain);
        pFields[i]->FwdTrans_IterPerExp(pFields[i]->GetPhys(),
                                        pFields[i]->UpdateCoeffs());
    }
}

/**
 * Evaluates a function defined in the xml session file at each quadrature point.
 *
 * @param   pFieldName   The name of the field to evaluate the function for.
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::Evaluate(std::string pFieldName,
                               Array<OneD, NekDouble> &pArray,
                               const NekDouble &pTime,
                               const int domain)
{
    LibUtilities::FunctionType vType =
        m_session->GetFunctionType(m_name, pFieldName, domain);

    unsigned int nq = m_field->GetNpoints();

    std::pair<std::string, int> key(pFieldName, domain);
    // sorry
    if ((m_arrays.find(key) != m_arrays.end()) &&
            (vType == LibUtilities::eFunctionTypeFile ||
             ((m_lastCached.find(key) != m_lastCached.end()) &&
              (pTime - m_lastCached[key] < NekConstants::kNekZeroTol))))
    {
        // found cached field
        if (pArray.size() < nq)
        {
            pArray = Array<OneD, NekDouble>(nq);
        }
        Vmath::Vcopy(nq, m_arrays[key], 1, pArray, 1);

        return;
    }

    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        EvaluateExp(pFieldName, pArray, pTime, domain);
    }
    else if (vType == LibUtilities::eFunctionTypeFile ||
             vType == LibUtilities::eFunctionTypeTransientFile)
    {
        std::string filename =
            m_session->GetFunctionFilename(m_name, pFieldName, domain);

        if (boost::filesystem::path(filename).extension() == ".pts" ||
            boost::filesystem::path(filename).extension() == ".csv")
        {
            EvaluatePts(pFieldName, pArray, pTime, domain);
        }
        else
        {
            EvaluateFld(pFieldName, pArray, pTime, domain);
        }
    }
    else
    {
        ASSERTL0(false, "unknown eFunctionType");
    }

    if (m_toCache)
    {
        m_arrays[key] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, pArray, 1, m_arrays[key], 1);

        m_lastCached[key] = pTime;
    }
}

/**
 * @brief Provide a description of a function for a given field name.
 *
 * @param pFieldName     Field name.
 * @param domain         The domain to evaluate the function in.
 */
std::string SessionFunction::Describe(std::string pFieldName, const int domain)
{
    std::string retVal;

    LibUtilities::FunctionType vType =
        m_session->GetFunctionType(m_name, pFieldName, domain);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        LibUtilities::EquationSharedPtr ffunc =
            m_session->GetFunction(m_name, pFieldName, domain);
        retVal = ffunc->GetExpression();
    }
    else if (vType == LibUtilities::eFunctionTypeFile ||
             vType == LibUtilities::eFunctionTypeTransientFile)
    {
        std::string filename =
            m_session->GetFunctionFilename(m_name, pFieldName, domain);
        retVal = "from file " + filename;
    }
    else
    {
        ASSERTL0(false, "unknown eFunctionType");
    }

    return retVal;
}

/**
 * Evaluates a function from expression
 *
 * @param   pFieldName   The name of the field to evaluate the function for.
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::EvaluateExp(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.size() < nq)
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
        m_session->GetFunction(m_name, pFieldName, domain);

    ffunc->Evaluate(x0, x1, x2, pTime, pArray);
}

/**
 * Evaluates a function from fld file
 *
 * @param   pFieldName   The name of the field to evaluate the function for.
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::EvaluateFld(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.size() < nq)
    {
        pArray = Array<OneD, NekDouble>(nq);
    }

    std::string filename =
        m_session->GetFunctionFilename(m_name, pFieldName, domain);
    std::string fileVar =
        m_session->GetFunctionFilenameVariable(m_name, pFieldName, domain);

    if (fileVar.length() == 0)
    {
        fileVar = pFieldName;
    }

    //  In case of eFunctionTypeTransientFile, generate filename from
    //  format string
    LibUtilities::FunctionType vType =
        m_session->GetFunctionType(m_name, pFieldName, domain);
    if (vType == LibUtilities::eFunctionTypeTransientFile)
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
                     "Invalid Filename in function \"" + m_name +
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

/**
 * Evaluates a function from pts file
 *
 * @param   pFieldName   The name of the field to evaluate the function for.
 * @param   pArray       The array into which to write the values.
 * @param   pTime        The time at which to evaluate the function.
 * @param   domain       The domain to evaluate the function in.
 */
void SessionFunction::EvaluatePts(string pFieldName,
                                  Array<OneD, NekDouble> &pArray,
                                  const NekDouble &pTime,
                                  const int domain)
{
    unsigned int nq = m_field->GetNpoints();
    if (pArray.size() < nq)
    {
        pArray = Array<OneD, NekDouble>(nq);
    }

    std::string filename =
        m_session->GetFunctionFilename(m_name, pFieldName, domain);
    std::string fileVar =
        m_session->GetFunctionFilenameVariable(m_name, pFieldName, domain);

    if (fileVar.length() == 0)
    {
        fileVar = pFieldName;
    }

    //  In case of eFunctionTypeTransientFile, generate filename from
    //  format string
    LibUtilities::FunctionType vType =
        m_session->GetFunctionType(m_name, pFieldName, domain);
    if (vType == LibUtilities::eFunctionTypeTransientFile)
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
                     "Invalid Filename in function \"" + m_name +
                     "\", variable \"" + fileVar + "\"")
        }
    }

    LibUtilities::PtsFieldSharedPtr outPts;
    // check if we already loaded this file. For transient files,
    // funcFilename != filename so we can make sure we only keep the
    // latest pts field per funcFilename.
    std::string funcFilename =
        m_session->GetFunctionFilename(m_name, pFieldName, domain);

    LibUtilities::PtsFieldSharedPtr inPts;
    if (boost::filesystem::path(filename).extension() == ".pts")
    {
        LibUtilities::PtsIO ptsIO(m_session->GetComm());
        ptsIO.Import(filename, inPts);
    }
    else if (boost::filesystem::path(filename).extension() == ".csv")
    {
        LibUtilities::CsvIO csvIO(m_session->GetComm());
        csvIO.Import(filename, inPts);
    }
    else
    {
        ASSERTL1(false, "Unsupported file type");
    }


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

    FieldUtils::Interpolator interp;
    if (m_interpolators.find(funcFilename) != m_interpolators.end())
    {
        interp = m_interpolators[funcFilename];
    }
    else
    {
        interp = FieldUtils::Interpolator(LibUtilities::eShepard);
        if (m_session->GetComm()->GetRank() == 0)
        {
            interp.SetProgressCallback(&SessionFunction::PrintProgressbar,
                                       this);
        }
        interp.CalcWeights(inPts, outPts);
        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << endl;
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                interp.PrintStatistics();
            }
        }
    }

    if (m_toCache)
    {
        m_interpolators[funcFilename] = interp;
    }

    // TODO: only interpolate the field we actually want
    interp.Interpolate(inPts, outPts);

    int fieldInd;
    vector<string> fieldNames = outPts->GetFieldNames();
    for (fieldInd = 0; fieldInd < fieldNames.size(); ++fieldInd)
    {
        if (outPts->GetFieldName(fieldInd) == fileVar)
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
