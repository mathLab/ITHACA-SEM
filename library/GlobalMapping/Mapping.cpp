///////////////////////////////////////////////////////////////////////////////
//
// File: Mapping.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Abstract base class for mappings.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <GlobalMapping/Mapping.h>

#include <boost/algorithm/string/predicate.hpp>

using namespace std;

namespace Nektar
{
namespace GlobalMapping
{

MappingSharedPtr Mapping::m_mappingPtr = MappingSharedPtr();
bool             Mapping::m_init       = false;
bool             Mapping::m_isDefined  = false;

MappingFactory& GetMappingFactory()
{
    static MappingFactory instance;
    return instance;
}

Mapping::Mapping(const LibUtilities::SessionReaderSharedPtr& pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields)
    : m_session(pSession), m_fields(pFields)
{
    switch (m_fields[0]->GetExpType())
    {
        case MultiRegions::e1D:
        {
            m_nConvectiveFields = 1;
        }
        break;

        case MultiRegions::e2D:
        {
            m_nConvectiveFields = 2;
        }
        break;

        case MultiRegions::e3D:
        case MultiRegions::e3DH1D:
        case MultiRegions::e3DH2D:
        {
            m_nConvectiveFields = 3;
        }
        break;

        default:
            ASSERTL0(0,"Dimension not supported");
        break;
    }

    m_fld = LibUtilities::FieldIO::CreateDefault(pSession);
}

/**
 * This function initialises the Mapping object. It computes the coordinates
 * and velocity coordinates, initialises the workspace variables, and calls
 * UpdateGeomInfo, which will perform the calculations specific for each type
 * of Mapping.
 * @param pFields ExpList array used in the mapping
 * @param pMapping xml element describing the mapping
 */
void Mapping::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement                                *pMapping)
{
    boost::ignore_unused(pFields);

    int phystot         = m_fields[0]->GetTotPoints();
    m_fromFunction      = true;
    // Initialise variables
    m_coords    = Array<OneD, Array<OneD, NekDouble> > (3);
    m_coordsVel = Array<OneD, Array<OneD, NekDouble> > (3);
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    for (int i = 0; i < 3; i++)
    {
        m_coords[i]    = Array<OneD, NekDouble> (phystot);
        m_coordsVel[i] = Array<OneD, NekDouble> (phystot);
        coords[i]      = Array<OneD, NekDouble> (phystot);
    }

    // Check if mapping is defined as time-dependent
    const TiXmlElement* timeDep = pMapping->
                                    FirstChildElement("TIMEDEPENDENT");
    if (timeDep)
    {
        string sTimeDep = timeDep->GetText();
        m_timeDependent = ( boost::iequals(sTimeDep,"true")) ||
                          ( boost::iequals(sTimeDep,"yes"));
    }
    else
    {
        m_timeDependent     = false;
    }

    // Load coordinates
    string fieldNames[3] = {"x", "y", "z"};
    const TiXmlElement* funcNameElmt = pMapping->FirstChildElement("COORDS");
    if (funcNameElmt)
    {
        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                "Function '" + m_funcName + "' not defined.");

        // Get coordinates in the domain
        m_fields[0]->GetCoords(coords[0], coords[1], coords[2]);

        std::string s_FieldStr;
        // Check if function from session file defines each component
        //      and evaluate them, otherwise use trivial transformation
        for(int i = 0; i < 3; i++)
        {
            s_FieldStr = fieldNames[i];
            if ( m_session->DefinesFunction(m_funcName, s_FieldStr))
            {
                EvaluateFunction(m_fields, m_session, s_FieldStr, m_coords[i],
                                        m_funcName);
                if ( i==2 && m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    ASSERTL0 (false,
                        "3DH1D does not support mapping in the z-direction.");
                }
            }
            else
            {
                // This coordinate is not defined, so use (x^i)' = x^i
                Vmath::Vcopy(phystot, coords[i], 1, m_coords[i], 1);
            }
        }
    }
    else
    {
        m_fields[0]->GetCoords(coords[0], coords[1], coords[2]);
        for(int i = 0; i < 3; i++)
        {
            // Use (x^i)' = x^i as default. This can be useful if we
            //    have a time-dependent mapping, and then only the
            //    initial mapping will be trivial
            Vmath::Vcopy(phystot, coords[i], 1, m_coords[i], 1);
        }
    }

    // Load coordinate velocity if they are defined,
    //      otherwise use zero to make it general
    string velFieldNames[3] = {"vx", "vy", "vz"};
    const TiXmlElement* velFuncNameElmt = pMapping->FirstChildElement("VEL");
    if (velFuncNameElmt)
    {
        m_velFuncName = velFuncNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_velFuncName),
                "Function '" + m_velFuncName + "' not defined.");

        std::string s_FieldStr;
        // Check if function from session file defines each component
        //      and evaluate them, otherwise use 0
        for(int i = 0; i < 3; i++)
        {
            s_FieldStr = velFieldNames[i];
            if ( m_session->DefinesFunction(m_velFuncName, s_FieldStr))
            {
                EvaluateFunction(m_fields, m_session, s_FieldStr,
                                    m_coordsVel[i], m_velFuncName);
                if ( i==2 && m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    ASSERTL0 (false,
                        "3DH1D does not support mapping in the z-direction.");
                }
            }
            else
            {
                // This coordinate velocity is not defined, so use 0
                Vmath::Zero(phystot, m_coordsVel[i], 1);
            }
        }
    }
    else
    {
        for(int i = 0; i < 3; i++)
        {
            Vmath::Zero(phystot, m_coordsVel[i], 1);
        }
    }

    // Initialise workspace variables
    int nvel = m_nConvectiveFields;
    m_wk1 = Array<OneD, Array<OneD, NekDouble> > (nvel*nvel);
    m_wk2 = Array<OneD, Array<OneD, NekDouble> > (nvel*nvel);
    m_tmp = Array<OneD, Array<OneD, NekDouble> > (nvel);
    for (int i=0; i< nvel; i++)
    {
        m_tmp[i] = Array<OneD, NekDouble>(phystot,0.0);
        for (int j=0; j< nvel; j++)
        {
            m_wk1[i*nvel+j] = Array<OneD, NekDouble>(phystot,0.0);
            m_wk2[i*nvel+j] = Array<OneD, NekDouble>(phystot,0.0);
        }
    }

    // Calculate information required by the particular mapping
    UpdateGeomInfo();
}

/**
 * This function replaces the expansion list contained in m_fields, and then
 * proceeds reinitialising the Mapping with this new field.
 * @param pFields New field to be used by the Mapping
 */
void Mapping::ReplaceField(
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields)
{
    m_fields = pFields;

    TiXmlElement* vMapping = NULL;

    if (m_session->DefinesElement("Nektar/Mapping"))
    {
        vMapping = m_session->GetElement("Nektar/Mapping");
    }
    InitObject(pFields, vMapping);
}

/**
 * This function is responsible for loading the Mapping, guaranteeing that a
 * single instance of this class exists. When it is first called, it creates a
 * Mapping and returns a pointer to it. On subsequent calls, it just returns
 * the pointer.
 * @param pSession Session reader object
 * @param pFields  Fields which will be used by the Mapping
 * @return Pointer to the Mapping
 */
MappingSharedPtr Mapping::Load(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields)
{
    if (!m_init)
    {
        TiXmlElement* vMapping = NULL;
        string vType;
        if (pSession->DefinesElement("Nektar/Mapping"))
        {
            vMapping = pSession->GetElement("Nektar/Mapping");
            vType = vMapping->Attribute("TYPE");
            m_isDefined = true;
        }
        else
        {
            vType = "Translation";
        }

        m_mappingPtr =   GetMappingFactory().CreateInstance(
                            vType, pSession, pFields,
                            vMapping);

        m_init = true;
    }

    return m_mappingPtr;
}

/**
 * This function should be called before writing a fld or chk file. It updates
 * the metadata with information from the Mapping. Also, if the mapping is not
 * defined by a function, it writes a .map file containing the coordinates and
 * velocity of the coordinates.
 * @param fieldMetaDataMap Metadata of the output file
 * @param outname          Name of the output file
 */
void Mapping::Output(
            LibUtilities::FieldMetaDataMap  &fieldMetaDataMap,
            const std::string                    &outname)
{
    // Only do anything if mapping exists
    if (m_isDefined)
    {
        fieldMetaDataMap["MappingCartesianVel"] = std::string("False");
        if (m_fromFunction)
        {
            // Add metadata
            fieldMetaDataMap["MappingType"] = std::string("Expression");
            fieldMetaDataMap["MappingExpression"] = m_funcName;
            if (m_timeDependent)
            {
                fieldMetaDataMap["MappingVelExpression"] = m_velFuncName;
            }
        }
        else
        {
            int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
            string fieldNames[3] = {"x", "y", "z"};
            string velFieldNames[3] = {"vx", "vy", "vz"};

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = m_fields[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

            int ncoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, NekDouble> fieldcoeffs(ncoeffs);

            bool wavespace = m_fields[0]->GetWaveSpace();
            m_fields[0]->SetWaveSpace(false);
            // copy coordinates Data into FieldData and set variable
            for(int j = 0; j < expdim; ++j)
            {
                m_fields[0]->FwdTrans_IterPerExp(m_coords[j], fieldcoeffs);

                for(int i = 0; i < FieldDef.size(); ++i)
                {
                    // Could do a search here to find correct variable
                    FieldDef[i]->m_fields.push_back(fieldNames[j]);
                    m_fields[0]->AppendFieldData(FieldDef[i], FieldData[i],
                                                                fieldcoeffs);
                }
            }
            if (m_timeDependent)
            {
                //copy coordinates velocity Data into FieldData and set variable
                for(int j = 0; j < expdim; ++j)
                {
                    m_fields[0]->FwdTrans_IterPerExp(m_coordsVel[j],
                                                        fieldcoeffs);

                    for(int i = 0; i < FieldDef.size(); ++i)
                    {
                        // Could do a search here to find correct variable
                        FieldDef[i]->m_fields.push_back(velFieldNames[j]);
                        m_fields[0]->AppendFieldData(FieldDef[i],
                                                    FieldData[i],
                                                    fieldcoeffs);
                    }
                }
            }

            std::string outfile = outname;
            outfile.erase(outfile.end()-4, outfile.end());
            outfile += ".map";

            m_fld->Write(outfile,FieldDef,FieldData,fieldMetaDataMap);

            // Write metadata to orginal output
            fieldMetaDataMap["MappingType"] = std::string("File");
            fieldMetaDataMap["FileName"] = outfile;

            m_fields[0]->SetWaveSpace(wavespace);
        }
    }
}

void Mapping::EvaluateTimeFunction(
        LibUtilities::SessionReaderSharedPtr              pSession,
        std::string                                       pFieldName,
        Array<OneD, NekDouble>&                           pArray,
        const std::string&                                pFunctionName,
        NekDouble                                         pTime)
{
    ASSERTL0(pSession->DefinesFunction(pFunctionName),
             "Function '" + pFunctionName + "' does not exist.");

    LibUtilities::EquationSharedPtr ffunc =
        pSession->GetFunction(pFunctionName, pFieldName);

    Array<OneD, NekDouble> x0(1,0.0);
    Array<OneD, NekDouble> x1(1,0.0);
    Array<OneD, NekDouble> x2(1,0.0);

    ffunc->Evaluate(x0, x1, x2, pTime, pArray);
}


void Mapping::EvaluateFunction(
        Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
        LibUtilities::SessionReaderSharedPtr              pSession,
        std::string                                       pFieldName,
        Array<OneD, NekDouble>&                           pArray,
        const std::string&                                pFunctionName,
        NekDouble                                         pTime)
{
    ASSERTL0(pSession->DefinesFunction(pFunctionName),
             "Function '" + pFunctionName + "' does not exist.");

    unsigned int nq = pFields[0]->GetNpoints();
    if (pArray.size() != nq)
    {
        pArray = Array<OneD, NekDouble> (nq);
    }

    LibUtilities::FunctionType vType;
    vType = pSession->GetFunctionType(pFunctionName, pFieldName);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        Array<OneD, NekDouble> x0(nq);
        Array<OneD, NekDouble> x1(nq);
        Array<OneD, NekDouble> x2(nq);

        pFields[0]->GetCoords(x0, x1, x2);
        LibUtilities::EquationSharedPtr ffunc =
            pSession->GetFunction(pFunctionName, pFieldName);

        ffunc->Evaluate(x0, x1, x2, pTime, pArray);
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        std::string filename = pSession->GetFunctionFilename(
                                            pFunctionName,
                                            pFieldName);

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::vector<NekDouble> > FieldData;
        Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
        Vmath::Zero(vCoeffs.size(), vCoeffs, 1);

        LibUtilities::FieldIOSharedPtr fld =
            LibUtilities::FieldIO::CreateForFile(pSession, filename);
        fld->Import(filename, FieldDef, FieldData);

        int idx = -1;
        for (int i = 0; i < FieldDef.size(); ++i)
        {
            for (int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
            {
                if (FieldDef[i]->m_fields[j] == pFieldName)
                {
                    idx = j;
                }
            }

            if (idx >= 0)
            {
                pFields[0]->ExtractDataToCoeffs(
                                            FieldDef[i],
                                            FieldData[i],
                                            FieldDef[i]->m_fields[idx],
                                            vCoeffs);
            }
            else
            {
                cout << "Field " + pFieldName + " not found." << endl;
            }
        }
        pFields[0]->BwdTrans_IterPerExp(vCoeffs, pArray);
    }
}

/**
 * This function converts a contravariant vector in transformed space
 * \f$v^{j}\f$ to the corresponding \f$\bar{v}^{i}\f$ in Cartesian (physical)
 * space using the relation
 * \f[\bar{v}^{i} = \frac{\partial \bar{x}^i}{\partial x^j}v^{j}\f]
 *
 * @param inarray  Components of the vector in transformed space (\f$v^{j}\f$)
 * @param outarray Components of the vector in Cartesian space (\f$\bar{v}^{i}\f$)
 */
void Mapping::ContravarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_ContravarToCartesian( m_tmp, outarray);
    }
    else
    {
        v_ContravarToCartesian( inarray, outarray);
    }
}

/**
 * This function converts a covariant vector in transformed space
 * \f$v_{j}\f$ to the corresponding \f$\bar{v}_{i}\f$ in Cartesian (physical)
 * space using the relation
 * \f[\bar{v}_{i} = \frac{\partial x^j}{\partial \bar{x}^i}v_{j}\f]
 *
 * @param inarray  Components of the vector in transformed space (\f$v_{j}\f$)
 * @param outarray Components of the vector in Cartesian space (\f$\bar{v}_{i}\f$)
 */
void Mapping::CovarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_CovarToCartesian( m_tmp, outarray);
    }
    else
    {
        v_CovarToCartesian( inarray, outarray);
    }
}

/**
 * This function converts a contravariant vector in Cartesian space
 * \f$\bar{v}^{i}\f$ to the corresponding \f$v^{j}\f$ in transformed
 * space using the relation
 * \f[v^{j} = \frac{\partial x^j}{\partial \bar{x}^i}\bar{v}^{i}\f]
 *
 * @param inarray  Components of the vector in Cartesian space (\f$\bar{v}^{i}\f$)
 * @param outarray Components of the vector in transformed space (\f$v^{j}\f$)
 */
void Mapping::ContravarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_ContravarFromCartesian( m_tmp, outarray);
    }
    else
    {
        v_ContravarFromCartesian( inarray, outarray);
    }
}

/**
 * This function converts a covariant vector in Cartesian space
 * \f$\bar{v}_{i}\f$ to the corresponding \f$v_{j}\f$ in transformed
 * space using the relation
 * \f[\bar{v}_{j} = \frac{\partial \bar{x}^i}{\partial x^j}v_{i}\f]
 *
 * @param inarray  Components of the vector in Cartesian space (\f$\bar{v}_{i}\f$)
 * @param outarray Components of the vector in transformed space (\f$v_{j}\f$)
 */
void Mapping::CovarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_CovarFromCartesian( m_tmp, outarray);
    }
    else
    {
        v_CovarFromCartesian( inarray, outarray);
    }
}

/**
 * This function lowers the index of the contravariant vector \f$v^{i}\f$,
 * transforming it into its associated covariant vector \f$v_{j}\f$
 *  according to the relation
 * \f[v_{j} = g_{ij}v^{i}\f]
 * where \f$g_{ij}\f$ is the metric tensor.
 *
 * @param inarray  Components of the contravariant vector \f$v^{i}\f$
 * @param outarray Components of the covariant vector \f$v_{j}\f$
 */
void Mapping::LowerIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_LowerIndex( m_tmp, outarray);
    }
    else
    {
        v_LowerIndex( inarray, outarray);
    }
}

/**
 * This function raises the index of the covariant vector \f$v_{j}\f$,
 * transforming it into its associated contravariant vector \f$v^{i}\f$
 *  according to the relation
 * \f[v^{i} = g^{ij}v_{j}\f]
 * where \f$g^{ij}\f$ is the inverse of the metric tensor.
 *
 * @param inarray  Components of the contravariant vector \f$v^{i}\f$
 * @param outarray Components of the covariant vector \f$v_{j}\f$
 */
void Mapping::RaiseIndex(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    if(inarray == outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;

        for (int i=0; i< nvel; i++)
        {
            Vmath::Vcopy(physTot, inarray[i], 1, m_tmp[i], 1);
        }
        v_RaiseIndex( m_tmp, outarray);
    }
    else
    {
        v_RaiseIndex( inarray, outarray);
    }
}

void Mapping::v_GetCartesianCoordinates(
            Array<OneD, NekDouble>               &out0,
            Array<OneD, NekDouble>               &out1,
            Array<OneD, NekDouble>               &out2)
{
    int physTot = m_fields[0]->GetTotPoints();

    out0 = Array<OneD, NekDouble>(physTot, 0.0);
    out1 = Array<OneD, NekDouble>(physTot, 0.0);
    out2 = Array<OneD, NekDouble>(physTot, 0.0);

    Vmath::Vcopy(physTot, m_coords[0], 1, out0, 1);
    Vmath::Vcopy(physTot, m_coords[1], 1, out1, 1);
    Vmath::Vcopy(physTot, m_coords[2], 1, out2, 1);
}

void Mapping::v_GetCoordVelocity(
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    for(int i = 0; i < m_nConvectiveFields; ++i)
    {
        outarray[i] = Array<OneD, NekDouble>(physTot, 0.0);
        Vmath::Vcopy(physTot, m_coordsVel[i], 1, outarray[i], 1);
    }
}

void Mapping::v_DotGradJacobian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, NekDouble>                            &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();

    outarray = Array<OneD, NekDouble>(physTot, 0.0);
    if ( !HasConstantJacobian() )
    {
        // Set wavespace to false and store current value
        bool wavespace = m_fields[0]->GetWaveSpace();
        m_fields[0]->SetWaveSpace(false);

        // Get Mapping Jacobian
        Array<OneD, NekDouble> Jac(physTot, 0.0);
        GetJacobian(Jac);

        // Calculate inarray . grad(Jac)
        Array<OneD, NekDouble> wk(physTot, 0.0);
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                    Jac, wk);
            Vmath::Vvtvp(physTot, inarray[i], 1, wk, 1,
                                    outarray, 1, outarray, 1);
        }
        m_fields[0]->SetWaveSpace(wavespace);
    }
}

void Mapping::v_LowerIndex(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    Array<OneD, Array<OneD, NekDouble> > g(nvel*nvel);

    GetMetricTensor(g);

    for (int i=0; i< nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
        for (int j=0; j< nvel; j++)
        {
            Vmath::Vvtvp(physTot, g[i*nvel+j], 1, inarray[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }
    }
}

void Mapping::v_RaiseIndex(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    Array<OneD, Array<OneD, NekDouble> > g(nvel*nvel);

    GetInvMetricTensor(g);

    for (int i=0; i< nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
        for (int j=0; j< nvel; j++)
        {
            Vmath::Vvtvp(physTot, g[i*nvel+j], 1, inarray[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }
    }
}

void Mapping::v_Divergence(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, NekDouble>                            &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> wk(physTot, 0.0);

    Vmath::Zero(physTot, outarray, 1);

    // Set wavespace to false and store current value
    bool wavespace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);

    // Get Mapping Jacobian
    Array<OneD, NekDouble> Jac(physTot, 0.0);
    GetJacobian(Jac);

    for(int i = 0; i < m_nConvectiveFields; ++i)
    {
        Vmath::Vmul(physTot,Jac, 1, inarray[i], 1, wk, 1); // J*Ui
        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i],
                               wk, wk);  // (J*Ui)_i
        Vmath::Vadd(physTot, wk, 1, outarray, 1, outarray, 1);
    }
    Vmath::Vdiv(physTot,outarray,1,Jac,1,outarray,1); //1/J*(J*Ui)_i

    m_fields[0]->SetWaveSpace(wavespace);
}

void Mapping::v_VelocityLaplacian(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const NekDouble                                    alpha)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    Array<OneD, Array<OneD, NekDouble> > tmp(nvel);

    // Set wavespace to false and store current value
    bool wavespace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);

    // Calculate vector gradient wk2 = u^i_(,k) = du^i/dx^k + {i,jk}*u^j
    ApplyChristoffelContravar(inarray, m_wk1);
    for (int i=0; i< nvel; i++)
    {
        if(nvel == 2)
        {
            m_fields[0]->PhysDeriv(inarray[i],
                                m_wk2[i*nvel+0],
                                m_wk2[i*nvel+1]);
        }
        else
        {
            m_fields[0]->PhysDeriv(inarray[i],
                                m_wk2[i*nvel+0],
                                m_wk2[i*nvel+1],
                                m_wk2[i*nvel+2]);
        }
        for (int k=0; k< nvel; k++)
        {
            Vmath::Vadd(physTot,m_wk1[i*nvel+k],1,m_wk2[i*nvel+k],1,
                                                m_wk1[i*nvel+k], 1);
        }
    }
    // Calculate wk1 = A^(ij) = g^(jk)*u^i_(,k)
    for (int i=0; i< nvel; i++)
    {
        for (int k=0; k< nvel; k++)
        {
            tmp[k] = m_wk1[i*nvel+k];
        }
        RaiseIndex(tmp, m_tmp);
        for (int j=0; j<nvel; j++)
        {
            Vmath::Vcopy(physTot, m_tmp[j], 1, m_wk1[i*nvel+j], 1);
        }
    }
    //
    // Calculate L(U)^i = (A^(ij))_(,j) - alpha*d^2(u^i)/dx^jdx^j
    //

    // Step 1 :
    //      d(A^(ij) - alpha*du^i/dx^j)/d(x^j)
    for (int i=0; i< nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble>(physTot,0.0);
        for (int j=0; j< nvel; j++)
        {
            Vmath::Smul(physTot, alpha, m_wk2[i*nvel+j], 1,
                                        m_tmp[0], 1);
            Vmath::Vsub(physTot, m_wk1[i*nvel+j], 1, m_tmp[0], 1,
                                                    m_tmp[0], 1);

            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                    m_tmp[0],
                                    m_tmp[1]);
            Vmath::Vadd(physTot,outarray[i],1,m_tmp[1],1,outarray[i], 1);
        }
    }

    // Step 2: d(A^(ij))/d(x^j) + {j,pj}*A^(ip)
    for (int i=0; i< nvel; i++)
    {
        for (int p=0; p< nvel; p++)
        {
            tmp[p] = m_wk1[i*nvel+p];
        }
        ApplyChristoffelContravar(tmp, m_wk2);
        for (int j=0; j< nvel; j++)
        {
                Vmath::Vadd(physTot,outarray[i],1,m_wk2[j*nvel+j],1,
                                                    outarray[i], 1);
        }
    }

    // Step 3: d(A^(ij))/d(x^j) + {j,pj}*A^(ip) + {i,pj} A^(pj)
    for (int j=0; j< nvel; j++)
    {
        for (int p=0; p< nvel; p++)
        {
            tmp[p] = m_wk1[p*nvel+j];
        }
        ApplyChristoffelContravar(tmp, m_wk2);
        for (int i=0; i< nvel; i++)
        {
                Vmath::Vadd(physTot,outarray[i], 1, m_wk2[i*nvel+j], 1,
                                                        outarray[i], 1);
        }
    }

    // Restore value of wavespace
    m_fields[0]->SetWaveSpace(wavespace);
}

void Mapping::v_gradgradU(
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    // Declare variables
    outarray = Array<OneD, Array<OneD, NekDouble> > (nvel*nvel*nvel);
    Array<OneD, Array<OneD, NekDouble> > tmp(nvel);
    for (int i=0; i< nvel*nvel*nvel; i++)
    {
        outarray[i] = Array<OneD, NekDouble>(physTot,0.0);
    }

    // Set wavespace to false and store current value
    bool wavespace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);

    // Calculate vector gradient u^i_(,j) = du^i/dx^j + {i,pj}*u^p
    ApplyChristoffelContravar(inarray, m_wk1);
    for (int i=0; i< nvel; i++)
    {
        if (nvel == 2)
        {
            m_fields[0]->PhysDeriv(inarray[i],
                                    m_wk2[i*nvel+0],
                                    m_wk2[i*nvel+1]);
        }
        else
        {
            m_fields[0]->PhysDeriv(inarray[i],
                                    m_wk2[i*nvel+0],
                                    m_wk2[i*nvel+1],
                                    m_wk2[i*nvel+2]);
        }
        for (int j=0; j< nvel; j++)
        {
            Vmath::Vadd(physTot,m_wk1[i*nvel+j],1,m_wk2[i*nvel+j],1,
                                                m_wk1[i*nvel+j], 1);
        }
    }

    //
    // Calculate (u^i_,j),k
    //

    // Step 1 : d(u^i_,j))/d(x^k)
    for (int i=0; i< nvel; i++)
    {
        for (int j=0; j< nvel; j++)
        {
            if (nvel == 2)
            {
                m_fields[0]->PhysDeriv(m_wk1[i*nvel+j],
                                        outarray[i*nvel*nvel+j*nvel+0],
                                        outarray[i*nvel*nvel+j*nvel+1]);
            }
            else
            {
                m_fields[0]->PhysDeriv(m_wk1[i*nvel+j],
                                        outarray[i*nvel*nvel+j*nvel+0],
                                        outarray[i*nvel*nvel+j*nvel+1],
                                        outarray[i*nvel*nvel+j*nvel+2]);
            }
        }
    }

    // Step 2: d(u^i_,j)/d(x^k) - {p,jk}*u^i_,p
    for (int i=0; i< nvel; i++)
    {
        for (int p=0; p< nvel; p++)
        {
            tmp[p] = m_wk1[i*nvel+p];
        }
        ApplyChristoffelCovar(tmp, m_wk2);
        for (int j=0; j< nvel; j++)
        {
            for (int k=0; k< nvel; k++)
            {
                Vmath::Vsub(physTot,outarray[i*nvel*nvel+j*nvel+k],1,
                                    m_wk2[j*nvel+k],1,
                                    outarray[i*nvel*nvel+j*nvel+k], 1);
            }
        }
    }

    // Step 3: d(u^i_,j)/d(x^k) - {p,jk}*u^i_,p + {i,pk} u^p_,j
    for (int j=0; j< nvel; j++)
    {
        for (int p=0; p< nvel; p++)
        {
            tmp[p] = m_wk1[p*nvel+j];
        }
        ApplyChristoffelContravar(tmp, m_wk2);
        for (int i=0; i< nvel; i++)
        {
            for (int k=0; k< nvel; k++)
            {
                Vmath::Vadd(physTot,outarray[i*nvel*nvel+j*nvel+k],1,
                                    m_wk2[i*nvel+k],1,
                                    outarray[i*nvel*nvel+j*nvel+k], 1);
            }
        }
    }

    // Restore value of wavespace
    m_fields[0]->SetWaveSpace(wavespace);
}

void Mapping::v_CurlCurlField(
    Array<OneD, Array<OneD, NekDouble> >              &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const bool                                        generalized)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;

    // Set wavespace to false and store current value
    bool wavespace = m_fields[0]->GetWaveSpace();
    m_fields[0]->SetWaveSpace(false);

    // For implicit treatment of viscous terms, we want the generalized
    //   curlcurl and for explicit treatment, we want the cartesian one.
    if (generalized)
    {
        // Get the second derivatives u^i_{,jk}
        Array<OneD, Array<OneD, NekDouble> > tmp(nvel);
        Array<OneD, Array<OneD, NekDouble> > ddU(nvel*nvel*nvel);
        gradgradU(inarray, ddU);

        // Raise index to obtain A^{ip}_{k} = g^pj u^i_{,jk}
        for (int i = 0; i < nvel; ++i)
        {
            for (int k = 0; k < nvel; ++k)
            {
                // Copy to wk
                for (int j = 0; j < nvel; ++j)
                {
                    tmp[j] =  ddU[i*nvel*nvel+j*nvel+k];
                }
                RaiseIndex(tmp, m_tmp);
                for (int p=0; p<nvel; ++p)
                {
                   Vmath::Vcopy(physTot, m_tmp[p], 1,
                                         ddU[i*nvel*nvel+p*nvel+k], 1);
                }
            }
        }
        // The curlcurl is g^ji u^k_{kj} - g^jk u^i_kj = A^{ki}_k - A^{ik}_k
        for (int i = 0; i < nvel; ++i)
        {
            outarray[i] = Array<OneD, NekDouble> (physTot, 0.0);
            for (int k = 0; k < nvel; ++k)
            {
                Vmath::Vadd(physTot, outarray[i], 1,
                                     ddU[k*nvel*nvel+i*nvel+k], 1,
                                     outarray[i], 1);
                Vmath::Vsub(physTot, outarray[i], 1,
                                     ddU[i*nvel*nvel+k*nvel+k], 1,
                                     outarray[i], 1);
            }
        }
    }
    else
    {
        m_fields[0]->CurlCurl(inarray, outarray);
    }

    // Restore value of wavespace
    m_fields[0]->SetWaveSpace(wavespace);
}

void Mapping::v_UpdateBCs( const NekDouble time)
{
    int physTot = m_fields[0]->GetTotPoints();
    int nvel = m_nConvectiveFields;
    int nfields = m_fields.size();
    int nbnds    = m_fields[0]->GetBndConditions().size();

    // Declare variables
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

    Array<OneD, bool>  isDirichlet(nfields);

    Array<OneD, Array<OneD, NekDouble> > values(nfields);
    for (int i=0; i < nfields; i++)
    {
        values[i] = Array<OneD, NekDouble> (physTot, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble> > tmp(nvel);
    Array<OneD, Array<OneD, NekDouble> > tmp2(nvel);
    Array<OneD, Array<OneD, NekDouble> > coordVel(nvel);
    for (int i = 0; i< nvel; i++)
    {
        tmp[i] = Array<OneD, NekDouble> (physTot, 0.0);
        tmp2[i] = Array<OneD, NekDouble> (physTot, 0.0);
        coordVel[i] = Array<OneD, NekDouble> (physTot, 0.0);
    }

    // Get coordinates velocity in transformed system (for MovingBody regions)
    GetCoordVelocity(tmp);
    ContravarFromCartesian(tmp, coordVel);

    // Get  Cartesian coordinates for evaluating boundary conditions
    Array<OneD, Array<OneD, NekDouble> > coords(3);
    for (int dir=0; dir < 3; dir++)
    {
        coords[dir] = Array<OneD, NekDouble> (physTot, 0.0);
    }
    GetCartesianCoordinates(coords[0],coords[1],coords[2]);

    // Loop boundary conditions looking for Dirichlet bc's
    for(int n = 0 ; n < nbnds ; ++n)
    {
        // Evaluate original Dirichlet boundary conditions in whole domain
        for (int i = 0; i < nfields; ++i)
        {
            BndConds   = m_fields[i]->GetBndConditions();
            BndExp     = m_fields[i]->GetBndCondExpansions();
            if ( BndConds[n]->GetBoundaryConditionType() ==
                                SpatialDomains::eDirichlet)
            {
                isDirichlet[i] = true;
                // If we have the a velocity component
                //      check if all vel bc's are also Dirichlet
                if ( i<nvel )
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        ASSERTL0(m_fields[j]->GetBndConditions()[n]->
                                              GetBoundaryConditionType() ==
                                                SpatialDomains::eDirichlet,
                            "Mapping only supported when all velocity components have the same type of boundary conditions");
                    }
                }
                // Check if bc is time-dependent
                ASSERTL0( !BndConds[n]->IsTimeDependent(),
                    "Time-dependent Dirichlet boundary conditions not supported with mapping yet.");

                // Get boundary condition
                LibUtilities::Equation condition =
                    std::static_pointer_cast<
                        SpatialDomains::DirichletBoundaryCondition>
                            (BndConds[n])->
                                m_dirichletCondition;
                // Evaluate
                condition.Evaluate(coords[0], coords[1], coords[2],
                                                time, values[i]);
            }
            else
            {
                isDirichlet[i] = false;
            }
        }
        // Convert velocity vector to transformed system
        if ( isDirichlet[0])
        {
            for (int i = 0; i < nvel; ++i)
            {
                Vmath::Vcopy(physTot, values[i], 1, tmp[i], 1);
            }
            ContravarFromCartesian(tmp, tmp2);
            for (int i = 0; i < nvel; ++i)
            {
                Vmath::Vcopy(physTot, tmp2[i], 1, values[i], 1);
            }
        }

        // Now, project result to boundary
        for (int i = 0; i < nfields; ++i)
        {
            BndConds   = m_fields[i]->GetBndConditions();
            BndExp     = m_fields[i]->GetBndCondExpansions();
            if( BndConds[n]->GetUserDefined() =="" ||
                BndConds[n]->GetUserDefined() =="MovingBody")
            {
                m_fields[i]->ExtractPhysToBnd(n,
                        values[i], BndExp[n]->UpdatePhys());

                // Apply MovingBody correction
                if (  (i<nvel) &&
                      BndConds[n]->GetUserDefined() ==
                      "MovingBody" )
                {
                    // Get coordinate velocity on boundary
                    Array<OneD, NekDouble> coordVelBnd(BndExp[n]->GetTotPoints());
                    m_fields[i]->ExtractPhysToBnd(n, coordVel[i], coordVelBnd);

                    // Apply correction
                    Vmath::Vadd(BndExp[n]->GetTotPoints(),
                                        coordVelBnd, 1,
                                        BndExp[n]->UpdatePhys(), 1,
                                        BndExp[n]->UpdatePhys(), 1);
                }
            }
        }
    }

    // Finally, perform FwdTrans in all fields
    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Get boundary condition information
        BndConds   = m_fields[i]->GetBndConditions();
        BndExp     = m_fields[i]->GetBndCondExpansions();
        for(int n = 0 ; n < BndConds.size(); ++n)
        {
            if ( BndConds[n]->GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet)
            {
                if( BndConds[n]->GetUserDefined() =="" ||
                    BndConds[n]->GetUserDefined() =="MovingBody")
                {
                    if (m_fields[i]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        BndExp[n]->SetWaveSpace(false);
                    }

                    BndExp[n]->FwdTrans_BndConstrained(
                        BndExp[n]->GetPhys(), BndExp[n]->UpdateCoeffs());
                }
            }
        }
    }
}

void Mapping::v_UpdateMapping(
        const NekDouble time,
        const Array<OneD, Array<OneD, NekDouble> > &coords  ,
        const Array<OneD, Array<OneD, NekDouble> > &coordsVel)
{
    if (m_fromFunction)
    {
        std::string s_FieldStr;
        string fieldNames[3] = {"x", "y", "z"};
        string velFieldNames[3] = {"vx", "vy", "vz"};
        // Check if function from session file defines each component
        //      and evaluate them, otherwise there is no need to update
        //          coords
        for(int i = 0; i < 3; i++)
        {
            s_FieldStr = fieldNames[i];
            if ( m_session->DefinesFunction(m_funcName, s_FieldStr))
            {
                EvaluateFunction(m_fields, m_session, s_FieldStr, m_coords[i],
                                        m_funcName, time);
                if ( i==2 && m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    ASSERTL0 (false,
                        "3DH1D does not support mapping in the z-direction.");
                }
            }
            s_FieldStr = velFieldNames[i];
            if ( m_session->DefinesFunction(m_velFuncName, s_FieldStr))
             {
                 EvaluateFunction(m_fields, m_session, s_FieldStr,
                                     m_coordsVel[i], m_velFuncName, time);
                 if ( i==2 && m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                 {
                     ASSERTL0 (false,
                        "3DH1D does not support mapping in the z-direction.");
                 }
             }
        }
    }
    else
    {
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_nConvectiveFields;
        // Copy coordinates
        for(int i = 0; i < 3; i++)
        {
            Vmath::Vcopy(physTot, coords[i], 1, m_coords[i], 1);
        }

        for(int i = 0; i < nvel; i++)
        {
            Vmath::Vcopy(physTot, coordsVel[i], 1, m_coordsVel[i], 1);
        }
    }

    // Update the information required by the specific mapping
    UpdateGeomInfo();
}

}
}
