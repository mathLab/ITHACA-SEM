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
// Description: Abstract base class for mappings.
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/Mapping.h>

namespace Nektar
{
    namespace GlobalMapping
    {
        MappingSharedPtr Mapping::m_mappingPtr;
        bool             Mapping::m_init      = false;
        bool             Mapping::m_isDefined = false;
        
        MappingFactory& GetMappingFactory()
        {
            typedef Loki::SingletonHolder<MappingFactory,
                                          Loki::CreateUsingNew,
                                          Loki::NoDestroy > Type;
            return Type::Instance();
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
            
            m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(pSession->GetComm());
        }

        /**
         *
         */
        void Mapping::v_InitObject(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const TiXmlElement                                *pMapping)
        {   
            int phystot         = m_fields[0]->GetTotPoints();
            m_timeDependent     = false;
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

            // Calculate information required by the particular mapping
            UpdateGeomInfo();           
        }
      
        /**
         *
         */
        MappingSharedPtr Mapping::Load(
                            const LibUtilities::SessionReaderSharedPtr& pSession,
                            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields)
        {            
            if (!m_init)
            {
                TiXmlElement* vMapping;
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
                            m_fields[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs);
                        }
                    }
                    if (m_timeDependent)
                    {
                        // copy coordinates velocity Data into FieldData and set variable
                        for(int j = 0; j < expdim; ++j)
                        {                
                            m_fields[0]->FwdTrans_IterPerExp(m_coordsVel[j], fieldcoeffs);

                            for(int i = 0; i < FieldDef.size(); ++i)
                            {
                                // Could do a search here to find correct variable
                                FieldDef[i]->m_fields.push_back(velFieldNames[j]);
                                m_fields[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs);
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
            if (pArray.num_elements() != nq)
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
                Vmath::Zero(vCoeffs.num_elements(), vCoeffs, 1);

                LibUtilities::FieldIOSharedPtr fld =
                    MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(m_session->GetComm());
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
        
        void Mapping::ContravarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_ContravarToCartesian( tmp, outarray);
            }
            else
            {
                v_ContravarToCartesian( inarray, outarray);
            }          
        }
        
        void Mapping::CovarToCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_CovarToCartesian( tmp, outarray);
            }
            else
            {
                v_CovarToCartesian( inarray, outarray);
            }          
        }
        
        void Mapping::ContravarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_ContravarFromCartesian( tmp, outarray);
            }
            else
            {
                v_ContravarFromCartesian( inarray, outarray);
            }          
        }
        
        void Mapping::CovarFromCartesian(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_CovarFromCartesian( tmp, outarray);
            }
            else
            {
                v_CovarFromCartesian( inarray, outarray);
            }          
        }
        
        void Mapping::LowerIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_LowerIndex( tmp, outarray);
            }
            else
            {
                v_LowerIndex( inarray, outarray);
            }          
        }
        
        void Mapping::RaiseIndex(
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            if(inarray == outarray)
            {
                int physTot = m_fields[0]->GetTotPoints();
                int nvel = m_nConvectiveFields;
                
                Array<OneD, Array<OneD, NekDouble> >  tmp(nvel);
                for (int i=0; i< nvel; i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot);
                    Vmath::Vcopy(physTot, inarray[i], 1, tmp[i], 1);
                }
                v_RaiseIndex( tmp, outarray);
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
            Array<OneD, NekDouble>               &outarray) 
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
            Array<OneD, Array<OneD, NekDouble> >              &outarray)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;

            Array<OneD, Array<OneD, NekDouble> > wk1(nvel*nvel);
            Array<OneD, Array<OneD, NekDouble> > wk2(nvel*nvel);
            Array<OneD, Array<OneD, NekDouble> > tmp1(nvel);   
            Array<OneD, Array<OneD, NekDouble> > tmp2(nvel);
            for (int i=0; i< nvel; i++)
            {
                tmp1[i] = Array<OneD, NekDouble>(physTot,0.0);
                tmp2[i] = Array<OneD, NekDouble>(physTot,0.0);
                for (int j=0; j< nvel; j++)
                {
                    wk1[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
                    wk2[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
                }
            }

            // Set wavespace to false and store current value
            bool wavespace = m_fields[0]->GetWaveSpace();
            m_fields[0]->SetWaveSpace(false);

            // Calculate vector gradient wk2 = u^i_(,k) = du^i/dx^k + {i,jk}*u^j
            ApplyChristoffelContravar(inarray, wk1);        
            for (int i=0; i< nvel; i++)
            {
                for (int k=0; k< nvel; k++)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],inarray[i],
                                            wk2[i*nvel+k]);

                    Vmath::Vadd(physTot,wk1[i*nvel+k],1,wk2[i*nvel+k],1,wk2[i*nvel+k], 1);               
                }
            }
            // Calculate wk1 = A^(ij) = g^(jk)*u^i_(,k)
            for (int i=0; i< nvel; i++)
            {
                for (int k=0; k< nvel; k++)
                {
                    Vmath::Vcopy(physTot, wk2[i*nvel+k], 1, tmp1[k], 1);
                }
                RaiseIndex(tmp1, tmp2);
                for (int j=0; j<nvel; j++)
                {
                    Vmath::Vcopy(physTot, tmp2[j], 1, wk1[i*nvel+j], 1);
                }
            }
            //
            // Calculate L(U)^i = (A^(ij))_(,j)
            //

            // Step 1 : d(A^(ij))/d(x^j)
            for (int i=0; i< nvel; i++)
            {
                outarray[i] = Array<OneD, NekDouble>(physTot,0.0); 
                for (int j=0; j< nvel; j++)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j],wk1[i*nvel+j],
                                            wk2[0]);  
                    Vmath::Vadd(physTot,outarray[i],1,wk2[0],1,outarray[i], 1);
                }
            }

            // Step 2: d(A^(ij))/d(x^j) + {j,pj}*A^(ip)
            for (int i=0; i< nvel; i++)
            {
                for (int p=0; p< nvel; p++)
                {
                    Vmath::Vcopy(physTot, wk1[i*nvel+p], 1, tmp1[p], 1);
                }
                ApplyChristoffelContravar(tmp1, wk2);
                for (int j=0; j< nvel; j++)
                {
                        Vmath::Vadd(physTot,outarray[i],1,wk2[j*nvel+j],1,outarray[i], 1);                 
                }
            }        

            // Step 3: d(A^(ij))/d(x^j) + {j,pj}*A^(ip) + {i,pj} A^(pj)
            for (int j=0; j< nvel; j++)
            {
                for (int p=0; p< nvel; p++)
                {
                    Vmath::Vcopy(physTot, wk1[p*nvel+j], 1, tmp1[p], 1);
                }
                ApplyChristoffelContravar(tmp1, wk2);
                for (int i=0; i< nvel; i++)
                {
                        Vmath::Vadd(physTot,outarray[i],1,wk2[i*nvel+j],1,outarray[i], 1);               
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
            Array<OneD, Array<OneD, NekDouble> > wk1(nvel*nvel);
            Array<OneD, Array<OneD, NekDouble> > wk2(nvel*nvel);
            Array<OneD, Array<OneD, NekDouble> > tmp(nvel);   
            for (int i=0; i< nvel; i++)
            {
                tmp[i] = Array<OneD, NekDouble>(physTot,0.0);
                for (int j=0; j< nvel; j++)
                {
                    wk1[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
                    wk2[i*nvel+j] = Array<OneD, NekDouble>(physTot,0.0);
                    for (int k=0; k< nvel; k++)
                    {
                        outarray[i*nvel*nvel+j*nvel+k] = Array<OneD, NekDouble>(physTot,0.0);
                    }
                }
            }
            
            // Set wavespace to false and store current value
            bool wavespace = m_fields[0]->GetWaveSpace();
            m_fields[0]->SetWaveSpace(false);
            
            // Calculate vector gradient u^i_(,j) = du^i/dx^j + {i,pj}*u^p
            ApplyChristoffelContravar(inarray, wk1);        
            for (int i=0; i< nvel; i++)
            {
                for (int j=0; j< nvel; j++)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j],inarray[i],
                                            wk2[i*nvel+j]);

                    Vmath::Vadd(physTot,wk1[i*nvel+j],1,wk2[i*nvel+j],1,wk1[i*nvel+j], 1);               
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
                    for (int k=0; k< nvel; k++)
                    {
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                wk1[i*nvel+j], outarray[i*nvel*nvel+j*nvel+k]); 
                    }
                }
            }

            // Step 2: d(u^i_,j)/d(x^k) - {p,jk}*u^i_,p
            for (int i=0; i< nvel; i++)
            {
                for (int p=0; p< nvel; p++)
                {
                    Vmath::Vcopy(physTot, wk1[i*nvel+p], 1, tmp[p], 1);
                }
                ApplyChristoffelCovar(tmp, wk2);
                for (int j=0; j< nvel; j++)
                {
                    for (int k=0; k< nvel; k++)
                    {
                        Vmath::Vsub(physTot,outarray[i*nvel*nvel+j*nvel+k],1,
                                            wk2[j*nvel+k],1,
                                            outarray[i*nvel*nvel+j*nvel+k], 1); 
                    }

                }
            }        

            // Step 3: d(u^i_,j)/d(x^k) - {p,jk}*u^i_,p + {i,pk} u^p_,j
            for (int j=0; j< nvel; j++)
            {
                for (int p=0; p< nvel; p++)
                {
                    Vmath::Vcopy(physTot, wk1[p*nvel+j], 1, tmp[p], 1);
                }
                ApplyChristoffelContravar(tmp, wk2);
                for (int i=0; i< nvel; i++)
                {
                    for (int k=0; k< nvel; k++)
                    {
                        Vmath::Vadd(physTot,outarray[i*nvel*nvel+j*nvel+k],1,
                                            wk2[i*nvel+k],1,
                                            outarray[i*nvel*nvel+j*nvel+k], 1);
                    }                                   
                }
            }            
                       
            // Restore value of wavespace 
            m_fields[0]->SetWaveSpace(wavespace);            
        }
        
        void Mapping::v_CurlCurlField(
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
            Array<OneD, Array<OneD, NekDouble> >              &outarray,
            const bool                                        generalized)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            Array<OneD, NekDouble> tmp (physTot, 0.0);
            Array<OneD, Array<OneD, NekDouble> > wk1(nvel);
            Array<OneD, Array<OneD, NekDouble> > wk2(nvel);
            for (int i = 0; i < nvel; ++i)
            {
                wk1[i] = Array<OneD, NekDouble> (physTot, 0.0);
                wk2[i] = Array<OneD, NekDouble> (physTot, 0.0);
            }
            
            // Set wavespace to false and store current value
            bool wavespace = m_fields[0]->GetWaveSpace();
            m_fields[0]->SetWaveSpace(false);
            
            // For implicit treatment of viscous terms, we want the generalized curlcurl
            //     and for explicit treatment, we want the cartesian one.
            if (generalized)
            {
                // Get the second derivatives u^i_{,jk}
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
                            Vmath::Vcopy(physTot, ddU[i*nvel*nvel+j*nvel+k], 1, wk1[j], 1);
                        }
                        RaiseIndex(wk1, wk2);
                        for (int p=0; p<nvel; ++p)
                        {
                           Vmath::Vcopy(physTot, wk2[p], 1, ddU[i*nvel*nvel+p*nvel+k], 1);
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
                switch(nvel)
                {
                    case 2:
                    {
                        Array<OneD,NekDouble> Vx(physTot);
                        Array<OneD,NekDouble> Uy(physTot);
                        Array<OneD,NekDouble> Dummy(physTot);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], inarray[1], Vx);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], inarray[0], Uy);

                        Vmath::Vsub(physTot, Vx, 1, Uy, 1, Dummy, 1);

                        m_fields[0]->PhysDeriv(Dummy,outarray[1],outarray[0]);

                        Vmath::Smul(physTot, -1.0, outarray[1], 1, outarray[1], 1);
                    }
                    break;
                    
                    case 3:
                    {
                        // Declare variables
                        Array<OneD,NekDouble> Ux(physTot);
                        Array<OneD,NekDouble> Uy(physTot);
                        Array<OneD,NekDouble> Uz(physTot);

                        Array<OneD,NekDouble> Vx(physTot);
                        Array<OneD,NekDouble> Vy(physTot);
                        Array<OneD,NekDouble> Vz(physTot);

                        Array<OneD,NekDouble> Wx(physTot);
                        Array<OneD,NekDouble> Wy(physTot);
                        Array<OneD,NekDouble> Wz(physTot);

                        Array<OneD,NekDouble> Dummy1(physTot);
                        Array<OneD,NekDouble> Dummy2(physTot);

                        // Calculate gradient
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], inarray[0], Ux);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], inarray[0], Uy);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], inarray[0], Uz);        

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], inarray[1], Vx);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], inarray[1], Vy);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], inarray[1], Vz);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], inarray[2], Wx);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], inarray[2], Wy);
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], inarray[2], Wz);

                        // x-component
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Vx, Dummy1); //Vxy
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Uy, Dummy2); //Uyy
                        Vmath::Vsub(physTot, Dummy1, 1, Dummy2, 1, outarray[0],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], Uz, Dummy1); //Uzz
                        Vmath::Vsub(physTot, outarray[0],   1, Dummy1, 1,  outarray[0],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Wz, Dummy2); //Wxz
                        Vmath::Vadd(physTot, outarray[0],   1, Dummy2, 1,  outarray[0],   1);

                        // y-component
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Wz, Dummy1); //Wzy
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2], Vz, Dummy2); //Vzz       
                        Vmath::Vsub(physTot, Dummy1, 1, Dummy2, 1, outarray[1],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Vx, Dummy1); //Vxx
                        Vmath::Vsub(physTot, outarray[1],   1, Dummy1, 1,  outarray[1],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Uy, Dummy2); //Uyx
                        Vmath::Vadd(physTot, outarray[1],   1, Dummy2, 1,  outarray[1],   1);  

                        // z-component
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Uz, Dummy1); //Uxz
                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Wy, Dummy2); //Wyy       
                        Vmath::Vsub(physTot, Dummy1, 1, Dummy2, 1, outarray[2],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0], Wx, Dummy1); //Wxx
                        Vmath::Vsub(physTot, outarray[2],   1, Dummy1, 1,  outarray[2],   1);

                        m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1], Vz, Dummy2); //Vyz
                        Vmath::Vadd(physTot, outarray[2],   1, Dummy2, 1,  outarray[2],   1);                        
                    }
                    break;
                }
            }

            // Restore value of wavespace 
            m_fields[0]->SetWaveSpace(wavespace);            
        }
        
        void Mapping::v_UpdateBCs( const NekDouble time)
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            int nfields = m_fields.num_elements();
            int nbnds    = m_fields[0]->GetBndConditions().num_elements();

            // Declare variables
            Array<OneD, int> BCtoElmtID;
            Array<OneD, int> BCtoTraceID;
            Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
            Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
            StdRegions::StdExpansionSharedPtr elmt;
            StdRegions::StdExpansionSharedPtr Bc;

            Array<OneD, NekDouble>  ElmtVal(physTot, 0.0);
            Array<OneD, NekDouble>  BndVal(physTot, 0.0);
            Array<OneD, NekDouble>  coordVelElmt(physTot, 0.0);
            Array<OneD, NekDouble>  coordVelBnd(physTot, 0.0);
            Array<OneD, NekDouble>  Vals(physTot, 0.0);
            
            Array<OneD, bool>  isDirichlet(nfields);
            
            
            Array<OneD, Array<OneD, NekDouble> > values(nfields);
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
                values[dir] = Array<OneD, NekDouble> (physTot, 0.0);
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
                                ASSERTL0(m_fields[j]->GetBndConditions()[n]->GetBoundaryConditionType() == 
                                                        SpatialDomains::eDirichlet,
                                    "Mapping only supported when all velocity components have the same type of boundary conditions");
                            }                            
                        }
                        // Check if bc is time-dependent
                        ASSERTL0( !BndConds[n]->IsTimeDependent(),
                            "Time-dependent Dirichlet boundary conditions not supported with mapping yet.");

                        // Get boundary condition 
                       LibUtilities::Equation condition =
                            boost::static_pointer_cast<
                                SpatialDomains::DirichletBoundaryCondition>
                                    (BndConds[n])->
                                        m_dirichletCondition;
                        // Evaluate
                        condition.Evaluate(coords[0], coords[1], coords[2],
                                                        0.0, values[i]);
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
                    
                    switch (m_fields[0]->GetExpType())
                    {
                        case MultiRegions::e2D:
                        case MultiRegions::e3D:
                        {                            
                            // Loop boundary conditions again to get correct
                            //    values for cnt
                            int cnt = 0;
                            for(int m = 0 ; m < nbnds; ++m)
                            {
                                int exp_size = BndExp[m]->GetExpSize();
                                if (m==n && isDirichlet[i])
                                {
                                    for (int j = 0; j < exp_size; ++j, cnt++)
                                    {                        
                                        m_fields[i]->GetBoundaryToElmtMap(BCtoElmtID,BCtoTraceID);
                                        /// Casting the boundary expansion to the specific case
                                        Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                                                (BndExp[n]->GetExp(j));
                                        // Get element expansion
                                        elmt = m_fields[i]->GetExp(BCtoElmtID[cnt]);
                                        // Get values on the element
                                        ElmtVal  = values[i] + m_fields[i]->GetPhys_Offset(
                                                                             BCtoElmtID[cnt]);                                        
                                        // Get values on boundary
                                        switch (m_fields[i]->GetExpType())
                                        {
                                            case MultiRegions::e2D:
                                            {
                                                elmt->GetEdgePhysVals(BCtoTraceID[cnt], Bc,
                                                                      ElmtVal, BndVal);
                                            }
                                            break;

                                            case MultiRegions::e3D:
                                            {
                                                elmt->GetFacePhysVals(BCtoTraceID[cnt], Bc, 
                                                                      ElmtVal, BndVal);
                                            }
                                            break;

                                            default:
                                                ASSERTL0(0,"Dimension not supported");
                                            break;
                                        }
                                        // Pointer to value that should be updated
                                        Vals = BndExp[n]->UpdatePhys()
                                                    + BndExp[n]->GetPhys_Offset(j);                            

                                        // Copy result
                                        Vmath::Vcopy(Bc->GetTotPoints(), BndVal, 1, Vals, 1);
                                        
                                        // Apply MovingBody correction
                                        if (  (i<nvel) && BndConds[n]->GetUserDefined() == "MovingBody" )
                                        {
                                            // get coordVel in the element
                                            coordVelElmt  = coordVel[i] + m_fields[i]->GetPhys_Offset(
                                                                BCtoElmtID[cnt]);
                                            
                                            // Get values on boundary
                                            switch (m_fields[i]->GetExpType())
                                            {
                                                case MultiRegions::e2D:
                                                {
                                                    elmt->GetEdgePhysVals(BCtoTraceID[cnt], Bc,
                                                        coordVelElmt, coordVelBnd);
     
                                                }
                                                break;

                                                case MultiRegions::e3D:
                                                {
                                                    elmt->GetFacePhysVals(BCtoTraceID[cnt], Bc, 
                                                        coordVelElmt, coordVelBnd);
                                                }
                                                break;

                                                default:
                                                    ASSERTL0(0,"Dimension not supported");
                                                break;
                                            }
                                            // Apply correction
                                            Vmath::Vadd(Bc->GetTotPoints(), coordVelBnd, 1, 
                                                                    Vals, 1, Vals, 1);
                                        }                       
                                    }   
                                }
                                else // setting if m!=n
                                {
                                    cnt += exp_size;
                                }
                            }                            
                        }
                        break;

                        case MultiRegions::e3DH1D:
                        {
                            Array<OneD, unsigned int> planes;
                            planes = m_fields[0]->GetZIDs();
                            int num_planes = planes.num_elements();            
                            // Loop boundary conditions again to get correct
                            //    values for cnt
                            int cnt = 0;
                            for(int k = 0; k < num_planes; k++)
                            {
                                for(int m = 0 ; m < nbnds; ++m)
                                {
                                    int exp_size = BndExp[m]->GetExpSize();
                                    int exp_size_per_plane = exp_size/num_planes;
                                    if (m==n && isDirichlet[i])
                                    {
                                        for (int j = 0; j < exp_size_per_plane; ++j, cnt++)
                                        {                        
                                            int bndElmtOffset = j+k*exp_size_per_plane;

                                            BndConds   = m_fields[i]->GetBndConditions();
                                            BndExp     = m_fields[i]->GetBndCondExpansions();
                                            m_fields[i]->GetBoundaryToElmtMap(BCtoElmtID,BCtoTraceID);
                                            /// Casting the boundary expansion to the specific case
                                            Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                                                    (BndExp[n]->GetExp(bndElmtOffset));
                                            // Get element expansion
                                            elmt = m_fields[i]->GetExp(BCtoElmtID[cnt]);
                                            // Get velocity on the element
                                            ElmtVal  = values[i] + m_fields[i]->GetPhys_Offset(
                                                                                 BCtoElmtID[cnt]);

                                            // Get values on boundary
                                            elmt->GetEdgePhysVals(BCtoTraceID[cnt], Bc,
                                                                  ElmtVal, BndVal);
    
                                            // Pointer to value that should be updated
                                            Vals = BndExp[n]->UpdatePhys()
                                                        + BndExp[n]->GetPhys_Offset(bndElmtOffset);                            

                                            // Copy result
                                            Vmath::Vcopy(Bc->GetTotPoints(), BndVal, 1, Vals, 1);
                                            
                                            // Apply coordinate velocity correction
                                            if (  (i<nvel) && BndConds[n]->GetUserDefined() == "MovingBody" )
                                            {
                                                // Get coordinate velocity on the element
                                                coordVelElmt  = coordVel[i] + m_fields[i]->GetPhys_Offset(
                                                                        BCtoElmtID[cnt]);
                                                // Get values on boundary
                                                elmt->GetEdgePhysVals(BCtoTraceID[cnt], Bc,
                                                                  coordVelElmt, coordVelBnd);
                                                // Apply correction
                                                Vmath::Vadd(Bc->GetTotPoints(), coordVelBnd, 1, 
                                                                        Vals, 1, Vals, 1);
                                            }                       
                                        }   
                                    }
                                    else // setting if m!=n
                                    {
                                        cnt += exp_size_per_plane;
                                    }
                                }                                
                            }                            
                        }
                        break;

                        default:
                            ASSERTL0(0,"Dimension not supported");
                        break;                        
                    }                    
                }                
            }
            
            // Finally, perform FwdTrans in all fields
            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                // Get boundary condition information
                BndConds   = m_fields[i]->GetBndConditions();
                BndExp     = m_fields[i]->GetBndCondExpansions();                                
                for(int n = 0 ; n < BndConds.num_elements(); ++n)
                {   
                    if ( BndConds[n]->GetBoundaryConditionType() == 
                            SpatialDomains::eDirichlet)
                    {
                        BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),
                                                    BndExp[n]->UpdateCoeffs());
                        if (m_fields[i]->GetExpType() == MultiRegions::e3DH1D)
                        {
                            BndExp[n]->HomogeneousFwdTrans(BndExp[n]->GetCoeffs(),
                                                    BndExp[n]->UpdateCoeffs());
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
