////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFld.cpp
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
//  Description: FLD file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFld.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputFld::m_className[2] = {
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "fld"), OutputFld::create,
                "Writes a Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "chk"), OutputFld::create,
                "Writes a Fld file."),
        };
        
        OutputFld::OutputFld(FieldSharedPtr f) : OutputModule(f)
        {
        }
        
        OutputFld::~OutputFld()
        {
        }
        
        void OutputFld::Process(po::variables_map &vm)
        {         
            // Extract the output filename and extension
            string filename = m_config["outfile"].as<string>();
            
            if(vm.count("boundary-region"))
            {
                vector<string> tmp1;
                ModuleKey      module;
                
                string boptions = vm["boundary-region"].as<string>().c_str();
                boost::split(tmp1, boptions, boost::is_any_of(":"));
                
                vector<unsigned int> values;
                ASSERTL0(ParseUtils::GenerateOrderedVector(tmp1[0].c_str(),values),
                         "Failed to interpret range string");
                
                if(tmp1.size() == 2)
                {
                    // Extract data to boundaryconditions/ 
                    if(tmp1[1].compare("FldToBoundary") <=0)
                    {
                        for(int i = 0; i < m_f->m_exp.size(); ++i)
                        {
                            m_f->m_exp[i]->FillBndCondFromField();
                        }
                    }
                    
                }
                else
                {
                    ASSERTL0(false,"boundary-region option is not recognised");
                }
                
                if (m_f->m_verbose)
                {
                    cout << "OutputFld: Writing boundary file(s): "; 
                    for(int i = 0; i < values.size(); ++i)
                    {
                        cout << values[i];
                        if(i < values.size()-1) 
                        {
                            cout << ",";
                        }
                    }
                    cout << endl;
                }
                
                int nfields = m_f->m_exp.size();
                Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr> > BndExp(nfields);
                for(int i = 0; i < nfields; ++i)
                {
                    BndExp[i] = m_f->m_exp[i]->GetBndCondExpansions();
                }

                // find ending of output file and insert _b1, _b2
                int    dot = filename.find_last_of('.') + 1;
                string ext = filename.substr(dot, filename.length() - dot);
                string name = filename.substr(0, dot-1);

                for(int i = 0; i < values.size(); ++i)
                {
                    string outname = name  + "_b" + boost::lexical_cast<string>(i) + 
                        "." + ext;
                    
                    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                        = BndExp[0][values[i]]->GetFieldDefinitions();
                    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

                    for(int j = 0; j < nfields; ++j)
                    {
                        for(int k = 0; k < FieldDef.size(); ++k)
                        {
                            BndExp[j][values[i]]->AppendFieldData(FieldDef[k], 
                                                                  FieldData[k]);
                            FieldDef[k]->m_fields.push_back(m_f->m_fielddef[0]->
                                                            m_fields[j]);
                        }
                    }
                    
                    m_f->m_fld->Write(outname,FieldDef,FieldData);


                    // output error for regression checking. 
                    if (vm.count("error"))
                    {
                        int rank = m_f->m_session->GetComm()->GetRank();
                        
                        for (int j = 0; j < nfields; ++j)
                        {
                            BndExp[j][values[i]]->BwdTrans(BndExp[j][values[i]]->GetCoeffs(),BndExp[j][values[i]]->UpdatePhys());
                            NekDouble l2err   = BndExp[j][values[i]]->L2(BndExp[j][values[i]]->GetPhys());
                            NekDouble linferr = BndExp[j][values[i]]->Linf(BndExp[j][values[i]]->GetPhys());
                            
                            if(rank == 0)
                            {
                                cout << "L 2 error (variable "<< FieldDef[0]->m_fields[j] << ") : " << l2err  << endl;
                                cout << "L inf error (variable "<< FieldDef[0]->m_fields[j] << ") : " << linferr << endl;
                            }
                        }
                    }

                }
            }
            else
            {
                if (m_f->m_verbose)
                {
                    cout << "OutputFld: Writing file..." << endl;
                }

                // Write the output file
                m_f->m_fld->Write(filename, m_f->m_fielddef, m_f->m_data);


                // output error for regression checking. 
                if (vm.count("error"))
                {
                    int rank = m_f->m_session->GetComm()->GetRank();
                    
                    for (int j = 0; j < m_f->m_exp.size(); ++j)
                    {
                        if(m_f->m_exp[j]->GetPhysState() == false)
                        {
                            m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(),
                                                    m_f->m_exp[j]->UpdatePhys());
                        }

                        NekDouble l2err   = m_f->m_exp[j]->L2(m_f->m_exp[j]->GetPhys());
                        NekDouble linferr = m_f->m_exp[j]->Linf(m_f->m_exp[j]->GetPhys());
                        if(rank == 0)
                        {
                            cout << "L 2 error (variable "<< m_f->m_fielddef[0]->m_fields[j] << ") : " << l2err  << endl;
                            cout << "L inf error (variable "<< m_f->m_fielddef[0]->m_fields[j] << ") : " << linferr << endl;
                        }
                    }
                }

            }
        }        
    }
}
