////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessVorticity.cpp
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
//  Description: Computes vorticity field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessVorticity.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessVorticity::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "vorticity"), 
                ProcessVorticity::create, "Computes vorticity field.");

        ProcessVorticity::ProcessVorticity(FieldSharedPtr f) : ProcessModule(f)
        {
        }

        ProcessVorticity::~ProcessVorticity()
        {
        }

        void ProcessVorticity::Process(po::variables_map &vm)
        {
            if (m_f->m_verbose)
            {
                cout << "ProcessVorticity: Calculating vorticity..." << endl;
            }
            
            int i, j;
            int expdim    = m_f->m_graph->GetMeshDimension();
            int spacedim  = expdim;
            if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
                (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
            {
                spacedim = 3;
            }
            int nfields = m_f->m_fielddef[0]->m_fields.size();
            if (spacedim == 1)
            {
                ASSERTL0(false, "Error: Vorticity for a 1D problem cannot "
                                "be computed")
            }
            int addfields = (spacedim == 2)? 1:3;
            
            int npoints = m_f->m_exp[0]->GetNpoints();
            Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
            Array<OneD, Array<OneD, NekDouble> > outfield(addfields);
            m_f->m_exp.resize(nfields+addfields);

            
            for (i = 0; i < nfields*nfields; ++i)
            {
                grad[i] = Array<OneD, NekDouble>(npoints);
            }
            
            for (i = 0; i < addfields; ++i)
            {
                outfield[i] = Array<OneD, NekDouble>(npoints);
            }
            
            // Calculate Gradient & Vorticity
            if (spacedim == 2)
            {
                for (i = 0; i < nfields; ++i)
                {
                    m_f->m_exp[i]->PhysDeriv(m_f->m_exp[i]->GetPhys(), 
                                             grad[i*nfields], 
                                             grad[i*nfields+1]);
                }
                // W_z = Vx - Uy
                Vmath::Vsub(npoints, grad[1*nfields+0], 1, 
                            grad[0*nfields+1], 1, 
                            outfield[0], 1);
            }
            else
            {
                for (i = 0; i < nfields; ++i)
                {

                    m_f->m_exp[i]->PhysDeriv(m_f->m_exp[i]->GetPhys(), 
                                             grad[i*nfields], 
                                             grad[i*nfields+1],
                                             grad[i*nfields+2]);
                }
                
                // W_x = Wy - Vz
                Vmath::Vsub(npoints, grad[2*nfields+1], 1, grad[1*nfields+2], 1, 
                            outfield[0],1);
                // W_y = Uz - Wx
                Vmath::Vsub(npoints, grad[0*nfields+2], 1, grad[2*nfields+0], 1, 
                            outfield[1], 1);
                // W_z = Vx - Uy
                Vmath::Vsub(npoints, grad[1*nfields+0], 1, grad[0*nfields+1], 1, 
                            outfield[2], 1);
            }

            for (i = 0; i < addfields; ++i)
            {
                m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
                m_f->m_exp[nfields + i]->UpdatePhys() = outfield[i];
                m_f->m_exp[nfields + i]->FwdTrans_IterPerExp(outfield[i],
                                    m_f->m_exp[nfields + i]->UpdateCoeffs());
            }
            
            vector<string > outname;
            if (addfields == 1)
            {
                outname.push_back("W_z");
            }
            else
            {
                outname.push_back("W_x");
                outname.push_back("W_y");
                outname.push_back("W_z");
            }
            
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = m_f->m_exp[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
            
            for (j = 0; j < nfields + addfields; ++j)
            {
                for (i = 0; i < FieldDef.size(); ++i)
                {   
                    if (j >= nfields)
                    {
                        FieldDef[i]->m_fields.push_back(outname[j-nfields]);
                    }
                    else
                    {
                        FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
                    }
                    m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
                }
            }
            
            m_f->m_fielddef = FieldDef;
            m_f->m_data     = FieldData;


        }
    }
}
