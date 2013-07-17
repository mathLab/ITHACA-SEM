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
//  Description: Refine prismatic boundary layer elements.
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
            /*
            // Vorticity configuration.
            config["layers"]     = ConfigOption(false, "2",
                "Number of layers to refine.");
            config["npoints"]         = ConfigOption(false, "5",
                "Number of points in high order elements.");
            config["surf"]       = ConfigOption(false, "",
                "Tag identifying surface connected to prism.");
            config["r"]          = ConfigOption(false, "2.0",
                "Ratio to use in geometry progression.");
             */
        }

        ProcessVorticity::~ProcessVorticity()
        {

        }

        void ProcessVorticity::Process()
        {
            if (f->verbose)
            {
                cout << "ProcessVorticity: Calculating vorticity..." << endl;
            }
            
            int i, j;
            int expdim  = f->graph->GetMeshDimension();
            int spacedim  = f->graph->GetSpaceDimension();
            int nfields = f->fielddef[0]->m_fields.size();
            if (spacedim == 1)
            {
                ASSERTL0(false, "Error: Vorticity for a 1D problem cannot "
                                "be computed")
            }
            int addfields = (spacedim == 2)? 1:3;
            
            cout << addfields << endl;
            
            int npoints = f->exp[0]->GetNpoints();
            int ncoeffs = f->exp[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
            Array<OneD, Array<OneD, NekDouble> > outfield(addfields);
            f->exp.resize(nfields+addfields);

            
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
                cout << "ciao ciao "<<endl;

                for (i = 1; i < nfields; ++i)
                {
                    f->exp[i]->PhysDeriv(f->exp[i]->GetPhys(), grad[i*nfields], 
                                         grad[i*nfields+1]);
                }
                // W_z = Vx - Uy
                Vmath::Vsub(npoints, grad[1*nfields+0], 1, grad[0*nfields+1], 1, 
                            outfield[0], 1);
                cout << "ciao ciao "<<endl;

            }
            else
            {
                for (i = 1; i < nfields; ++i)
                {

                    f->exp[i]->PhysDeriv(f->exp[i]->GetPhys(), grad[i*nfields], 
                                         grad[i*nfields+1], grad[i*nfields+2]);
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
                cout << "ciao ciao"<<endl;
                Array<OneD, NekDouble>  tmp(ncoeffs);
                
                //f->exp[nfields + i] = f->AppendExpList();
                //f->AppendExpList("W_z");
                f->exp[nfields + i]->FwdTrans(outfield[i], f->exp[nfields + i]->UpdateCoeffs());
                
                cout << "ciao ciaoooo"<<endl;

            }

            /*
            string out(argv[argc-1]);
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = f->exp[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
            */
            
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
            /*
            for (j = 0; j < nfields + addfields; ++j)
            {
                for (i = 0; i < f->fielddef.size(); ++i)
                {   
                    if (j >= nfields)
                    {
                        f->fielddef[i]->m_fields.push_back(outname[j-nfields]);
                    }
                    else
                    {
                        f->fielddef[i]->m_fields.push_back(f->fielddef[i]->
                                                           m_fields[j]);
                    }
                    f->exp[j]->AppendFieldData(f->fielddef[i], f->data[i]);
                }
            }
            LibUtilities::Write("Flate-Plate-coarse.fld", f->fielddef, f->data);
             */
        }
    }
}
