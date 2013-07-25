////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.cpp
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
//  Description: DAT file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputTecplot.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputTecplot::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "dat"), OutputTecplot::create,
                "Writes a Tecplot file.");

        OutputTecplot::OutputTecplot(FieldSharedPtr f) : OutputModule(f)
        {
        }

        OutputTecplot::~OutputTecplot()
        {
        }
        
        void OutputTecplot::Process()
        {
            int i, j;
            if (f->verbose)
            {
                cout << "OutputTecplot: Writing file..." << endl;
            }
            
            int expdim  = f->graph->GetMeshDimension();

            // Extract the output filename and extension
            string filename = config["outfile"].as<string>();
            
            // Write solution.
            ofstream outfile(filename.c_str());
    
            std::string var = "";
                
            for (int j = 0; j < f->fielddef[0]->m_fields.size(); ++j)
            {
                var = var + ", " + f->fielddef[0]->m_fields[j];
            }
                
            if (expdim == 3)
            {
                f->exp[0]->WriteTecplotHeader(outfile,var);
                f->exp[0]->WriteTecplotZone(outfile);
                for(int j = 0; j < f->fielddef[0]->m_fields.size(); ++j)
                {
                    f->exp[j]->WriteTecplotField(outfile);
                }
                    f->exp[0]->WriteTecplotConnectivity(outfile);
            }
            else
            {
                f->exp[0]->WriteTecplotHeader(outfile,var);
                for (int i = 0; i < f->exp[0]->GetNumElmts(); ++i)
                {
                    f->exp[0]->WriteTecplotZone(outfile,i);
                    for (int j = 0; j < f->fielddef[0]->m_fields.size(); ++j)
                    {
                        f->exp[j]->WriteTecplotField(outfile,i);
                    }
                }
            }
            
            cout << "Written file: " << filename << endl;
        }        
    }
}
