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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputVtk.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey OutputVtk::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eOutputModule, "vtu"), OutputVtk::create,
                "Writes a VTU file.");

        OutputVtk::OutputVtk(FieldSharedPtr f) : OutputModule(f)
        {
            
        }

        OutputVtk::~OutputVtk()
        {

        }
        
        void OutputVtk::Process()
        {
            int i, j;
            if (f->verbose)
            {
                cout << "OutputVtk: Writing file..." << endl;
            }

            // Write solution.
            string fname = std::string("flate-plate.vtu");
            ofstream outfile(fname.c_str());
            f->exp[0]->WriteVtkHeader(outfile);
            
            
            
            // For each field write out field data for each expansion.
            for (i = 0; i < f->exp[0]->GetNumElmts(); ++i)
            {
                f->exp[0]->WriteVtkPieceHeader(outfile,i);
                // For this expansion write out each field.
                for (j = 0; j < f->fielddef[0]->m_fields.size(); ++j)
                {
                    f->exp[j]->WriteVtkPieceData(outfile, i, 
                                                 f->fielddef[0]->m_fields[j]);
                }
                f->exp[0]->WriteVtkPieceFooter(outfile, i);
            }
            f->exp[0]->WriteVtkFooter(outfile);
            cout << "Written file: " << fname << endl;
        }        
    }
}
