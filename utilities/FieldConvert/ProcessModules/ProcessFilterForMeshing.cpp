////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessQCriterion.cpp
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
//  Description: Computes Q Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessFilterForMeshing.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessFilterForMeshing::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "filterformeshing"),
        ProcessFilterForMeshing::create, "Extracts points for use in meshing");

ProcessFilterForMeshing::ProcessFilterForMeshing(FieldSharedPtr f)
    : ProcessModule(f)
{
}

ProcessFilterForMeshing::~ProcessFilterForMeshing()
{
}

void ProcessFilterForMeshing::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessFilterForMeshing: obtaining points for adaptive meshing..." << endl;
    }

    int ct = 0;

    ASSERTL0(m_f->m_exp.size() == 3,"not correct number of fields");

    Array<OneD, NekDouble> U = m_f->m_exp[0]->UpdatePhys();
    Array<OneD, NekDouble> V = m_f->m_exp[1]->UpdatePhys();
    Array<OneD, NekDouble> W = m_f->m_exp[2]->UpdatePhys();

    Array<OneD, NekDouble> x(m_f->m_exp[0]->GetTotPoints()),
                           y(m_f->m_exp[0]->GetTotPoints()),
                           z(m_f->m_exp[0]->GetTotPoints());
    m_f->m_exp[0]->GetCoords(x,y,z);

    ofstream file;
    file.open("refinment_pts.txt");

    for(int i = 0; i < m_f->m_exp[0]->GetNumElmts(); i++)
    {
        int offset = m_f->m_exp[0]->GetPhys_Offset(i);
        LocalRegions::ExpansionSharedPtr el = m_f->m_exp[0]->GetExp(i);
        if(el->GetGeom()->GetShapeType() == LibUtilities::ePrism)
        {
            continue;
        }

        for(int j = 0; j < el->GetTotPoints(); j++)
        {
            if(fabs(U[j+offset]*U[j+offset] +
                    V[j+offset]*V[j+offset] +
                    W[j+offset]*W[j+offset] - 2.0) < 1e-2)
            {
                file << x[j+offset] << " " << y[j+offset] << " " << z[j+offset] << endl;
                ct++;
            }
        }
    }

    file.close();

    cout << ct << " points" << endl;

    exit(-1);


}

}
}
