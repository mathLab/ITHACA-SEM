////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJacobianEnergy.cpp
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
//  Description: Compute energy of Jacobian.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessJacobianEnergy.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessJacobianEnergy::className =
        GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "jacobianenergy"),
                ProcessJacobianEnergy::create,
                "Show high frequency energy of Jacobian.");

ProcessJacobianEnergy::ProcessJacobianEnergy(FieldSharedPtr f) :
    ProcessModule(f)
{
    m_config["topmodes"] = ConfigOption(false, "1",
                                        "how many top modes to keep ");
}

ProcessJacobianEnergy::~ProcessJacobianEnergy()
{
}

void ProcessJacobianEnergy::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessJacobianEnergy: Process Jacobian fld" << endl;
    }

    Array<OneD, NekDouble> phys   = m_f->m_exp[0]->UpdatePhys();
    Array<OneD, NekDouble> coeffs = m_f->m_exp[0]->UpdateCoeffs();
    Array<OneD, NekDouble> tmp,tmp1;

    for(int i =0; i < m_f->m_exp[0]->GetExpSize(); ++i)
    {
        // copy Jacobian into field
        StdRegions::StdExpansionSharedPtr Elmt = m_f->m_exp[0]->GetExp(i);

        int ncoeffs = Elmt->GetNcoeffs();
        int nquad   = Elmt->GetTotPoints();
        int coeffoffset = m_f->m_exp[0]->GetCoeff_Offset(i);
        Array<OneD, NekDouble> coeffs1(ncoeffs);
        Array<OneD, const NekDouble> Jac =
                        Elmt->GetMetricInfo()->GetJac(Elmt->GetPointsKeys());
        if(Elmt->GetMetricInfo()->GetGtype() == SpatialDomains::eRegular)
        {
            Vmath::Fill(nquad,Jac[0],phys,1);
        }
        else
        {
            Vmath::Vcopy(nquad,Jac,1,phys,1);
        }

        if(Elmt->GetMetricInfo()->GetGtype() == SpatialDomains::eDeformed)
        {
            NekDouble jacmax = Vmath::Vmax(nquad,Jac,1);
            NekDouble jacmin = Vmath::Vmin(nquad,Jac,1);

            NekDouble jacmeasure = jacmax/jacmin -1.0;
            Vmath::Fill(nquad,jacmeasure,phys,1);
        }
        else
        {
            Vmath::Fill(nquad,0.0,phys,1);
        }

        Elmt->FwdTrans(phys,tmp = coeffs + coeffoffset);
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (int i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("JacobianEnergy");
        m_f->m_exp[0]->AppendFieldData(FieldDef[i], FieldData[i]);
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}
