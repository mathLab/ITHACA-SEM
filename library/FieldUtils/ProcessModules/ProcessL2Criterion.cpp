////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessL2Criterion.cpp
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
//  Description: Computes Lambda 2 Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessL2Criterion.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

namespace Nektar
{
    namespace FieldUtils
    {

        ModuleKey ProcessL2Criterion::className =
                GetModuleFactory().RegisterCreatorFunction(
                        ModuleKey(eProcessModule, "L2Criterion"),
                        ProcessL2Criterion::create,
                        "Computes Lambda 2 Criterion.");

        ProcessL2Criterion::ProcessL2Criterion(FieldSharedPtr f) : ProcessModule(f)
        {
        }

        ProcessL2Criterion::~ProcessL2Criterion()
        {
        }

        //Eigen values of a symmetric 3x3 matrix
        void MatSymEVals(NekDouble d1, NekDouble d2, NekDouble d3,
                         NekDouble a, NekDouble b, NekDouble c,
                         NekDouble& l1, NekDouble& l2, NekDouble& l3)
        {
            NekDouble p = a*a+b*b+c*c;
            if(p == 0)
            {
                l1=d1; l2=d2; l3=d3;
                if (l1 > l3)
                    swap(l1, l3);
                if (l1 > l2)
                    swap(l1, l2);
                if (l2 > l3)
                    swap(l2, l3);
                //cout << "Diagonal ";
            } else {
                NekDouble q = (d1 + d2 + d3)/3.0;
                p = (d1 - q) * (d1 - q) + (d2 - q) * (d2 - q) + (d3 - q) * (d3 - q) + 2.0 * p;
                p = sqrt(p / 6.0);
                NekDouble r = -0.5 * (a * a * d3 - a * a * q - 2.0 * a * b * c + b * b * d2 - b * b * q + c * c * d1
                                      - c * c * q - d1 * d2 * d3 + d1 * d2 * q + d1 * d3 * q - d1 * q * q + d2 * d3 * q
                                      - d2 * q * q - d3 * q * q + q * q * q) / (p * p * p);

                NekDouble phi = 0;
                if (r <= -1)
                    phi = PI / 3.0;
                else if (r >= 1)
                    phi = 0.0;
                else
                    phi = acos(r) / 3.0;

                // the eigenvalues satisfy eig3 >= eig2 >= eig1
                l3 = q + 2.0 * p * cos(phi);
                l1 = q + 2.0 * p * cos(phi + (2.0 * PI / 3.0));
                l2 = 3.0 * q - l1 - l3;   // since trace(A) = eig1 + eig2 + eig3
            }
            //cout << l1 << " " << l2 << " " << l3 << endl;
        }

        void ProcessL2Criterion::Process(po::variables_map &vm)
        {
            auto nfields = m_f->m_variables.size();
            //m_f->m_variables.push_back("L1");
            m_f->m_variables.push_back("L2");
            //m_f->m_variables.push_back("L3");

            // Skip in case of empty partition
            if (m_f->m_exp[0]->GetNumElmts() == 0)
            {
                return;
            }

            int i, s;
            int expdim   = m_f->m_graph->GetMeshDimension();
            int spacedim = expdim + (m_f->m_numHomogeneousDir);

            ASSERTL0(spacedim == 3,
                     "ProcessL2Criterion must be computed for a 3D (or quasi-3D) case.");

            int npoints = m_f->m_exp[0]->GetNpoints();

            Array<OneD, Array<OneD, NekDouble> > grad(spacedim * spacedim);

            // Will store the Lambdas
            NekDouble a00, a11, a22, a01, a02, a12;
            NekDouble t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
            //Array<OneD, NekDouble> outfield1 (npoints);
            NekDouble outfield1, outfield3;
            Array<OneD, NekDouble> outfield2 (npoints);
            //Array<OneD, NekDouble> outfield3 (npoints);

            int nstrips;
            m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

            for (i = 0; i < spacedim * spacedim; ++i)
            {
                grad[i] = Array<OneD, NekDouble>(npoints);
            }

            MultiRegions::ExpListSharedPtr Exp;

            for (s = 0; s < nstrips; ++s) // homogeneous strip varient
            {
                for (i = 0; i < spacedim; ++i)
                {
                    m_f->m_exp[s * nfields + i]->PhysDeriv(
                            m_f->m_exp[s * nfields + i]->GetPhys(), grad[i * spacedim],
                            grad[i * spacedim + 1], grad[i * spacedim + 2]);
                }

                for(int j=0; j<npoints; ++j)
                {
                    t1 = grad[0 * spacedim + 1][j] + grad[1 * spacedim + 0][j];//u01 + u10;
                    t2 = grad[0 * spacedim + 2][j] + grad[2 * spacedim + 0][j];//u02 + u20;
                    t3 = grad[0 * spacedim + 1][j] - grad[1 * spacedim + 0][j];//u01 - u10;
                    t4 = grad[0 * spacedim + 2][j] - grad[2 * spacedim + 0][j];//u02 - u20;
                    t5 = t2*t2;
                    t6 = t4*t4;
                    t7 = t3*t3;
                    t8 = t1*t1;
                    t9 = 0.25;
                    t10 = grad[2 * spacedim + 1][j] + grad[1 * spacedim + 2][j];//u21 + u12;
                    t11 = grad[2 * spacedim + 1][j] - grad[1 * spacedim + 2][j];//u21 - u12;
                    t12 = 0.5;
                    t13 = t9 * (t10 * t2 + t11 * t4) + t12 * t1 * (grad[0 * spacedim + 0][j]+grad[1 * spacedim + 1][j]);//(u00 + u11);
                    t14 = t12 * t2 * (grad[0 * spacedim + 0][j]+grad[2 * spacedim + 2][j]) + t9 * (t1 * t10 - t11 * t3);//(u00 + u22) + t9 * (t1 * t10 - t11 * t3);
                    t15 = t10*t10;
                    t11 = t11*t11;
                    t1 = t12 * t10 * (grad[1 * spacedim + 1][j]+grad[2 * spacedim + 2][j]) - t9 * (-t1 * t2 + t3 * t4);

                    a00 = t9 * (t5 - t6 - t7 + t8) + grad[0 * spacedim + 0][j] * grad[0 * spacedim + 0][j];
                    a01 = t13;
                    a02 = t14;
                    a11 = t9 * (-t7 + t8 + t15 - t11) + grad[1 * spacedim + 1][j] * grad[1 * spacedim + 1][j];
                    a12 = t1;
                    a22 = t9 * (t5 - t6 + t15 - t11) + grad[2 * spacedim + 2][j] * grad[2 * spacedim + 2][j];

                    MatSymEVals(a00, a11, a22, a01, a02, a12,
                                outfield1,  outfield2[j],  outfield3);
                }

                Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
                Vmath::Vcopy(npoints, outfield2, 1, Exp->UpdatePhys(), 1);
                Exp->FwdTrans_IterPerExp(outfield2, Exp->UpdateCoeffs());
                auto it = m_f->m_exp.begin() + s * (nfields + 1) + nfields;
                m_f->m_exp.insert(it, Exp);

                /*Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
                Vmath::Vcopy(npoints, outfield2, 1, Exp->UpdatePhys(), 1);
                Exp->FwdTrans_IterPerExp(outfield2, Exp->UpdateCoeffs());
                it = m_f->m_exp.begin() + s * (nfields + 2) + nfields;
                m_f->m_exp.insert(it, Exp);

                Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
                Vmath::Vcopy(npoints, outfield3, 1, Exp->UpdatePhys(), 1);
                Exp->FwdTrans_IterPerExp(outfield3, Exp->UpdateCoeffs());
                it = m_f->m_exp.begin() + s * (nfields + 3) + nfields;
                m_f->m_exp.insert(it, Exp);*/
            }
        }

    }
}
