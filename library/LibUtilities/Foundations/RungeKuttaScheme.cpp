///////////////////////////////////////////////////////////////////////////////
//
// File RungeKuttaScheme.cpp
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
// Description: Implementation of forward multi-step time integration schemes
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/LibUtilities.h>
#include <LibUtilities/Foundations/RungeKuttaScheme.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {  
        RungeKuttaScheme::RungeKuttaScheme(const TimeIntegrationSchemeKey &tkey):
            TimeIntegrationScheme(tkey),
            m_a(tkey.GetOrder(),tkey.GetOrder(),0.0),
            m_b(tkey.GetOrder()),
            m_c(tkey.GetOrder())
        {
            ASSERTL0(tkey.GetIntegrationSchemeType() == eRungeKutta,"Invalid integration type");

            switch(tkey.GetOrder())
            {
            case 1:
                {
                    m_a[0][0] = 0.0;
                    m_b[0] = 1.0;
                    m_c[0] = 0.0;
                }
                break;
            case 2:
                {
                    m_a[1][0] = 1.0;

                    m_b[0] = 0.5;
                    m_b[1] = 0.5;

                    m_c[0] = 0.0;
                    m_c[1] = 1.0;
                }
                break;
            case 3:
                {
                    m_a[1][0] = 0.5;
                    m_a[2][0] = -1.0;
                    m_a[2][1] = 2.0;

                    m_b[0] = 1.0/6.0;
                    m_b[1] = 2.0/3.0;
                    m_b[2] = 1.0/6.0;

                    m_c[0] = 0.0;
                    m_c[1] = 0.5;
                    m_c[2] = 1.0;
                }
                break;
            case 4:
                {
                    m_a[1][0] = 0.5;
                    m_a[2][0] = 0.0;
                    m_a[3][0] = 0.0;
                    m_a[2][1] = 0.5;
                    m_a[3][1] = 0.0;
                    m_a[3][2] = 1.0;

                    m_b[0] = 1.0/6.0;
                    m_b[1] = 1.0/3.0;
                    m_b[2] = 1.0/3.0;
                    m_b[3] = 1.0/6.0;

                    m_c[0] = 0.0;
                    m_c[1] = 0.5;
                    m_c[2] = 0.5;
                    m_c[3] = 1.0;
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,"The Adams-Bashforth scheme is only supported for orders 1,2 and 3");
                }
            }
        }

        boost::shared_ptr<TimeIntegrationScheme> RungeKuttaScheme::Create(const TimeIntegrationSchemeKey &tkey)
        {
            boost::shared_ptr<TimeIntegrationScheme> returnval(new RungeKuttaScheme(tkey));
            return returnval;
        }
      
    } // end of namespace LibUtilities
} // end of namespace Nektar



