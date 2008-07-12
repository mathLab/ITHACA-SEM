///////////////////////////////////////////////////////////////////////////////
//
// File ForwardMultiStepIntegrationScheme.cpp
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
#include <LibUtilities/Foundations/ForwardMultiStepIntegrationScheme.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {  
        ForwardMultiStepIntegrationScheme::ForwardMultiStepIntegrationScheme(const TimeIntegrationSchemeKey &tkey):
            TimeIntegrationScheme(tkey),
            m_beta(4)
        {
            switch(tkey.GetIntegrationSchemeType())
            {
            case eAdamsBashforth:
                {

                    switch(tkey.GetOrder())
                    {
                    case 1:
                        {
                            m_beta[0] = 0.0;
                            m_beta[1] = 1.0;
                            m_beta[2] = 0.0;
                            m_beta[3] = 0.0;
                        }
                        break;
                    case 2:
                        {
                            m_beta[0] = 0.0;
                            m_beta[1] = 3.0/2.0;
                            m_beta[2] = -1.0/2.0;
                            m_beta[3] = 0.0;
                        }
                        break;
                    case 3:
                        {
                            m_beta[0] = 0.0;
                            m_beta[1] = 23.0/12.0;
                            m_beta[2] = -16.0/12.0;
                            m_beta[3] = 5.0/12.0;
                        }
                        break;
                    default:
                        {
                            NEKERROR(ErrorUtil::efatal,"The Adams-Bashforth scheme is only supported for orders 1,2 and 3");
                        }
                    }

                }
                break;
            case eAdamsMoulton:
                {

                    switch(tkey.GetOrder())
                    {
                    case 1:
                        {
                            m_beta[0] = 1.0;
                            m_beta[1] = 0.0;
                            m_beta[2] = 0.0;
                            m_beta[3] = 0.0;
                        }
                        break;
                    case 2:
                        {
                            m_beta[0] = 1.0/2.0;
                            m_beta[1] = 1.0/2.0;
                            m_beta[2] = 0.0;
                            m_beta[3] = 0.0;
                        }
                        break;
                    case 3:
                        {
                            m_beta[0] = 5.0/12.0;
                            m_beta[1] = 8.0/12.0;
                            m_beta[2] = -1.0/12.0;
                            m_beta[3] = 0.0;
                        }
                        break;
                    default:
                        {
                            NEKERROR(ErrorUtil::efatal,"The Adams-Moulton scheme is only supported for orders 1,2 and 3");
                        }
                    }

                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,"Invalid integration type");
                }
            }
        }

        boost::shared_ptr<TimeIntegrationScheme> ForwardMultiStepIntegrationScheme::Create(const TimeIntegrationSchemeKey &tkey)
        {
            boost::shared_ptr<TimeIntegrationScheme> returnval(new ForwardMultiStepIntegrationScheme(tkey));
            return returnval;
        }
      
    } // end of namespace LibUtilities
} // end of namespace Nektar



