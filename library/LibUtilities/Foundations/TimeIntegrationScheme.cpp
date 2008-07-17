///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class 
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/LibUtilities.h>
#include <iostream>
#include <LibUtilities/Foundations/TimeIntegrationScheme.h>
#include <LibUtilities/Foundations/Foundations.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {
            return (lhs.m_integrationtype == rhs.m_integrationtype);
        }

        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {            
            return (lhs.m_integrationtype < rhs.m_integrationtype);
        }
        
        bool TimeIntegrationSchemeKey::opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const
        {
            return (lhs.m_integrationtype < rhs.m_integrationtype);
        }

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs)
        {
            os << "Time Integration Scheme: " << TimeIntegrationTypeMap[rhs.GetIntegrationSchemeType()] << endl;

            return os;
        }

        TimeIntegrationSchemeSharedPtr TimeIntegrationScheme::Create(const TimeIntegrationSchemeKey &key)
        {
            TimeIntegrationSchemeSharedPtr returnval(new TimeIntegrationScheme(key));
            return returnval;
        }

        TimeIntegrationScheme::TimeIntegrationScheme(const TimeIntegrationSchemeKey &key):
            m_schemeKey(key)
        {
            switch(key.GetIntegrationSchemeType())
            {
            case eForwardEuler:
            case eAdamsBashforthOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eAdamsBashforthOrder2:
                {
                    m_numsteps = 2;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,0.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,0.0);

                    m_B[0][0] = 3.0/2.0;
                    m_B[1][0] = 1.0;

                    m_U[0][0] = 1.0;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = -0.5;
                }
                break;
            case eBackwardEuler:
            case eAdamsMoultonOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,1.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eAdamsMoultonOrder2:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.5);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eClassiscalRungeKutta4:
                {
                    m_numsteps = 1;
                    m_numstages = 4;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,0.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);

                    m_A[1][0] = 0.5;
                    m_A[2][1] = 0.5;
                    m_A[3][2] = 1.0;

                    m_B[0][0] = 1.0/6.0;
                    m_B[0][1] = 1.0/3.0;
                    m_B[0][2] = 1.0/3.0;
                    m_B[0][3] = 1.0/6.0;
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,"Invalid Time Integration Scheme Type");
                }
            }
        }

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs)
        {
            int i,j;
            int r = rhs.GetNsteps();
            int s = rhs.GetNstages();

            int oswidth = 8;

            os << "Time Integration Scheme: " << TimeIntegrationTypeMap[rhs.GetIntegrationSchemeType()] << endl;
            os << "- number of steps:  " << r << endl;
            os << "- number of stages: " << s << endl;
            os << "General linear method tableau: " << endl;

            for(i = 0; i < s; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetA())[i][j] << " ";
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetU())[i][j] << " ";
                }
                os << endl;
            }
            for(int i = 0; i < (r+s)*oswidth+2; i++)
            {
                os << "-";
            }
            os << endl;
            for(i = 0; i < r; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetB())[i][j] << " ";
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetV())[i][j] << " ";
                }
                os << endl;
            }
            return os;
        }

    }
}

