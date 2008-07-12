///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationScheme.h
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
// Description: Header file of time integration scheme base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
#define NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H

#include <math.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        class TimeIntegrationSchemeKey
        {
        public:

            struct opLess
            {
                bool operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;
            };
            
            TimeIntegrationSchemeKey(const TimeIntegrationType &integrationtype, const int &order): 
                m_order(order), 
                m_integrationtype(integrationtype)
            {
            }
            
            virtual ~TimeIntegrationSchemeKey()
            {
            }

            TimeIntegrationSchemeKey(const TimeIntegrationSchemeKey &key)
            {
                *this = key; // defer to assignment operator
            }

            TimeIntegrationSchemeKey& operator=(const TimeIntegrationSchemeKey &key)
            {
                m_order = key.m_order;
                m_integrationtype  = key.m_integrationtype;

                return *this;
            }

            inline unsigned int GetOrder() const
            {
                return m_order;
            }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_integrationtype;
            }

            inline bool operator==(const TimeIntegrationSchemeKey &key)
            {
                return (m_order == key.m_order &&
                    m_integrationtype == key.m_integrationtype);
            }

            inline bool operator== (const TimeIntegrationSchemeKey *y)
            {
                return (*this == *y);
            }

            inline bool operator != (const TimeIntegrationSchemeKey& y)
            {
                return (!(*this == y));
            }

            inline bool operator != (const TimeIntegrationSchemeKey *y)
            {
                return (!(*this == *y));
            }

            friend bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            friend bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
            friend bool opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const;

        protected:
            unsigned int m_order;                   //!< order of the integration scheme
            TimeIntegrationType m_integrationtype;  //!< Type of the integration scheme

        private:
            // This should never be called
            TimeIntegrationSchemeKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationSchemeKey class should not be called");
            }
        };

        static const TimeIntegrationSchemeKey NullTimeIntegrationSchemeKey(eNoTimeIntegrationType, 0);

        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs);

        class TimeIntegrationScheme
        {
        public:

            virtual ~TimeIntegrationScheme()
            {
            }

            inline unsigned int GetOrder() const
            {
                return m_schemeKey.GetOrder();
            }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_schemeKey.GetIntegrationSchemeType();
            }

            virtual const Array<OneD, const NekDouble>& GetBetaCoefficients(void) const
            {
                NEKERROR(ErrorUtil::efatal,"This function is not defined for this type of time integration.");
            }  

            virtual const Array<TwoD, const NekDouble>& GetAcoefficients(void) const
            {
                NEKERROR(ErrorUtil::efatal,"This function is not defined for this type of time integration.");
            }       

            virtual const Array<OneD, const NekDouble>& GetBcoefficients(void) const
            {
                NEKERROR(ErrorUtil::efatal,"This function is not defined for this type of time integration.");
            }       

            virtual const Array<OneD, const NekDouble>& GetCcoefficients(void) const
            {
                NEKERROR(ErrorUtil::efatal,"This function is not defined for this type of time integration.");
            }

        protected:
            TimeIntegrationSchemeKey m_schemeKey;

            TimeIntegrationScheme(const TimeIntegrationSchemeKey &key):m_schemeKey(key)
            {
            }

            // These should never be called
            TimeIntegrationScheme(const TimeIntegrationScheme &in):m_schemeKey(in.m_schemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Copy Constructor for the TimeIntegrationScheme class should not be called");
            }

            TimeIntegrationScheme():m_schemeKey(NullTimeIntegrationSchemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationScheme class should not be called");
            }

        private:
        };

        typedef boost::shared_ptr<TimeIntegrationScheme> TimeIntegrationSchemeSharedPtr;

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
