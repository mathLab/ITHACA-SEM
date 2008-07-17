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
            
        TimeIntegrationSchemeKey(const TimeIntegrationType &integrationtype): 
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
            m_integrationtype  = key.m_integrationtype;

            return *this;
        }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_integrationtype;
            }

            inline bool operator==(const TimeIntegrationSchemeKey &key)
            {
                return (m_integrationtype == key.m_integrationtype);
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
            TimeIntegrationType m_integrationtype;  //!< Type of the integration scheme

        private:
            // This should never be called
            TimeIntegrationSchemeKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationSchemeKey class should not be called");
            }
        };

        static const TimeIntegrationSchemeKey NullTimeIntegrationSchemeKey(eNoTimeIntegrationType);

        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs);
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs);

        class TimeIntegrationScheme
        {
        public:

            static boost::shared_ptr<TimeIntegrationScheme> Create(const TimeIntegrationSchemeKey &key);

            virtual ~TimeIntegrationScheme()
            {
            }

            inline TimeIntegrationType GetIntegrationSchemeType() const
            {
                return m_schemeKey.GetIntegrationSchemeType();
            }

            inline const Array<TwoD, const NekDouble>& GetA(void) const
            {
                return m_A;
            }  

            inline const Array<TwoD, const NekDouble>& GetB(void) const
            {
                return m_B;
            }       

            inline const Array<TwoD, const NekDouble>& GetU(void) const
            {
                return m_U;
            }       

            inline const Array<TwoD, const NekDouble>& GetV(void) const
            {
                return m_V;
            }

            inline unsigned int GetNsteps(void) const
            {
                return m_numsteps;
            }

            inline unsigned int GetNstages(void) const
            {
                return m_numstages;
            }

        protected:
            TimeIntegrationSchemeKey m_schemeKey;
            unsigned int             m_numsteps;
            unsigned int             m_numstages;

            Array<TwoD,NekDouble>    m_A;
            Array<TwoD,NekDouble>    m_B;
            Array<TwoD,NekDouble>    m_U;
            Array<TwoD,NekDouble>    m_V;

        private: 
            
            TimeIntegrationScheme(const TimeIntegrationSchemeKey &key);
            
            // These should never be called
            TimeIntegrationScheme():m_schemeKey(NullTimeIntegrationSchemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for the TimeIntegrationScheme class should not be called");
            }
            
            TimeIntegrationScheme(const TimeIntegrationScheme &in):m_schemeKey(NullTimeIntegrationSchemeKey)
            {
                NEKERROR(ErrorUtil::efatal,"Copy Constructor for the TimeIntegrationScheme class should not be called");
            }
        };

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs);

        typedef boost::shared_ptr<TimeIntegrationScheme> TimeIntegrationSchemeSharedPtr;
        typedef std::vector< TimeIntegrationSchemeSharedPtr > TimeIntegrationSchemeVector; 
        typedef std::vector< TimeIntegrationSchemeSharedPtr >::iterator TimeIntegrationSchemeVectorIter; 

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_TIMEINTEGRATIONSCHEME_H
