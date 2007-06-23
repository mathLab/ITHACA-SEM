////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/Equation.hpp,v $
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
//  Description:  Template class for holding and evaluating classes.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_EQUATION_HPP
#define NEKTAR_SPATIALDOMAINS_EQUATION_HPP

#include <string>

namespace Nektar
{
    namespace SpatialDomains
    {
        template <typename T>
        class ConstantEvaluator
        {
        public:
            static T Evaluate(const std::string &eqn, T x=0, T y=0, T z=0)
            {
                return atof(eqn.c_str());
            }
        };

        template <typename T, template<typename> class Evaluator = ConstantEvaluator>
        class Equation
        {
        public:
            Equation(const std::string &eqn):
              m_eqn(eqn)
              {
              }

              T Evaluate(T x=0, T y=0, T z=0)
              {
                  return Evaluator<T>::Evaluate(m_eqn, x, y, z);
              }

              std::string GetEquation(void)
              {
                  return m_eqn;
              }

              void SetEquation(std::string eqn)
              {
                  m_eqn = eqn;
              }

              void SetEquation(const char *eqn)
              {
                  m_eqn = eqn;
              }

        private:
            std::string m_eqn;
        };
    }
}

#endif NEKTAR_SPATIALDOMAINS_EQUATION_HPP