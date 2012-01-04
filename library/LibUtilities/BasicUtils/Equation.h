////////////////////////////////////////////////////////////////////////////////
//
//  File:  Equation.h
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
#ifndef NEKTAR_LIBUTILITIES_EQUATION_HPP
#define NEKTAR_LIBUTILITIES_EQUATION_HPP

#include <string>
#include <LibUtilities/Interpreter/ExpressionEvaluator.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
namespace Nektar
{
    namespace LibUtilities
    {
        class Equation
        {
        public: 

            Equation()
            {
            }

            Equation(const Equation &eqn):
              m_eqn(eqn.m_eqn)
            {
            }

            Equation(const std::string &eqn):
              m_eqn(eqn)
            {
            }

            NekDouble Evaluate(NekDouble x=0, NekDouble y=0, NekDouble z=0, NekDouble t=0) const
            {
                try
                {
                    m_evaluator.DefineFunction("x y z t", m_eqn.c_str());
                    return m_evaluator.Evaluate(x, y, z, t);
                }
                catch (const std::runtime_error& e)
                {
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return 0;
                }
                catch (const std::string& e)
                {
                    std::cout << "ERROR: " << e << std::endl;
                    return 0;
                }
            }

            static void SetConstParameters(const std::map<std::string, NekDouble> &constants)
            {
                m_evaluator.AddConstants(constants);
            }

            std::string GetEquation(void) const
            {
              return m_eqn;
            }

            void SetEquation(std::string eqn)
            {
              m_eqn = eqn;
            }

        private:
            std::string m_eqn;
            LIB_UTILITIES_EXPORT static LibUtilities::ExpressionEvaluator m_evaluator;
        };

        typedef boost::shared_ptr<Equation> EquationSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_EQUATION_HPP
