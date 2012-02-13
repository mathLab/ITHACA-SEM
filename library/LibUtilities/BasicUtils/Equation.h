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
//  Description:  Wrapper to ExpressionEvaluator class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIBUTILITIES_EQUATION_HPP
#define NEKTAR_LIBUTILITIES_EQUATION_HPP

#include <string>
#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
namespace Nektar
{
    namespace LibUtilities
    {
        class Equation
        {
        public: 
            Equation(const Equation &src):
              m_expr   (src.m_expr),
              m_expr_id(src.m_expr_id)
            {
            }

            Equation(const std::string& expr = ""):
              m_expr(expr),
              m_expr_id(-1)
            {
                try
                {
                    if (!expr.empty())
                    {
                        m_expr_id = m_evaluator.DefineFunction("x y z t", m_expr);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    m_expr_id = -1;
                    std::cout << "Equation::Constructor fails on expression [" << m_expr << "]" << std::endl;

                    // this instanse is wrongly used by {DisContField2D,BoundaryCondition,...}
                    // classes as a wrapper to a string container holding the link to the file
                    // with boundary conditions or type modifier for the solver-dependent
                    // type of boundary conditions.
                    //
                    // AnalyticExpressionEvaluator cannot parse the expression of the form
                    // e.g. "FILE:whatever.bc" and throws an instance of std::runtime_error.
                    // In case of "TimeDependent" and others it throws parse exception
                    // 'illegal parameter specified' since it looks like an unspecified
                    // parameter in an analytic expression.
                    //
                    // In order not to destroy the currently existing code we catch
                    // this exception in order to ignore it. Hope it does not affect
                    // the performance that much, let the run continue; we'll fix wrong
                    // initialisers later

                    if (expr.find("FILE:") != std::string::npos)
                    {
                        return;
                    }
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return;
                }
                catch (const std::string& e)
                {
                    m_expr_id = -1;
                    std::cout << "Equation::Constructor fails on expression [" << m_expr << "]" << std::endl;

                    std::cout << "ERROR: " << e << std::endl;
                    return;
                }
            }

            Equation& operator=(const Equation &src)
            {
                return *this;
            }

            NekDouble Evaluate0() const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator.Evaluate0(m_expr_id);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return 0;
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    std::cout << "ERROR: " << e << std::endl;
                    return 0;
                }
            }

            NekDouble Evaluate(NekDouble x=0, NekDouble y=0, NekDouble z=0, NekDouble t=0) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator.Evaluate4(m_expr_id, x,y,z,t);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return 0;
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    std::cout << "ERROR: " << e << std::endl;
                    return 0;
                }
            }

            Array<OneD, NekDouble> Evaluate4Array(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const Array<OneD, const NekDouble>& t) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator.Evaluate4Array(m_expr_id, x,y,z,t);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    Array<OneD, NekDouble> empty;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return empty;
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    Array<OneD, NekDouble> empty;
                    std::cout << "ERROR: " << e << std::endl;
                    return empty;
                }
            }

            static void SetConstParameters(const std::map<std::string, NekDouble> &constants)
            {
                m_evaluator.AddConstants(constants);
            }

            std::string GetEquation(void) const
            {
              return m_expr;
            }

        private:
            std::string  m_expr;
            int          m_expr_id;
            LIB_UTILITIES_EXPORT static LibUtilities::AnalyticExpressionEvaluator m_evaluator;
        };

        typedef boost::shared_ptr<Equation> EquationSharedPtr;
    }
}

#endif //NEKTAR_LIBUTILITIES_EQUATION_HPP
