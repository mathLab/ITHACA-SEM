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
//  Description:  Wrapper to AnalyticExpressionEvaluator class.
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
//#include <loki/Singleton.h>

namespace Nektar
{
    namespace LibUtilities
    {
        class Equation
        {
  //      typedef Loki::SingletonHolder<LibUtilities::AnalyticExpressionEvaluator, Loki::CreateStatic, Loki::DefaultLifetime> SingleExpressionEvaluator;


        public: 
            Equation(const Equation &src):
              m_expr     (src.m_expr),
              m_expr_id  (src.m_expr_id)//,
    //          m_evaluator(src.m_evaluator)
            {
            }

            Equation(const std::string& expr = ""):
              m_expr(expr),
              m_expr_id(-1)//,
//              m_evaluator(SingleExpressionEvaluator::Instance())
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

                    // this instanse is wrongly used by DisContField2D
                    // classes as a wrapper to a string container holding the link to the file
                    // with boundary conditions
                    //
                    // AnalyticExpressionEvaluator cannot parse the expression of the form
                    // e.g. "FILE:whatever.bc" and throws an instance of std::runtime_error.
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

            NekDouble Evaluate() const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator.Evaluate(m_expr_id);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    std::cout << "ERROR: " << e << std::endl;
                }
                return 0;
            }

            NekDouble Evaluate(NekDouble x, NekDouble y=0, NekDouble z=0, NekDouble t=0) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator.Evaluate(m_expr_id, x,y,z,t);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    std::cout << "ERROR: " << e << std::endl;
                }
                return 0;
            }

//            void Evaluate(
//                    const Array<OneD, const NekDouble>& x,
//                    const Array<OneD, const NekDouble>& y,
//                    Array<OneD, NekDouble>& result)
//            {
//                Array<OneD, NekDouble>  zero(x.num_elements(), 0.0);
//                Evaluate(x,y,zero,zero, result);
//            }

            void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    Array<OneD, NekDouble>& result)
            {
                Array<OneD, NekDouble>  zero(x.num_elements(), 0.0);
                Evaluate(x,y,z,zero, result);
            }

            void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const NekDouble t,
                    Array<OneD, NekDouble>& result) const
            {
                Array<OneD, NekDouble>  time(x.num_elements(), t);
                Evaluate(x,y,z,time, result);
            }


            void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const Array<OneD, const NekDouble>& t,
                    Array<OneD, NekDouble>& result) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        m_evaluator.Evaluate(m_expr_id, x,y,z,t, result);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    ASSERTL0(false, std::string("ERROR: ") + e.what());
                    return;
                }
                catch (const std::string& e)
                {
                    std::cout << "Equation::Evaluate fails on expression [" << m_expr << "]" << std::endl;
                    std::cout << "ERROR: " << e << std::endl;
                    return;
                }
            }

            void SetParameter(const std::string& name, double value)
            {
                m_evaluator.SetParameter(name, value);
            }

            void SetConstants(const std::map<std::string, NekDouble> &constants)
            {
                m_evaluator.AddConstants(constants);
            }

            std::string GetExpression(void) const
            {
              return m_expr;
            }

            double GetTime() const
            {
                return m_evaluator.GetTime();
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
