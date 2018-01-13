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
#include <map>
#include <stdexcept>
#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <boost/algorithm/string/trim.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        class Equation
        {

        public: 
            LIB_UTILITIES_EXPORT Equation(const Equation &src):
              m_vlist     (src.m_vlist),
              m_expr      (src.m_expr),
              m_expr_id   (src.m_expr_id),
              m_evaluator (src.m_evaluator)
            {
            }

            LIB_UTILITIES_EXPORT Equation(ExpressionEvaluatorShPtr evaluator,
                                          const std::string& expr = "",
                                          const std::string& vlist = ""):
              m_vlist     (vlist),
              m_expr      (expr),
              m_expr_id   (-1),
              m_evaluator (evaluator)
            {
                boost::algorithm::trim(m_expr);
                boost::algorithm::trim(m_vlist);

                if (m_vlist.empty())
                {
                    m_vlist = "x y z t";
                }

                try
                {
                    if (!m_expr.empty())
                    {
                        m_expr_id = m_evaluator->DefineFunction(m_vlist, m_expr);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    m_expr_id = -1;
                    std::string msg(std::string("Equation::Equation() fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg);
                    throw e;
                    return;
                }
                catch (const std::string& e)
                {
                    m_expr_id = -1;
                    std::string msg(std::string("Equation::Equation() fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg);
                    throw e;
                    return;
                }
            }

            LIB_UTILITIES_EXPORT Equation& operator=(const Equation &src)
            {
                return *this;
            }

            LIB_UTILITIES_EXPORT NekDouble Evaluate() const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator->Evaluate(m_expr_id);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e.what());
                }
                catch (const std::string& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e);
                }
                return 0;
            }

            LIB_UTILITIES_EXPORT NekDouble Evaluate(NekDouble x, NekDouble y=0, NekDouble z=0, NekDouble t=0) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        return m_evaluator->Evaluate(m_expr_id, x,y,z,t);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e.what());
                }
                catch (const std::string& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e);
                }
                return 0;
            }

            LIB_UTILITIES_EXPORT void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    Array<OneD, NekDouble>& result)
            {
                Array<OneD, NekDouble>  zero(x.num_elements(), 0.0);
                Evaluate(x,y,z,zero, result);
            }

            LIB_UTILITIES_EXPORT void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const NekDouble t,
                    Array<OneD, NekDouble>& result) const
            {
                Array<OneD, NekDouble>  time(x.num_elements(), t);
                Evaluate(x,y,z,time, result);
            }


            LIB_UTILITIES_EXPORT void Evaluate(
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const Array<OneD, const NekDouble>& t,
                    Array<OneD, NekDouble>& result) const
            {
                std::vector<Array<OneD, const NekDouble> > points;
                points.push_back(x);
                points.push_back(y);
                points.push_back(z);
                points.push_back(t);
                Evaluate(points, result);
            }

            LIB_UTILITIES_EXPORT void Evaluate(
                    const std::vector<Array<OneD, const NekDouble> > points,
                    Array<OneD, NekDouble>& result) const
            {
                try
                {
                    if (m_expr_id != -1)
                    {
                        m_evaluator->Evaluate(m_expr_id, points, result);
                    }
                }
                catch (const std::runtime_error& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e.what());
                    return;
                }
                catch (const std::string& e)
                {
                    std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
                    ASSERTL0(false, msg + std::string("ERROR: ") + e);
                    return;
                }
            }

            LIB_UTILITIES_EXPORT void SetParameter(const std::string& name, NekDouble value)
            {
                m_evaluator->SetParameter(name, value);
            }

            LIB_UTILITIES_EXPORT void SetConstants(const std::map<std::string, NekDouble> &constants)
            {
                m_evaluator->AddConstants(constants);
            }

            LIB_UTILITIES_EXPORT std::string GetExpression(void) const
            {
                return m_expr;
            }

            LIB_UTILITIES_EXPORT std::string GetVlist(void) const
            {
                return m_vlist;
            }

            /// Returns time spend on expression evaluation at
            /// points (it does not include parse/pre-processing time).
            LIB_UTILITIES_EXPORT NekDouble GetTime() const
            {
                return m_evaluator->GetTime();
            }

        private:
            std::string  m_vlist;
            std::string  m_expr;
            int          m_expr_id;
            ExpressionEvaluatorShPtr m_evaluator;
        };
    }
}

#endif //NEKTAR_LIBUTILITIES_EQUATION_HPP
