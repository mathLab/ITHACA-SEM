////////////////////////////////////////////////////////////////////////////////
//
//  File:  Equation.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description:  Wrapper to Interpreter class
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/Interpreter/Interpreter.h>

#include <boost/algorithm/string/trim.hpp>

namespace Nektar
{
namespace LibUtilities
{

Equation::Equation(InterpreterSharedPtr evaluator,
                   const std::string& expr,
                   const std::string& vlist):
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
        NEKERROR(ErrorUtil::efatal, msg);
        throw e;
        return;
    }
    catch (const std::string& e)
    {
        m_expr_id = -1;
        std::string msg(std::string("Equation::Equation() fails on expression [") + m_expr + std::string("]\n"));
        NEKERROR(ErrorUtil::efatal, msg);
        throw e;
        return;
    }
}

Equation& Equation::operator=(const Equation &src)
{
    m_vlist     = src.m_vlist;
    m_expr      = src.m_expr;
    m_expr_id   = src.m_expr_id;
    m_evaluator = src.m_evaluator;
    return *this;
}

NekDouble Equation::Evaluate() const
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
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e.what());
    }
    catch (const std::string& e)
    {
        std::string msg(std::string("Equation::Evaluate fails on expression [") + m_expr + std::string("]\n"));
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e);
    }
    return 0;
}

NekDouble Equation::Evaluate(
    NekDouble x, NekDouble y, NekDouble z, NekDouble t) const
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
        std::string msg =
            "Equation::Evaluate fails on expression [" + m_expr + "]\n";
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e.what());
    }
    catch (const std::string& e)
    {
        std::string msg =
            "Equation::Evaluate fails on expression [" + m_expr + "]\n";
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e);
    }

    return 0;
}

void Equation::Evaluate(
    const Array<OneD, const NekDouble>& x,
    const Array<OneD, const NekDouble>& y,
    const Array<OneD, const NekDouble>& z,
    Array<OneD, NekDouble>& result) const
{
    Array<OneD, NekDouble>  zero(x.size(), 0.0);
    Evaluate(x,y,z,zero, result);
}

void Equation::Evaluate(
    const Array<OneD, const NekDouble>& x,
    const Array<OneD, const NekDouble>& y,
    const Array<OneD, const NekDouble>& z,
    const NekDouble t,
    Array<OneD, NekDouble>& result) const
{
    Array<OneD, NekDouble>  time(x.size(), t);
    Evaluate(x,y,z,time, result);
}


void Equation::Evaluate(
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

void Equation::Evaluate(
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
        std::string msg =
            "Equation::Evaluate fails on expression [" + m_expr + "]\n";
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e.what());
        return;
    }
    catch (const std::string& e)
    {
        std::string msg =
            "Equation::Evaluate fails on expression [" + m_expr + "]\n";
        NEKERROR(ErrorUtil::efatal, msg + std::string("ERROR: ") + e);
        return;
    }
}

void Equation::SetParameter(const std::string& name, NekDouble value)
{
    m_evaluator->SetParameter(name, value);
}

void Equation::SetConstants(const std::map<std::string, NekDouble> &constants)
{
    m_evaluator->AddConstants(constants);
}

std::string Equation::GetExpression(void) const
{
    return m_expr;
}

std::string Equation::GetVlist(void) const
{
    return m_vlist;
}

/// Returns time spend on expression evaluation at
/// points (it does not include parse/pre-processing time).
NekDouble Equation::GetTime() const
{
    return m_evaluator->GetTime();
}

}
}
