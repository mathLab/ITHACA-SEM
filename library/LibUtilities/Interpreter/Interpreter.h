///////////////////////////////////////////////////////////////////////////////
//
// File: Interpreter.h
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
// Description: Parser and evaluator of analytic expressions.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_INTERPRETER_INTERPRETER_H
#define NEKTAR_LIBUTILITIES_INTERPRETER_INTERPRETER_H

#include <string>
#include <memory>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace LibUtilities
{

class Interpreter
{
public:
    Interpreter();
    ~Interpreter();
    Interpreter(Interpreter &&) = default;
    Interpreter &operator=(Interpreter &&) = default;

    LIB_UTILITIES_EXPORT void SetRandomSeed(unsigned int seed = 123u);

    /// Constants are evaluated and inserted into the function at the time it is
    /// parsed when calling the #DefineFunction function. After parsing, if a
    /// constant is changed, it will not be reflected in the function when
    /// Evaluate is called. This also means that if a function with an unknown
    /// constant is added, and then the constant is added, the function will not
    /// see the added constant and through an exception. This function will add
    /// all of the constants in the map argument to the global internal
    /// constants. If a constant was already loaded previously, it will throw an
    /// exception stating which constants in the map had this issue. It will add
    /// all of the constants it can and output the constants it couldn't add in
    /// the string exception.
    LIB_UTILITIES_EXPORT void AddConstants(std::map<std::string, NekDouble> const& constants);

    /// This function behaves in the same way as #AddConstants, but it only adds
    /// one constant at a time. If the constant existed previously, an exception
    /// will be thrown stating the fact. If it did not exist previously, it will
    /// be added to the global constants and will be used the next time
    /// #DefineFunction is called.
    LIB_UTILITIES_EXPORT int AddConstant(std::string const& name, NekDouble value);

    /// If a constant with the specified name exists, it returns the NekDouble
    /// value that the constant stores. If the constant doesn't exist, it throws
    /// an exception.
    LIB_UTILITIES_EXPORT NekDouble GetConstant(std::string const& name);

    /// Parameters are like constants, but they are inserted into the function
    /// at the time
    /// #Evaluate is called instead of when the function is parsed. This
    /// #function can
    /// be called at any time, and it will take effect in the next call to
    /// #Evaluate.  This function will delete all of the parameters, and replace
    /// all of them with only the ones in the map argument.
    LIB_UTILITIES_EXPORT void SetParameters(std::map<std::string, NekDouble> const& params);

    /// This function behaves in the same way as #SetParameters, but it only
    /// adds one parameter and it does not delete the others. If the parameter
    /// existed previously, it will be overridden and replaced with the new
    /// value. If it did not exist previously, it will be added to the current
    /// parameters.
    LIB_UTILITIES_EXPORT void SetParameter(std::string const& name, NekDouble value);

    /// If a parameter with the specified name exists, it returns the NekDouble
    /// value that the parameter stores. If the parameter doesn't exist, it
    /// throws an exception.
    LIB_UTILITIES_EXPORT NekDouble GetParameter(std::string const& name);

    /// Returns the total time spent in evaluation procedures, seconds.
    LIB_UTILITIES_EXPORT NekDouble GetTime() const;

    // ======================================================
    //  Parsing and evaluation methods
    // ======================================================

    ///  This function allows one to define a function to evaluate. The first argument (vlist)
    ///  is a list of variables (separated by spaces) that the second argument (function)
    ///  depends on. For example, if function = "x + y", then vlist should most likely be
    ///  "x y", unless you are defining x or y as parameters with #SetParameters.
    ///  \output   parsed expression ID. You will need this expression id to call evaluation
    ///            methods below.
    LIB_UTILITIES_EXPORT int DefineFunction(const std::string& vlist, const std::string& function);

    ///  Evaluation method for expressions depending on parameters only.
    LIB_UTILITIES_EXPORT NekDouble Evaluate(const int AnalyticExpression_id);

    ///  Evaluation method for expressions depending on 4 variables (+parameters).
    LIB_UTILITIES_EXPORT NekDouble Evaluate(
        const int AnalyticExpression_id,
        const NekDouble,
        const NekDouble,
        const NekDouble,
        const NekDouble);

    ///  Evaluation method for expressions depending on unspecified number of variables.
    ///  This suitable for expressions depending on more than 4 variables or for the dynamic
    ///  setting some variables as parameters (there is currently no interface method
    ///  for removing a variable from parameter map though).
    LIB_UTILITIES_EXPORT NekDouble EvaluateAtPoint(
        const int AnalyticExpression_id,
        std::vector<NekDouble> point);

    ///  Vectorized evaluation method for expressions depending on 4 variables.
    LIB_UTILITIES_EXPORT void Evaluate(
        const int expression_id,
        const Array<OneD, const NekDouble>&,
        const Array<OneD, const NekDouble>&,
        const Array<OneD, const NekDouble>&,
        const Array<OneD, const NekDouble>&,
        Array<OneD, NekDouble>& result);

    ///  Vectorized evaluation method for expressions depending on unspecified
    ///  number of variables.
    LIB_UTILITIES_EXPORT void Evaluate(
        const int expression_id,
        const std::vector<Array<OneD, const NekDouble> > &points,
        Array<OneD, NekDouble>& result);

private:
    /// Forward declaration of evaluator to avoid boost::spirit includes.
    class ExpressionEvaluator;
    /// Concrete implementation of the above API calls.
    std::unique_ptr<ExpressionEvaluator> m_impl;
};

typedef std::shared_ptr<Interpreter> InterpreterSharedPtr;

}
}

#endif
