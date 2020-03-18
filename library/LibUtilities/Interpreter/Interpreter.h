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

/**
 * @brief Interpreter class for the evaluation of mathematical expressions.
 *
 * The Interpreter class uses the boost::spirit parser framework in order to
 * construct a simple grammar for the evaluation of mathematical
 * expressions. Primarily this is used within LibUtilities::SessionReader to
 * evaluate functions defined using strings, but can be used generically
 * anywhere within the framework in order to provide parsing functionality.
 *
 * The interpreter supports:
 * - classical mathematical binary operators of addition `+`, subtraction `-`,
 *   multiplication `*`, division `/` and exponentiation `^`, evaluated in the
 *   appropriate manner and with use of parentheses;
 * - logical operators `&&`, `||`, `==`, `<=` and `>=`;
 * - several predefined constants such as `PI` and other related constants from
 *   the `cmath` header assigned by the #AddConstant function;
 * - various mathematical functions including trigonometric functions (`sin`,
 *   `cos`, etc), `pow`, `log` and `exp`;
 * - noise generation with the `awgn` function;
 * - the evaluation of spatially-constant parameters, controlled by the
 *   #SetParameters function;
 * - support for either single-value or multiple-value variables, in order to
 *   e.g. support the evaluation of spatially-variable expressions.
 *
 * Note that this class exposes the interface of the Interpreter class only in
 * order to reduce compile-time includes. The implementation itself can be found
 * in Interpreter::ExpressionEvaluator.
 */
class Interpreter
{
public:
    /**
     * @brief Default constructor.
     */
    LIB_UTILITIES_EXPORT Interpreter();

    /**
     * @brief Default destructor.
     */
    LIB_UTILITIES_EXPORT ~Interpreter();

    /**
     * @brief Default move constructor.
     */
    LIB_UTILITIES_EXPORT Interpreter(Interpreter &&);

    /**
     * @brief Default assignment operator.
     */
    LIB_UTILITIES_EXPORT Interpreter &operator=(Interpreter &&);

    /**
     * @brief Sets the random seed for the pseudorandom number generator.
     *
     * This allows for e.g. different ranks to be given different seeds to
     * ensure appropriate entropy in noise generation.
     */
    LIB_UTILITIES_EXPORT void SetRandomSeed(unsigned int seed = 123u);

    /**
     * @brief Set constants to be evaluated.
     *
     * Constants are evaluated and inserted into the function at the time it is
     * parsed when calling the #DefineFunction function. After parsing, if a
     * constant is changed, it will not be reflected in the function when
     * Evaluate is called. This also means that if a function with an unknown
     * constant is added, and then the constant is added, the function will not
     * see the added constant and through an exception.
     *
     * This function will add all of the constants in the @p constants parameter
     * to the global internal constants. If a constant was already loaded
     * previously, it will throw an exception stating which constants in the map
     * had this issue. It will add all of the constants it can from @p constants
     * and output the constants it couldn't add in the string exception.
     *
     * @param constants  A std::map with string names for the constants (which
     *                   will be evalauted in the expression) and their
     *                   NekDouble value.
     */
    LIB_UTILITIES_EXPORT void AddConstants(
        std::map<std::string, NekDouble> const& constants);

    /**
     * @brief Set constants to be evaluated.
     *
     * This function behaves in the same way as #AddConstants, but it only adds
     * one constant at a time. If the constant existed previously, an exception
     * will be thrown stating the fact. If it did not exist previously, it will
     * be added to the global constants and will be used the next time
     * #DefineFunction is called.
     *
     * @return Total number of constants after this one has been added.
     */
    LIB_UTILITIES_EXPORT int AddConstant(
        std::string const& name, NekDouble value);

    /**
     * @brief Return the value of a constant.
     *
     * If a constant with the specified name @p name exists, this function
     * returns the NekDouble value that the constant stores. If the constant
     * doesn't exist, this throws an exception.
     *
     * @param name  Name of constant to return.
     */
    LIB_UTILITIES_EXPORT NekDouble GetConstant(std::string const& name);

    /**
     * @brief Set parameter values.
     *
     * Parameters are functionally similar to constants, but they are inserted
     * into the function at the time that #Evaluate is called, instead of when
     * the function is parsed. This function can therefore be called at any
     * time, and it will take effect in the next call to #Evaluate.  This
     * function will delete all of the parameters, and replace all of them with
     * only the ones in the map argument.
     */
    LIB_UTILITIES_EXPORT void SetParameters(
        std::map<std::string, NekDouble> const& params);

    /**
     * @brief Set parameter values.
     *
     * This function behaves in the same way as #SetParameters, but it only adds
     * one parameter and it does not delete the others. If the parameter @p name
     * existed previously, it will be overridden and replaced with the new
     * value. If it did not exist previously, it will be added to the current
     * parameters.
     *
     * @param name   Name of the parameter to define.
     * @param value  The parameter's value.
     */
    LIB_UTILITIES_EXPORT void SetParameter(
        std::string const& name, NekDouble value);

    /**
     * @brief Get the value of a parameter
     *
     * If a parameter with the specified @p name exists, it returns the
     * NekDouble value that the parameter stores. If the parameter doesn't
     * exist, it throws an exception.
     *
     * @param name  Name of the parameter to query.
     */
    LIB_UTILITIES_EXPORT NekDouble GetParameter(std::string const& name);

    /**
     * @brief Returns the total walltime spent in evaluation procedures in
     * seconds.
     */
    LIB_UTILITIES_EXPORT NekDouble GetTime() const;

    // ======================================================
    //  Parsing and evaluation methods
    // ======================================================

    /**
     * @brief Defines a function for the purposes of evaluation.
     *
     * This function allows one to define a function to evaluate. The @p vlist
     * argument should define a list of space-separated variables that the
     * expression defined in @p function is dependent upon. For example, if @p
     * function is defined as the string `x + y`, then vlist should most likely
     * be `x y`, unless you are defining `x` or `y` as parameters with the
     * #SetParameters function.
     *
     * @param vlist     List of variable names separated with spaces.
     * @param expr      String definition of the function to be evaluated.
     *
     * @return  An integer denoting the unique ID of this function. This should
     *          be passed into #Evaluate functions for later evaluation.
     */
    LIB_UTILITIES_EXPORT int DefineFunction(
        const std::string& vlist, const std::string& expr);

    /**
     * @brief Evaluate a function which depends only on constants and/or
     * parameters.
     *
     * @param id  The ID returned from #DefineFunction representing the function
     *            to be evaluated.
     */
    LIB_UTILITIES_EXPORT NekDouble Evaluate(const int id);

    /**
     * @brief Evaluate a function which depends on four variables: typically
     * space \f$ (x, y, z) \f$ and time \f$ t\f$.
     *
     * @param id  The ID returned from #DefineFunction representing the function
     *            to be evaluated.
     * @param x   The value of variable 1 (typically \f$ x \f$).
     * @param y   The value of variable 2 (typically \f$ y \f$).
     * @param z   The value of variable 3 (typically \f$ z \f$).
     * @param t   The value of variable 4 (typically \f$ t \f$).
     */
    LIB_UTILITIES_EXPORT NekDouble Evaluate(
        const int id, const NekDouble x, const NekDouble y, const NekDouble z,
        const NekDouble t);

    /**
     * @brief Evaluate a function which depends on zero or more variables.
     *
     * This is suitable for expressions depending on more than 4 variables or
     * for the dynamic setting of some variables as parameters (there is
     * currently no interface method for removing a variable from parameter map
     * however).
     *
     * @param id     The ID returned from #DefineFunction representing the
     *               function to be evaluated.
     * @param point  A std::vector of points to be evaluated.
     */
    LIB_UTILITIES_EXPORT NekDouble EvaluateAtPoint(
        const int id,
        std::vector<NekDouble> point);

    /**
     * @brief Evaluate a function which depends on four variables: typically
     * space \f$ (x, y, z) \f$ and time \f$ t\f$ in a vectorised manner.
     *
     * This is a vectorised version of the evaluation method that will allow the
     * same function to be evaluated many times for each of the entries in the
     * input variable arrays @p x, @p y, @p z and @p t. Note that this is
     * typically far quicker than the use of a loop with many calls to the
     * single-variable variants in this class.
     *
     * @note All parameters @p x, @p y, @p z and @p t should be Arrays of the
     *       same length.
     *
     * @param id     The ID returned from #DefineFunction representing the
     *               function to be evaluated.
     * @param x      An Array of values for variable 1 (typically \f$ x\f$).
     * @param y      An Array of values for variable 2 (typically \f$ y\f$).
     * @param z      An Array of values for variable 3 (typically \f$ z\f$).
     * @param t      An Array of values for variable 4 (typically \f$ t\f$).
     * @param result An Array containing the evaluation of the function for each
     *               of the variable values.
     */
    LIB_UTILITIES_EXPORT void Evaluate(
        const int id,
        const Array<OneD, const NekDouble> &x,
        const Array<OneD, const NekDouble> &y,
        const Array<OneD, const NekDouble> &z,
        const Array<OneD, const NekDouble> &t,
        Array<OneD, NekDouble>& result);

    /**
     * @brief Evaluate a function which depends on zero or more variables.
     *
     * This is a vectorised version of the evaluation method that will allow the
     * same function to be evaluated many times for each of the entries in the
     * input variable arrays within @p points. Note that this is typically far
     * quicker than the use of a loop with many calls to the single-variable
     * variants in this class.
     *
     * @note All entries of the input array @p points should be of the same
     *       length.
     *
     * @param id     The ID returned from #DefineFunction representing the
     *               function to be evaluated.
     * @param points An Array in which the i-th entry corresponds to values of
     *               variable \f$ i\f$.
     * @param result An Array containing the evaluation of the function for each
     *               of the variable values.
     */
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
