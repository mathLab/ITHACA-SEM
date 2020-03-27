///////////////////////////////////////////////////////////////////////////////
//
// File: Interpreter.cpp
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

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Interpreter/Interpreter.h>

#define BOOST_SPIRIT_THREADSAFE
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_ast.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>
#include <boost/spirit/include/classic_symbols.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/algorithm/string/trim.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace bsp = boost::spirit::classic;

#include <map>
#include <string>
#include <vector>
#if defined(__INTEL_COMPILER)
#include <mathimf.h>
#else
#include <cmath>
#endif

namespace Nektar
{
namespace LibUtilities
{

// signum function
NekDouble sign(NekDouble arg)
{
    return (arg > 0.0) - (arg < 0.0);
}

// Additive white Gaussian noise function.  Arg: sigma of the zero-mean gaussian
// distribution Attention: this function is not actually used for evaluation
// purposes.
NekDouble awgn(NekDouble sigma)
{
    boost::mt19937 rng;
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
        _normal(rng, boost::normal_distribution<>(0, sigma));
    return _normal();
}

static NekDouble rad(NekDouble x, NekDouble y)
{
    return (x != 0.0 || y != 0.0) ? sqrt(x * x + y * y) : 0.0;
}

static NekDouble ang(NekDouble x, NekDouble y)
{
    return (x != 0.0 || y != 0.0) ? atan2(y, x) : 0.0;
}

// =========================================================================
//             AnalyticExpression definitions for Spirit Parser
// =========================================================================

typedef NekDouble (*PFD)();
typedef NekDouble (*PFD1)(NekDouble);
typedef NekDouble (*PFD2)(NekDouble, NekDouble);
typedef NekDouble (*PFD3)(NekDouble, NekDouble, NekDouble);
typedef NekDouble (*PFD4)(NekDouble, NekDouble, NekDouble, NekDouble);
struct func
{
    func(PFD1 p) : func1(p), size(1){};
    func(PFD2 p) : func2(p), size(2){};
    func(PFD3 p) : func3(p), size(3){};
    func(PFD4 p) : func4(p), size(4){};

    union // Pointer to a function
    {
        PFD1 func1;
        PFD2 func2;
        PFD3 func3;
        PFD4 func4;
    };
    size_t size;
};

/** This struct creates a parser that matches the function definitions from
    math.h. All of the functions accept one of more NekDoubles as arguments and
    returns a NekDouble. **/
static struct functions : bsp::symbols<func>
{
    functions()
    {
        // Add all of the functions from math.h
        add("abs", std::abs) // absolute value
            ("asin", asin)   // arcsin
            ("acos", acos)   // arccos
            ("atan", atan)   // arctan
            ("atan2", atan2) // arctan2
            ("ang", ang)     // angle calculation
            ("ceil", ceil)   // ceiling
            ("cos", cos)     // cosine
            ("cosh", cosh)   // hyperbolic cosine
            ("exp", exp)     // exponential
            ("fabs", fabs)   // absolute value
            ("floor", floor) // floor
            ("log", log)     // natural log
            ("log10", log10) // log base 10
            ("rad", rad)     // radians
            ("sin", sin)     // sine
            ("sinh", sinh)   // hyperbolic sine
            ("sqrt", sqrt)   // square root
            ("tan", tan)     // tangent
            ("tanh", tanh)   // hyperbolic tangent
            // and few more custom functions
            ("sign", sign) // sign
            ("awgn", awgn) // white noise
            ;
    }
} functions_p;

/**
 * @brief Concrete implementation of the interface defined in Interpreter.
 */
class Interpreter::ExpressionEvaluator
{
private:
    struct EvaluationStep;

public:
    typedef std::map<std::string, int> VariableMap;
    typedef std::map<std::string, int> ConstantMap;
    typedef std::map<std::string, int> ParameterMap;
    typedef std::map<std::string, int> ExpressionMap;
    typedef std::map<std::string, int> FunctionNameMap;
    typedef std::vector<EvaluationStep *> ExecutionStack;
    typedef std::pair<bool, NekDouble> PrecomputedValue;
    typedef NekDouble (*OneArgFunc)(NekDouble);
    typedef NekDouble (*TwoArgFunc)(NekDouble, NekDouble);

    typedef bsp::tree_parse_info<std::string::const_iterator,
                                 bsp::node_val_data_factory<NekDouble>>
        ParsedTreeInfo;
    typedef bsp::tree_match<std::string::const_iterator,
                            bsp::node_val_data_factory<NekDouble>>::
        tree_iterator ParsedTreeIterator;

    typedef std::vector<const Array<OneD, const NekDouble> *> VariableArray;

    /**
     * @brief Initializes the evaluator.
     *
     * This routine will initialize the evaluator with some basic default
     * constants,
     */
    ExpressionEvaluator() : m_timer(), m_total_eval_time(0)
    {
        m_state_size = 1;

        // Constant definitions.
        AddConstant("MEANINGLESS", 0.0);
        AddConstant("E", M_E);                        // Natural logarithm
        AddConstant("LOG2E", M_LOG2E);                // log_2 e
        AddConstant("LOG10E", M_LOG10E);              // log_10 e
        AddConstant("LN2", M_LN2);                    // log_e 2
        AddConstant("LN10", M_LN10);                  // log_e 10
        AddConstant("PI", M_PI);                      // pi
        AddConstant("PI_2", M_PI_2);                  // pi/2
        AddConstant("PI_4", M_PI_4);                  // pi/4
        AddConstant("1_PI", M_1_PI);                  // 1/pi
        AddConstant("2_PI", M_2_PI);                  // 2/pi
        AddConstant("2_SQRTPI", M_2_SQRTPI);          // 2/sqrt(pi)
        AddConstant("SQRT2", M_SQRT2);                // sqrt(2)
        AddConstant("SQRT1_2", M_SQRT1_2);            // 1/sqrt(2)
        AddConstant("GAMMA", 0.57721566490153286060); // Euler
        AddConstant("DEG", 57.2957795130823208768);   // deg/radian
        AddConstant("PHI", 1.61803398874989484820);   // golden ratio

        // Function definitions.
        m_functionMapNameToInstanceType["abs"]    = E_ABS;
        m_functionMapNameToInstanceType["asin"]   = E_ASIN;
        m_functionMapNameToInstanceType["acos"]   = E_ACOS;
        m_functionMapNameToInstanceType["atan"]   = E_ATAN;
        m_functionMapNameToInstanceType["atan2"]  = E_ATAN2;
        m_functionMapNameToInstanceType["ang"]    = E_ANG;
        m_functionMapNameToInstanceType["bessel"] = E_BESSEL;
        m_functionMapNameToInstanceType["ceil"]   = E_CEIL;
        m_functionMapNameToInstanceType["cos"]    = E_COS;
        m_functionMapNameToInstanceType["cosh"]   = E_COSH;
        m_functionMapNameToInstanceType["exp"]    = E_EXP;
        m_functionMapNameToInstanceType["fabs"]   = E_FABS;
        m_functionMapNameToInstanceType["floor"]  = E_FLOOR;
        m_functionMapNameToInstanceType["log"]    = E_LOG;
        m_functionMapNameToInstanceType["log10"]  = E_LOG10;
        m_functionMapNameToInstanceType["rad"]    = E_RAD;
        m_functionMapNameToInstanceType["sin"]    = E_SIN;
        m_functionMapNameToInstanceType["sinh"]   = E_SINH;
        m_functionMapNameToInstanceType["sqrt"]   = E_SQRT;
        m_functionMapNameToInstanceType["tan"]    = E_TAN;
        m_functionMapNameToInstanceType["tanh"]   = E_TANH;
        m_functionMapNameToInstanceType["sign"]   = E_SIGN;
        m_functionMapNameToInstanceType["awgn"]   = E_AWGN;

        m_function[E_ABS]     = std::abs;
        m_function[E_ASIN]    = asin;
        m_function[E_ACOS]    = acos;
        m_function[E_ATAN]    = atan;
        m_function[E_CEIL]    = ceil;
        m_function[E_COS]     = cos;
        m_function[E_COSH]    = cosh;
        m_function[E_EXP]     = exp;
        m_function[E_FABS]    = fabs;
        m_function[E_FLOOR]   = floor;
        m_function[E_LOG]     = log;
        m_function[E_LOG10]   = log10;
        m_function[E_SIN]     = sin;
        m_function[E_SINH]    = sinh;
        m_function[E_SQRT]    = sqrt;
        m_function[E_TAN]     = tan;
        m_function[E_TANH]    = tanh;
        m_function[E_SIGN]    = sign;
        m_function2[E_ATAN2]  = atan2;
        m_function2[E_ANG]    = ang;
        m_function2[E_RAD]    = rad;
        m_function2[E_BESSEL] = boost::math::cyl_bessel_j;

        // Note that there is no entry in m_function that corresponds to the
        // awgn function. This is intentional as this function need not be
        // pre-evaluated once!
    }

    /**
     * @brief Destructor that removes all entries from the execution stack.
     */
    ~ExpressionEvaluator(void)
    {
        for (auto &it_es : m_executionStack)
        {
            for (auto &it : it_es)
            {
                delete it;
            }
            it_es.clear();
        }
        m_executionStack.clear();
    }

    /**
     * @copydoc Interpreter::SetRandomSeed
     */
    void SetRandomSeed(unsigned int seed)
    {
        m_generator.seed(seed);
    }

    /**
     * @copydoc Interpreter::AddConstants
     */
    void AddConstants(std::map<std::string, NekDouble> const &constants)
    {
        for (auto const &it : constants)
        {
            AddConstant(it.first, it.second);
        }
    }

    /**
     * @copydoc Interpreter::AddConstant
     */
    int AddConstant(std::string const &name, NekDouble value)
    {
        ConstantMap::const_iterator it = m_constantMapNameToId.find(name);
        if (it == m_constantMapNameToId.end())
        {
            // We are trying to avoid duplicating entries in m_constantParser
            // and m_constants.
            m_constantsParser.add(name.c_str(), value);
            int index                   = m_constant.size();
            m_constantMapNameToId[name] = index;
            m_constant.push_back(value);
            return index;
        }
        else
        {
            if (m_constant[it->second] != value)
            {
                std::string errormsg =
                    "Attempt to add numerically different constants under the "
                    "same name: " +
                    name;
                std::cout << errormsg << std::endl;
            }
        }
        return it->second;
    }

    /**
     * @copydoc Interpreter::GetConstant
     */
    NekDouble GetConstant(std::string const &name)
    {
        NekDouble *value = find(m_constantsParser, name.c_str());
        ASSERTL1(value != NULL, "Constant variable not found: " + name);
        return *value;
    }

    /**
     * @copydoc Interpreter::SetParameters
     */
    void SetParameters(std::map<std::string, NekDouble> const &params)
    {
        for (auto const &it : params)
        {
            SetParameter(it.first, it.second);
        }
    }

    /**
     * @copydoc Interpreter::SetParameter
     */
    void SetParameter(std::string const &name, NekDouble value)
    {
        ParameterMap::const_iterator it = m_parameterMapNameToId.find(name);
        if (it == m_parameterMapNameToId.end())
        {
            m_parameterMapNameToId[name] = m_parameter.size();
            m_parameter.push_back(value);
        }
        else
        {
            // If parameter is known, change its value.
            m_parameter[it->second] = value;
        }
    }

    /**
     * @copydoc Interpreter::GetParameter
     */
    NekDouble GetParameter(std::string const &name)
    {
        ParameterMap::const_iterator it = m_parameterMapNameToId.find(name);
        ASSERTL1(it != m_parameterMapNameToId.end(),
                 "Parameter not found: " + name);
        return m_parameter[it->second];
    }

    /**
     * @copydoc Interpreter::GetTime
     */
    NekDouble GetTime() const
    {
        return m_total_eval_time;
    }

    /**
     * @copydoc Interpreter::DefineFunction
     */
    int DefineFunction(const std::string &vlist, const std::string &expr)
    {
        // If we have previously parsed an identical expression, return its ID.
        auto it = m_parsedMapExprToExecStackId.find(expr);
        if (it != m_parsedMapExprToExecStackId.end())
        {
            // If this function is already defined, don't do anything but return
            // its ID.
            return it->second;
        }

        // Otherwise, prepare an iterator that allows to walk along the string
        // representing an analytic expression in the order that respects its
        // recursive structure (thanks to boost::spirit).

        // Parse the vlist input and separate the variables into ordered entries
        // in a vector<string> object. These need to be ordered because this is
        // the order the variables will get assigned to in the Map when
        // Evaluate(...)  is called.
        std::vector<std::string> variableNames;
        bsp::parse((char *)vlist.c_str(),
                   (*bsp::space_p >>
                    *(+(+bsp::graph_p)[bsp::push_back_a(variableNames)] >>
                      +bsp::space_p)));

        // Set up our grammar.
        AnalyticExpression myGrammar(&m_constantsParser, variableNames);

        // Do the actual parsing with boost::spirit and alert the user if there
        // was an error with an exception.
        ParsedTreeInfo parseInfo =
            bsp::ast_parse<bsp::node_val_data_factory<NekDouble>,
                           std::string::const_iterator, AnalyticExpression,
                           bsp::space_parser>(expr.begin(), expr.end(),
                                              myGrammar, bsp::space_p);

        ASSERTL1(parseInfo.full != false,
                 "Unable to fully parse function. Stopped just before: " +
                     std::string(parseInfo.stop, parseInfo.stop + 15));

        if (!parseInfo.full)
        {
            throw std::runtime_error(
                "Unable to fully parse function at: " +
                std::string(parseInfo.stop, parseInfo.stop + 15));
        }

        // ----------------------------------------------
        // Data parsed, start setting up internal data structures.
        // ----------------------------------------------

        ExecutionStack stack;
        VariableMap variableMap;

        int stackId  = m_executionStack.size();
        m_state_size = 1;

        // Register all variables declared in the expression
        for (int i = 0; i < variableNames.size(); i++)
        {
            variableMap[variableNames[i]] = i;
        }

        // Then prepare an execution stack. This method also calculates a
        // length of internal state storage (m_state_size) for this function.
        PrecomputedValue v = PrepareExecutionAsYouParse(parseInfo.trees.begin(),
                                                        stack, variableMap, 0);

        // Constant expression, fully evaluated already.
        if (true == v.first)
        {
            ASSERTL1(stack.size() == 0,
                     "Constant expression yeilds non-empty execution stack. "
                     "Bug in PrepareExecutionAsYouParse()");

            int const_index = AddConstant(
                std::string("EXPRESSION_") + std::to_string(stackId), v.second);
            stack.push_back(makeStep<StoreConst>(0, const_index));
        }

        m_parsedMapExprToExecStackId[expr] = stackId;

        // The execution stack and its corresponding variable index map are two
        // parallel std::vectors that share their ids. This split helps to
        // achieve some performance improvement.
        m_executionStack.push_back(stack);
        m_stackVariableMap.push_back(variableMap);
        m_state_sizes.push_back(m_state_size);

        return stackId;
    }

    /**
     * @copydoc Interpreter::Evaluate
     */
    NekDouble Evaluate(const int id)
    {
        m_timer.Start();

        ASSERTL1(m_executionStack.size() > id,
                 "unknown analytic expression, it must first be defined "
                 "with DefineFunction(...)");

        ExecutionStack &stack = m_executionStack[id];

        m_state.resize(m_state_sizes[id]);
        for (int i = 0; i < stack.size(); i++)
        {
            (*stack[i]).run_once();
        }

        m_timer.Stop();
        m_total_eval_time += m_timer.TimePerTest(1);

        return m_state[0];
    }

    /**
     * @copydoc Interpreter::Evaluate
     */
    NekDouble Evaluate(const int id, const NekDouble x, const NekDouble y,
                       const NekDouble z, const NekDouble t)
    {
        m_timer.Start();

        ASSERTL1(m_executionStack.size() > id,
                 "unknown analytic expression, it must first be defined with "
                 "DefineFunction(...)");

        ExecutionStack &stack = m_executionStack[id];

        // initialise internal vector of variable values
        m_state.resize(m_state_sizes[id]);

        if (m_variable.size() < 4)
        {
            m_variable.resize(4);
        }

        // no flexibility, no change of variable ordering in m_variable
        // container depending on their names ordering in the input vlist
        // argument of DefineFunction. Ordering convention (x,y,z,t) is assumed.
        m_variable[0] = x;
        m_variable[1] = y;
        m_variable[2] = z;
        m_variable[3] = t;

        // main execution cycle is hidden here
        for (int i = 0; i < stack.size(); i++)
        {
            (*stack[i]).run_once();
        }

        m_timer.Stop();
        m_total_eval_time += m_timer.TimePerTest(1);

        return m_state[0];
    }

    /**
     * @copydoc Interpreter::EvaluateAtPoint
     */
    NekDouble EvaluateAtPoint(const int id, const std::vector<NekDouble> point)
    {
        m_timer.Start();

        ASSERTL1(m_executionStack.size() > id,
                 "unknown analytic expression, it must first be defined with "
                 "DefineFunction(...)");

        ExecutionStack &stack    = m_executionStack[id];
        VariableMap &variableMap = m_stackVariableMap[id];

        ASSERTL1(point.size() == variableMap.size(),
                 "The number of variables used to define this expression should"
                 " match the point dimensionality.");

        // initialise internal vector of variable values
        m_state.resize(m_state_sizes[id]);
        m_variable.resize(point.size());
        VariableMap::const_iterator it;

        for (it = variableMap.begin(); it != variableMap.end(); ++it)
        {
            m_variable[it->second] = point[it->second];
        }

        // main execution cycle is hidden here
        for (int i = 0; i < stack.size(); i++)
        {
            (*stack[i]).run_once();
        }

        m_timer.Stop();
        m_total_eval_time += m_timer.TimePerTest(1);

        return m_state[0];
    }

    /**
     * @copydoc Interpreter::Evaluate
     */
    void Evaluate(const int id, const Array<OneD, const NekDouble> &x,
                  const Array<OneD, const NekDouble> &y,
                  const Array<OneD, const NekDouble> &z,
                  const Array<OneD, const NekDouble> &t,
                  Array<OneD, NekDouble> &result)
    {
        m_timer.Start();
        std::vector<Array<OneD, const NekDouble>> points = {x, y, z, t};
        Evaluate(id, points, result);
        m_timer.Stop();
        m_total_eval_time += m_timer.TimePerTest(1);
    }

    /**
     * @copydoc Interpreter::Evaluate
     */
    void Evaluate(const int id,
                  const std::vector<Array<OneD, const NekDouble>> &points,
                  Array<OneD, NekDouble> &result)
    {
        m_timer.Start();

        const int num_points = points[0].size();
        ASSERTL1(m_executionStack.size() > id,
                 "unknown analytic expression, it must first be defined "
                 "with DefineFunction(...)");
        ASSERTL1(result.size() >= num_points,
                 "destination array must have enough capacity to store "
                 "expression values at each given point");

        ExecutionStack &stack = m_executionStack[id];

        /// If number of points tends to 10^6, one may end up with up to ~0.5Gb
        /// data allocated for m_state only.  Lets split the work into
        /// cache-sized chunks.  Ahtung, magic constant!
        const int max_chunk_size = 1024;
        const int nvals          = points.size();
        const int chunk_size     = (std::min)(max_chunk_size, num_points);

        if (m_state.size() < chunk_size * m_state_sizes[id])
        {
            m_state.resize(m_state_sizes[id] * chunk_size, 0.0);
        }
        if (m_variable.size() < nvals * chunk_size)
        {
            m_variable.resize(nvals * chunk_size, 0.0);
        }
        if (result.size() < num_points)
        {
            result = Array<OneD, NekDouble>(num_points, 0.0);
        }

        int offset    = 0;
        int work_left = num_points;
        while (work_left > 0)
        {
            const int this_chunk_size = (std::min)(work_left, 1024);
            for (int i = 0; i < this_chunk_size; i++)
            {
                for (int j = 0; j < nvals; ++j)
                {
                    m_variable[i + this_chunk_size * j] = points[j][offset + i];
                }
            }
            for (int i = 0; i < stack.size(); i++)
            {
                (*stack[i]).run_many(this_chunk_size);
            }
            for (int i = 0; i < this_chunk_size; i++)
            {
                result[offset + i] = m_state[i];
            }
            work_left -= this_chunk_size;
            offset += this_chunk_size;
        }
        m_timer.Stop();
        m_total_eval_time += m_timer.TimePerTest(1);
    }

    // ======================================================
    //  Private parsing and partial evaluation method
    // ======================================================

    /**
     * @brief Prepares an execution stack for the evaluation of a function.
     *
     * This method prepares the execution stack (an ordered sequence of
     * operators that perform the evaluation) for the parsed evaluation tree.
     * In order to do this, it unrolls the binary tree representing the
     * recursive evaluation into an ordered sequence of commands.  That ordered
     * sequence of commands is equivalent to a bottom-up walk up the evaluation
     * tree, but this allows not to form tree explicitly.
     *
     * This approach requires to introduce explicitly an execution state
     * (memory) shared by commands in the evaluation sequence: recursively
     * dependent commands need to pass data between each other. Such state for
     * the recursive evaluation is passed via return values of a recursive
     * evaluation function --- which is bad if one wants to implement vectorized
     * evaluator.
     *
     * On the other hand, to run through a sequential container of functors is
     * faster than to walk the tree and at each node to check the node type.
     *
     * @param root        Iterator generated by boost::spirit.
     * @param stack       Initially empty sequential container of evaluation
     *                    steps.
     * @param varMap      Maps variable names to their ids.
     * @param stateIndex  An index in the state[] array where an evaluation
     *                    step corresponding to the current tree node
     *                    is allowed to write.
     *
     * @return A std::pair<bool, NekDouble> which encodes fully pre-evaluated
     *         NekDouble values as `(true, value)` if all sub-tree down the
     *         current node evaluates to constant, or flags the opposite
     *         via `(false, 0)`.
     */
    LIB_UTILITIES_EXPORT PrecomputedValue PrepareExecutionAsYouParse(
        const ParsedTreeIterator &location, ExecutionStack &stack,
        VariableMap &variableMap, int stateIndex)
    {
        std::string valueStr(location->value.begin(), location->value.end());
        boost::algorithm::trim(valueStr);

        const bsp::parser_id parserID = location->value.id();
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
        const int num_children = location->children.size();
#endif

        if (parserID == AnalyticExpression::constantID)
        {
            ASSERTL1(num_children == 0,
                     "Illegal children under constant node: " + valueStr);

            auto it = m_constantMapNameToId.find(valueStr);
            ASSERTL1(it != m_constantMapNameToId.end(),
                     "Cannot find the value for the specified constant: " +
                         valueStr);

            return std::make_pair(true, m_constant[it->second]);
        }
        else if (parserID == AnalyticExpression::numberID)
        {
            ASSERTL1(num_children == 0,
                     "Illegal children under number node: " + valueStr);
            return std::make_pair(
                true, boost::lexical_cast<NekDouble>(valueStr.c_str()));
        }
        else if (parserID == AnalyticExpression::variableID)
        {
            ASSERTL1(num_children == 0,
                     "Illegal children under variable node: " + valueStr);

            VariableMap::const_iterator it = variableMap.find(valueStr);
            ASSERTL1(it != variableMap.end(),
                     "Unknown variable parsed: " + valueStr);

            // Variables are not defined at the time of this parse.
            stack.push_back(makeStep<StoreVar>(stateIndex, it->second));
            return std::make_pair(false, 0);
        }
        else if (parserID == AnalyticExpression::parameterID)
        {
            ASSERTL1(num_children == 0,
                     "Illegal children under parameter node: " + valueStr);

            auto it = m_parameterMapNameToId.find(valueStr);
            ASSERTL1(it != m_parameterMapNameToId.end(),
                     "Unknown parameter parsed: " + valueStr);

            // Parameters may change in between of evalutions.
            stack.push_back(makeStep<StorePrm>(stateIndex, it->second));
            return std::make_pair(false, 0);
        }
        else if (parserID == AnalyticExpression::functionID)
        {
            auto it = m_functionMapNameToInstanceType.find(valueStr);
            ASSERTL1(it != m_functionMapNameToInstanceType.end(),
                     "Invalid function specified: " + valueStr);
            ASSERTL1(num_children == 1 || num_children == 2,
                     "Function " + valueStr +
                         " has neither one or two "
                         "arguments: this is not implemented yet.");

            if (location->children.size() == 1)
            {
                PrecomputedValue v = PrepareExecutionAsYouParse(
                    location->children.begin(), stack, variableMap, stateIndex);

                // additive white gaussian noise function
                if (it->second == E_AWGN)
                {
                    int const_index =
                        AddConstant(std::string("SUB_EXPR_") +
                                        std::to_string(m_constant.size()),
                                    v.second);
                    stack.push_back(
                        makeStep<StoreConst>(stateIndex, const_index));
                    stack.push_back(makeStep<EvalAWGN>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                }

                if (true == v.first)
                {
                    return std::make_pair(true,
                                          m_function[it->second](v.second));
                }
            }
            else
            {
                PrecomputedValue v1 =
                    PrepareExecutionAsYouParse(location->children.begin() + 0,
                                               stack, variableMap, stateIndex);
                PrecomputedValue v2 = PrepareExecutionAsYouParse(
                    location->children.begin() + 1, stack, variableMap,
                    stateIndex + 1);
                m_state_size++;

                if (true == v1.first && true == v2.first)
                {
                    return std::make_pair(
                        true, m_function2[it->second](v1.second, v2.second));
                }
            }

            // if somewhere down the parse tree there is a variable or
            // parameter, set up an evaluation sequence.
            switch (it->second)
            {
                case E_ABS:
                    stack.push_back(makeStep<EvalAbs>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_ASIN:
                    stack.push_back(makeStep<EvalAsin>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_ACOS:
                    stack.push_back(makeStep<EvalAcos>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_ATAN:
                    stack.push_back(makeStep<EvalAtan>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_ATAN2:
                    stack.push_back(makeStep<EvalAtan2>(stateIndex, stateIndex,
                                                        stateIndex + 1));
                    return std::make_pair(false, 0);
                case E_ANG:
                    stack.push_back(makeStep<EvalAng>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case E_BESSEL:
                    stack.push_back(makeStep<EvalBessel>(stateIndex, stateIndex,
                                                         stateIndex + 1));
                    return std::make_pair(false, 0);
                case E_CEIL:
                    stack.push_back(makeStep<EvalCeil>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_COS:
                    stack.push_back(makeStep<EvalCos>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_COSH:
                    stack.push_back(makeStep<EvalCosh>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_EXP:
                    stack.push_back(makeStep<EvalExp>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_FABS:
                    stack.push_back(makeStep<EvalFabs>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_FLOOR:
                    stack.push_back(
                        makeStep<EvalFloor>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_LOG:
                    stack.push_back(makeStep<EvalLog>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_LOG10:
                    stack.push_back(
                        makeStep<EvalLog10>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_RAD:
                    stack.push_back(makeStep<EvalRad>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case E_SIN:
                    stack.push_back(makeStep<EvalSin>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_SINH:
                    stack.push_back(makeStep<EvalSinh>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_SQRT:
                    stack.push_back(makeStep<EvalSqrt>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_TAN:
                    stack.push_back(makeStep<EvalTan>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_TANH:
                    stack.push_back(makeStep<EvalTanh>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                case E_SIGN:
                    stack.push_back(makeStep<EvalSign>(stateIndex, stateIndex));
                    return std::make_pair(false, 0);
                default:
                    ASSERTL0(false, "Evaluation of " + valueStr +
                                        " is not implemented yet");
            }
            return std::make_pair(false, 0);
        }
        else if (parserID == AnalyticExpression::unaryID)
        {
            ASSERTL1(*valueStr.begin() == '-',
                     "Illegal factor - it can only be '-' and it was: " +
                         valueStr);
            ASSERTL1(num_children == 1,
                     "Illegal number of children under factor node: " +
                         valueStr);
            PrecomputedValue v = PrepareExecutionAsYouParse(
                location->children.begin(), stack, variableMap, stateIndex);

            // if precomputed value is valid, process it further.
            if (true == v.first)
            {
                return std::make_pair(true, -v.second);
            }
            stack.push_back(makeStep<EvalNeg>(stateIndex, stateIndex));
            return std::make_pair(false, 0);
        }
        else if (parserID == AnalyticExpression::operatorID)
        {
            ASSERTL1(
                num_children == 2,
                "Too few or too many arguments for mathematical operator: " +
                    valueStr);
            PrecomputedValue left = PrepareExecutionAsYouParse(
                location->children.begin() + 0, stack, variableMap, stateIndex);
            PrecomputedValue right =
                PrepareExecutionAsYouParse(location->children.begin() + 1,
                                           stack, variableMap, stateIndex + 1);
            m_state_size++;

            // if both precomputed values are valid, process them further.
            if ((left.first == true) && (right.first == true))
            {
                switch (*valueStr.begin())
                {
                    case '+':
                        return std::make_pair(true, left.second + right.second);
                    case '-':
                        return std::make_pair(true, left.second - right.second);
                    case '*':
                        return std::make_pair(true, left.second * right.second);
                    case '/':
                        return std::make_pair(true, left.second / right.second);
                    case '^':
                        return std::make_pair(
                            true, std::pow(left.second, right.second));
                    case '=':
                        return std::make_pair(true,
                                              left.second == right.second);
                    case '<':
                        if (*(valueStr.end() - 1) == '=')
                        {
                            return std::make_pair(true,
                                                  left.second <= right.second);
                        }
                        else
                        {
                            return std::make_pair(true,
                                                  left.second < right.second);
                        }
                        return std::make_pair(false, 0);
                    case '>':
                        if (*(valueStr.end() - 1) == '=')
                        {
                            return std::make_pair(true,
                                                  left.second >= right.second);
                        }
                        else
                        {
                            return std::make_pair(true,
                                                  left.second > right.second);
                        }
                        return std::make_pair(false, 0);
                    default:
                        ASSERTL0(false,
                                 "Invalid operator encountered: " + valueStr);
                }
                return std::make_pair(false, 0);
            }

            // either operator argument is not fully evaluated
            // add pre-evaluated value to the contaner of constants
            if (true == left.first)
            {
                int const_index =
                    AddConstant(std::string("SUB_EXPR_") +
                                    std::to_string(m_constant.size()),
                                left.second);
                stack.push_back(makeStep<StoreConst>(stateIndex, const_index));
            }
            if (true == right.first)
            {
                int const_index =
                    AddConstant(std::string("SUB_EXPR_") +
                                    std::to_string(m_constant.size()),
                                right.second);
                stack.push_back(
                    makeStep<StoreConst>(stateIndex + 1, const_index));
            }

            switch (*valueStr.begin())
            {
                case '+':
                    stack.push_back(makeStep<EvalSum>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case '-':
                    stack.push_back(makeStep<EvalSub>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case '*':
                    stack.push_back(makeStep<EvalMul>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case '/':
                    stack.push_back(makeStep<EvalDiv>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case '^':
                    stack.push_back(makeStep<EvalPow>(stateIndex, stateIndex,
                                                      stateIndex + 1));
                    return std::make_pair(false, 0);
                case '=':
                    stack.push_back(makeStep<EvalLogicalEqual>(
                        stateIndex, stateIndex, stateIndex + 1));
                    return std::make_pair(false, 0);
                case '<':
                    if (*(valueStr.end() - 1) == '=')
                    {
                        stack.push_back(makeStep<EvalLogicalLeq>(
                            stateIndex, stateIndex, stateIndex + 1));
                    }
                    else
                    {
                        stack.push_back(makeStep<EvalLogicalLess>(
                            stateIndex, stateIndex, stateIndex + 1));
                    }
                    return std::make_pair(false, 0);

                case '>':
                    if (*(valueStr.end() - 1) == '=')
                    {
                        stack.push_back(makeStep<EvalLogicalGeq>(
                            stateIndex, stateIndex, stateIndex + 1));
                    }
                    else
                    {
                        stack.push_back(makeStep<EvalLogicalGreater>(
                            stateIndex, stateIndex, stateIndex + 1));
                    }
                    return std::make_pair(false, 0);

                default:
                    ASSERTL0(false,
                             "Invalid operator encountered: " + valueStr);
            }
            return std::make_pair(false, 0);
        }
        else if (parserID == AnalyticExpression::operatorID)
        {
            ASSERTL1(
                false,
                "Too few or too many arguments for mathematical operator: " +
                    valueStr);
        }
        ASSERTL0(false, "Illegal expression encountered: " + valueStr);
        return std::make_pair(false, 0);
    }

    // ======================================================
    //  Boost::spirit related data structures
    // ======================================================

    /** This is a parser for spirit that parses the CONSTANT values. The default
        constants are those that are in math.h without the M_ prefix and they
       are initialized in the AnalyticExpressionEvaluator constructor. **/

    bsp::symbols<NekDouble> m_constantsParser;

    /** This is the class that is used as the grammar parser for the spirit
     * engine. **/
    class AnalyticExpression : public bsp::grammar<AnalyticExpression>
    {
    private:
        const bsp::symbols<NekDouble> *constants_p;

        /** Variables is a customized parser that will match the variables that
           the function depends on (the first argument of #DefineFunction). **/
        struct variables : bsp::symbols<NekDouble *>
        {
            variables(std::vector<std::string> const &vars)
            {
                for (auto const &it : vars)
                {
                    add(it.c_str(), 0);
                }
            }
        } variables_p;

    public:
        /** These constants are used to determine what parser was used to parse
           what value, which allows for type identification when analyzing the
           parsed AST. **/
        static const int constantID  = 1;
        static const int numberID    = 2;
        static const int variableID  = 3;
        static const int parameterID = 4;
        static const int functionID  = 5;
        static const int unaryID     = 6;
        static const int operatorID  = 7;

        AnalyticExpression(const bsp::symbols<NekDouble> *constants,
                           const std::vector<std::string> &variables)
            : bsp::grammar<AnalyticExpression>(), constants_p(constants),
              variables_p(variables)
        {
        }

        // Trivial constructor to avoid compiler warning with
        // constants_p.
        ~AnalyticExpression()
        {
            constants_p = NULL;
        }

        template <typename ScannerT> struct definition
        {
            /** This function specifies the grammar of the
             * MathAnalyticExpression parser. **/
            definition(AnalyticExpression const &self)
            {
                expression = logical_or;

                logical_or =
                    logical_and >>
                    *((bsp::root_node_d[bsp::str_p("||")] >> logical_and));

                logical_and =
                    equality >>
                    *((bsp::root_node_d[bsp::str_p("&&")] >> equality));

                equality =
                    lt_gt >> *((bsp::root_node_d[bsp::str_p("==")] >> lt_gt));

                lt_gt = add_sub >>
                        *((bsp::root_node_d[bsp::str_p("<=")] >> add_sub) |
                          (bsp::root_node_d[bsp::str_p(">=")] >> add_sub) |
                          (bsp::root_node_d[bsp::ch_p('<')] >> add_sub) |
                          (bsp::root_node_d[bsp::ch_p('>')] >> add_sub));

                add_sub =
                    negate >> *((bsp::root_node_d[bsp::ch_p('+')] >> negate) |
                                (bsp::root_node_d[bsp::ch_p('-')] >> negate));

                negate = !(bsp::root_node_d[bsp::ch_p('-')]) >> mult_div;

                mult_div = exponential >>
                           *((bsp::root_node_d[bsp::ch_p('*')] >> exponential) |
                             (bsp::root_node_d[bsp::ch_p('/')] >> exponential));

                exponential =
                    base >> !(bsp::root_node_d[bsp::ch_p('^')] >> exponential);

                base = number | function | variable | constant | parameter |
                       bsp::inner_node_d[bsp::ch_p('(') >> expression >>
                                         bsp::ch_p(')')];
                parameter =
                    bsp::leaf_node_d[bsp::lexeme_d[(bsp::alpha_p | '_' | '$') >>
                                                   *(bsp::alnum_p | '_' |
                                                     '$')]] >>
                    op;

                function =
                    bsp::root_node_d[functions_p] >>
                    bsp::infix_node_d[bsp::inner_node_d[bsp::ch_p('(') >>
                                                        expression >>
                                                        *(',' >> expression) >>
                                                        bsp::ch_p(')')]];

                variable =
                    bsp::leaf_node_d[bsp::lexeme_d[self.variables_p]] >> op;

                number = bsp::leaf_node_d[bsp::lexeme_d[bsp::real_p]] >> op;

                constant =
                    bsp::leaf_node_d[bsp::lexeme_d[*self.constants_p]] >> op;

                op = bsp::eps_p(bsp::end_p | "||" | "&&" | "==" | "<=" | ">=" |
                                '<' | '>' | '+' | '-' | '*' | '/' | '^' | ')' |
                                ',');
            }

            /** This holds the NekDouble value that is parsed by spirit so it
             * can be stored in the AST. **/
            NekDouble ParsedDouble;

            template <int N>
            using bsp_rule =
                bsp::rule<ScannerT, bsp::parser_context<>, bsp::parser_tag<N>>;

            bsp_rule<constantID> constant;
            bsp_rule<numberID> number;
            bsp_rule<variableID> variable;
            bsp_rule<parameterID> parameter;
            bsp_rule<functionID> function;
            bsp_rule<unaryID> negate;
            bsp_rule<operatorID> base;
            bsp_rule<operatorID> exponent;
            bsp_rule<operatorID> exponential;
            bsp_rule<operatorID> mult_div;
            bsp_rule<operatorID> add_sub;
            bsp_rule<operatorID> lt_gt;
            bsp_rule<operatorID> equality;
            bsp_rule<operatorID> logical_and;
            bsp_rule<operatorID> logical_or;
            bsp_rule<operatorID> expression;
            bsp_rule<operatorID> op;

            bsp_rule<operatorID> const &start() const
            {
                return expression;
            }
        };
    }; // class AnalyticExpression

private:
    // ======================================================
    //  Pre-processed expressions
    // ======================================================

    ///  These vector and map store pre-processed evaluation sequences
    ///  for the analytic expressions. Each ExecutionStack is an ordered
    ///  container of steps of sequential execution process which
    ///  evaluates an analytic expression.

    ExpressionMap m_parsedMapExprToExecStackId;
    std::vector<ExecutionStack> m_executionStack;

    ///  Keeping map of variables individually per each analytic expression
    ///  allows correctly handling expressions which depend on different
    ///  number of variables.

    std::vector<VariableMap> m_stackVariableMap;

    // ======================================================
    //  Execution state and data
    // ======================================================

    ///  The following data structures hold input data to be used on evaluation
    ///  stage. There are three types of input data:
    ///  - constants (never change their value)
    ///  - parameters are allowed to change their values between evaluations
    ///    (compared to constants)
    ///  - variables always change their values at every evaluation call.
    ///  First map looks like <parameter_name, parameter_id> while the second is
    ///  <parameter_id, parameter_value>. The map is used at a preparation
    ///  stage when the analytic expression is parsed. This associates an
    ///  integer id with a parameter name in its string form. On evaluation
    ///  stage the id and a std::vector constant lookup time make evaluation
    ///  faster compared to permanent std::map<std::string, NekDouble> lookup.

    ParameterMap m_parameterMapNameToId;
    ConstantMap m_constantMapNameToId;
    VariableMap m_expressionVariableMap;

    std::vector<NekDouble> m_parameter;
    std::vector<NekDouble> m_constant;
    std::vector<NekDouble> m_variable;

    ///  This vector stores the execution state (memory) used by the
    ///  sequential execution process.
    std::vector<NekDouble> m_state;

    ///  Vector of state sizes per each
    std::vector<int> m_state_sizes;

    ///  This counter is used by PrepareExecutionAsYouParse for finding
    ///  the minimal state size necessary for evaluation of function parsed.
    int m_state_size;

    ///  Timer and sum of evaluation times
    Timer m_timer;
    NekDouble m_total_eval_time;

    boost::mt19937 m_generator;

    // ======================================================
    //  A map of (external) mathematical functions
    // ======================================================

    FunctionNameMap m_functionMapNameToInstanceType;
    std::map<int, OneArgFunc> m_function;
    std::map<int, TwoArgFunc> m_function2;

    // ======================================================
    //  Internal representation of evaluation step
    // ======================================================

    ///  Short names to minimise the infractructural code mess in defining
    ///  functors below.
    typedef std::vector<NekDouble> &vr;
    typedef const std::vector<NekDouble> &cvr;
    typedef const int ci;
    typedef boost::mt19937 &rgt;

    ///  Factory method which makes code little less messy
    template <typename StepType>
    EvaluationStep *makeStep(ci dest, ci src_left = 0, ci src_right = 0)
    {
        return (new StepType(m_generator, m_state, m_constant, m_parameter,
                             m_variable, dest, src_left, src_right));
    }

    enum EvaluationStepType
    {
        E_ABS,
        E_ASIN,
        E_ACOS,
        E_ATAN,
        E_ATAN2,
        E_ANG,
        E_CEIL,
        E_COS,
        E_COSH,
        E_EXP,
        E_FABS,
        E_FLOOR,
        E_LOG,
        E_LOG10,
        E_POW,
        E_RAD,
        E_SIN,
        E_SINH,
        E_SQRT,
        E_TAN,
        E_TANH,
        E_SIGN,
        E_AWGN,
        E_BESSEL
    };

    ///  Function objects (functors)
    struct EvaluationStep
    {
        ///  reference to random number generator
        rgt rng;

        ///  references to arrays holding the common state
        vr state;
        cvr consts;
        cvr params;
        cvr vars;

        ///  indices in the above arrays uniquely defining actual command
        ///  arguments
        ci storeIdx;
        ci argIdx1;
        ci argIdx2;

        EvaluationStep(rgt rn, ci i, ci l, ci r, vr s, cvr c, cvr p, cvr v)
            : rng(rn), state(s), consts(c), params(p), vars(v), storeIdx(i),
              argIdx1(l), argIdx2(r){};

        virtual ~EvaluationStep()
        {
        }

        ///  declaring this guy pure virtual shortens virtual table. It saves
        ///  some execution time.
        virtual void run_many(ci n) = 0;
        virtual void run_once()     = 0;
    };
    struct CopyState : public EvaluationStep
    {
        CopyState(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = state[argIdx1];
        }
        virtual void run_once()
        {
            state[storeIdx] = state[argIdx1];
        }
    };
    struct StoreConst : public EvaluationStep
    {
        StoreConst(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = consts[argIdx1];
        }
        virtual void run_once()
        {
            state[storeIdx] = consts[argIdx1];
        }
    };
    struct StoreVar : public EvaluationStep
    {
        StoreVar(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = vars[argIdx1 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = vars[argIdx1];
        }
    };
    struct StorePrm : public EvaluationStep
    {
        StorePrm(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = params[argIdx1];
        }
        virtual void run_once()
        {
            state[storeIdx] = params[argIdx1];
        }
    };
    struct EvalSum : public EvaluationStep
    {
        EvalSum(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    state[argIdx1 * n + i] + state[argIdx2 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = state[argIdx1] + state[argIdx2];
        }
    };
    struct EvalSub : public EvaluationStep
    {
        EvalSub(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    state[argIdx1 * n + i] - state[argIdx2 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = state[argIdx1] - state[argIdx2];
        }
    };
    struct EvalMul : public EvaluationStep
    {
        EvalMul(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    state[argIdx1 * n + i] * state[argIdx2 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = state[argIdx1] * state[argIdx2];
        }
    };
    struct EvalDiv : public EvaluationStep
    {
        EvalDiv(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    state[argIdx1 * n + i] / state[argIdx2 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = state[argIdx1] / state[argIdx2];
        }
    };
    struct EvalPow : public EvaluationStep
    {
        EvalPow(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    std::pow(state[argIdx1 * n + i], state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::pow(state[argIdx1], state[argIdx2]);
        }
    };
    struct EvalNeg : public EvaluationStep
    {
        EvalNeg(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = -state[argIdx1 * n + i];
        }
        virtual void run_once()
        {
            state[storeIdx] = -state[argIdx1];
        }
    };
    struct EvalLogicalEqual : public EvaluationStep
    {
        EvalLogicalEqual(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    (state[argIdx1 * n + i] == state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = (state[argIdx1] == state[argIdx2]);
        }
    };
    struct EvalLogicalLeq : public EvaluationStep
    {
        EvalLogicalLeq(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    (state[argIdx1 * n + i] <= state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = (state[argIdx1] <= state[argIdx2]);
        }
    };
    struct EvalLogicalLess : public EvaluationStep
    {
        EvalLogicalLess(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    (state[argIdx1 * n + i] < state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = (state[argIdx1] < state[argIdx2]);
        }
    };
    struct EvalLogicalGeq : public EvaluationStep
    {
        EvalLogicalGeq(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    (state[argIdx1 * n + i] >= state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = (state[argIdx1] >= state[argIdx2]);
        }
    };
    struct EvalLogicalGreater : public EvaluationStep
    {
        EvalLogicalGreater(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    (state[argIdx1 * n + i] > state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = (state[argIdx1] > state[argIdx2]);
        }
    };
    struct EvalAbs : public EvaluationStep
    {
        EvalAbs(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::abs(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::abs(state[argIdx1]);
        }
    };
    struct EvalSign : public EvaluationStep
    {
        EvalSign(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = ((state[argIdx1 * n + i] > 0.0) -
                                           (state[argIdx1 * n + i] < 0.0));
        }
        virtual void run_once()
        {
            state[storeIdx] = ((state[argIdx1] > 0.0) - (state[argIdx1] < 0.0));
        }
    };
    struct EvalAsin : public EvaluationStep
    {
        EvalAsin(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::asin(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::asin(state[argIdx1]);
        }
    };
    struct EvalAcos : public EvaluationStep
    {
        EvalAcos(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::acos(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::acos(state[argIdx1]);
        }
    };
    struct EvalAtan : public EvaluationStep
    {
        EvalAtan(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::atan(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::atan(state[argIdx1]);
        }
    };
    struct EvalAtan2 : public EvaluationStep
    {
        EvalAtan2(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    std::atan2(state[argIdx1 * n + i], state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::atan2(state[argIdx1], state[argIdx2]);
        }
    };
    struct EvalAng : public EvaluationStep
    {
        EvalAng(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    ang(state[argIdx1 * n + i], state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = ang(state[argIdx1], state[argIdx2]);
        }
    };
    struct EvalBessel : public EvaluationStep
    {
        EvalBessel(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = boost::math::cyl_bessel_j(
                    state[argIdx1 * n + i], state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] =
                boost::math::cyl_bessel_j(state[argIdx1], state[argIdx2]);
        }
    };
    struct EvalCeil : public EvaluationStep
    {
        EvalCeil(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::ceil(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::ceil(state[argIdx1]);
        }
    };
    struct EvalCos : public EvaluationStep
    {
        EvalCos(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::cos(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::cos(state[argIdx1]);
        }
    };
    struct EvalCosh : public EvaluationStep
    {
        EvalCosh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::cosh(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::cosh(state[argIdx1]);
        }
    };
    struct EvalExp : public EvaluationStep
    {
        EvalExp(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::exp(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::exp(state[argIdx1]);
        }
    };
    struct EvalFabs : public EvaluationStep
    {
        EvalFabs(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::fabs(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::fabs(state[argIdx1]);
        }
    };
    struct EvalFloor : public EvaluationStep
    {
        EvalFloor(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::floor(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::floor(state[argIdx1]);
        }
    };
    struct EvalLog : public EvaluationStep
    {
        EvalLog(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::log(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::log(state[argIdx1]);
        }
    };
    struct EvalLog10 : public EvaluationStep
    {
        EvalLog10(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::log10(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::log10(state[argIdx1]);
        }
    };
    struct EvalRad : public EvaluationStep
    {
        EvalRad(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] =
                    rad(state[argIdx1 * n + i], state[argIdx2 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = rad(state[argIdx1], state[argIdx2]);
        }
    };
    struct EvalSin : public EvaluationStep
    {
        EvalSin(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::sin(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::sin(state[argIdx1]);
        }
    };
    struct EvalSinh : public EvaluationStep
    {
        EvalSinh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::sinh(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::sinh(state[argIdx1]);
        }
    };
    struct EvalSqrt : public EvaluationStep
    {
        EvalSqrt(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::sqrt(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::sqrt(state[argIdx1]);
        }
    };
    struct EvalTan : public EvaluationStep
    {
        EvalTan(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::tan(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::tan(state[argIdx1]);
        }
    };
    struct EvalTanh : public EvaluationStep
    {
        EvalTanh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            for (int i = 0; i < n; i++)
                state[storeIdx * n + i] = std::tanh(state[argIdx1 * n + i]);
        }
        virtual void run_once()
        {
            state[storeIdx] = std::tanh(state[argIdx1]);
        }
    };
    struct EvalAWGN : public EvaluationStep
    {
        EvalAWGN(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r)
            : EvaluationStep(rn, i, l, r, s, c, p, v)
        {
        }
        virtual void run_many(ci n)
        {
            // assuming the argument to AWGN does not depend on spatial
            // variables =>
            boost::variate_generator<boost::mt19937 &,
                                     boost::normal_distribution<>>
                _normal(rng,
                        boost::normal_distribution<>(0, state[storeIdx * n]));
            for (int i = 0; i < n; i++)
            {
                state[storeIdx * n + i] = _normal();
            }
        }
        virtual void run_once()
        {
            boost::variate_generator<boost::mt19937 &,
                                     boost::normal_distribution<>>
                _normal(rng, boost::normal_distribution<>(0, state[storeIdx]));
            state[storeIdx] = _normal();
        }
    };
};

Interpreter::Interpreter()
{
    m_impl = std::unique_ptr<ExpressionEvaluator>(new ExpressionEvaluator());
}

Interpreter::~Interpreter()
{
}

Interpreter::Interpreter(Interpreter &&r) : m_impl(std::move(r.m_impl))
{
}

Interpreter &Interpreter::operator=(Interpreter &&r)
{
    m_impl = std::move(r.m_impl);
    return *this;
}

void Interpreter::SetRandomSeed(unsigned int seed)
{
    m_impl->SetRandomSeed(seed);
}

void Interpreter::AddConstants(
    std::map<std::string, NekDouble> const &constants)
{
    m_impl->AddConstants(constants);
}

int Interpreter::AddConstant(std::string const &name, NekDouble value)
{
    return m_impl->AddConstant(name, value);
}

NekDouble Interpreter::GetConstant(std::string const &name)
{
    return m_impl->GetConstant(name);
}

void Interpreter::SetParameters(std::map<std::string, NekDouble> const &params)
{
    m_impl->SetParameters(params);
}

void Interpreter::SetParameter(std::string const &name, NekDouble value)
{
    m_impl->SetParameter(name, value);
}

NekDouble Interpreter::GetParameter(std::string const &name)
{
    return m_impl->GetParameter(name);
}

NekDouble Interpreter::GetTime() const
{
    return m_impl->GetTime();
}

int Interpreter::DefineFunction(const std::string &vlist,
                                const std::string &function)
{
    return m_impl->DefineFunction(vlist, function);
}

NekDouble Interpreter::Evaluate(const int AnalyticExpression_id)
{
    return m_impl->Evaluate(AnalyticExpression_id);
}

NekDouble Interpreter::Evaluate(const int AnalyticExpression_id,
                                const NekDouble x, const NekDouble y,
                                const NekDouble z, const NekDouble t)
{
    return m_impl->Evaluate(AnalyticExpression_id, x, y, z, t);
}

NekDouble Interpreter::EvaluateAtPoint(const int AnalyticExpression_id,
                                       std::vector<NekDouble> point)
{
    return m_impl->EvaluateAtPoint(AnalyticExpression_id, point);
}

void Interpreter::Evaluate(const int expression_id,
                           const Array<OneD, const NekDouble> &x,
                           const Array<OneD, const NekDouble> &y,
                           const Array<OneD, const NekDouble> &z,
                           const Array<OneD, const NekDouble> &t,
                           Array<OneD, NekDouble> &result)
{
    m_impl->Evaluate(expression_id, x, y, z, t, result);
}

void Interpreter::Evaluate(
    const int expression_id,
    const std::vector<Array<OneD, const NekDouble>> &points,
    Array<OneD, NekDouble> &result)
{
    m_impl->Evaluate(expression_id, points, result);
}

} // namespace LibUtilities
} // namespace Nektar
