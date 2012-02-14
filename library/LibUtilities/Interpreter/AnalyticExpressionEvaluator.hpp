///////////////////////////////////////////////////////////////////////////////
//
// File AnalyticExpressionEvaluator.hpp
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

#ifndef _ANALYTIC_EXPRESSION_EVALUATOR_HPP
#define _ANALYTIC_EXPRESSION_EVALUATOR_HPP

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <boost/version.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>


#if( BOOST_VERSION / 100 % 1000 >= 36 )
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_ast.hpp>
#include <boost/spirit/include/classic_symbols.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_push_back_actor.hpp>

namespace boost_spirit = boost::spirit::classic;
#else
#include <boost/spirit/core.hpp>
#include <boost/spirit/tree/ast.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>

namespace boost_spirit = boost::spirit;
#endif

#include <iostream>
#include <stdarg.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <errno.h>
#include <stdlib.h>

namespace Nektar
{
    namespace LibUtilities
    {

        ///  This class defines evaluator of analytic (symbolic)
        ///  mathematical expressions. Expressions are allowed to
        ///  depend on a number of spatial-time variables and
        ///  parameters. Pre-processing and evaluation stages are
        ///  split. At evaluation stage one specifies values for
        ///  each variable, resulting expression value is returned.
        ///  Vectorized evaluator (evaluate expression at a set of
        ///  points) is available.
        ///
        ///  Internally this class uses boost::spirit to parse analytic
        ///  expressions and unrolls their recursive bracketing structure
        ///  into a sequence of evaluation steps (aka execution stack)
        ///  with resolved data dependencies. Once an expression is
        ///  pre-processed, its execution stack is stored internally
        ///  in order to be re-used.

        class AnalyticExpressionEvaluator
        {

        private:

            struct EvaluationStep;

        public:


            typedef std::map<std::string, int>    VariableMap;
            typedef std::map<std::string, int>    ConstantMap;
            typedef std::map<std::string, int>    ParameterMap;
            typedef std::map<std::string, int>    ExpressionMap;
            typedef std::map<std::string, int>    FunctionNameMap;
            typedef std::vector<EvaluationStep*>  ExecutionStack;
            typedef std::pair<bool, double>       PrecomputedValue;
            typedef double (*OneArgFunc)(double);

            typedef boost_spirit::tree_parse_info<
                            std::string::const_iterator,
                            boost_spirit::node_val_data_factory<double>
                    > ParsedTreeInfo;
            typedef boost_spirit::tree_match<
                            std::string::const_iterator,
                            boost_spirit::node_val_data_factory<double>
                    >::tree_iterator ParsedTreeIterator;


            // ======================================
            //  Methods
            // ======================================


            ///  Initializes the evaluator to a state where it is ready to accept input
            ///  from the #DefineFunction function.
            LIB_UTILITIES_EXPORT AnalyticExpressionEvaluator(void);

            ///  Destroys the execution stack.
            LIB_UTILITIES_EXPORT ~AnalyticExpressionEvaluator(void);


            // ======================================
            //  Setting up methods
            // ======================================


            ///  Constants are evaluated and inserted into the function at the time it is parsed
            ///  when calling the #DefineFunction function. After parsing, if a constant is
            ///  changed, it will not be reflected in the function when Evaluate is called. This
            ///  also means that if a function with an unknown constant is added, and then the
            ///  constant is added, the function will not see the added constant and through an
            ///  exception. This function will add all of the constants in the map argument to
            ///  the global internal constants. If a constant was already loaded previously, it will
            ///  throw an exception stating which constants in the map had this issue. It will add
            ///  all of the constants it can and output the constants it couldn't add in the string
            ///  exception.
            LIB_UTILITIES_EXPORT void AddConstants(std::map<std::string, double> const& constants);

            ///  This function behaves in the same way as #AddConstants, but it only adds one
            ///    constant at a time. If the constant existed previously, an exception will be thrown
            ///    stating the fact. If it did not exist previously, it will be added to the global
            ///    constants and will be used the next time #DefineFunction is called.
            LIB_UTILITIES_EXPORT int AddConstant(std::string const& name, double value);

            ///  If a constant with the specified name exists, it returns the double value that the
            ///  constant stores. If the constant doesn't exist, it throws an exception.
            LIB_UTILITIES_EXPORT double GetConstant(std::string const& name);

            ///  Parameters are like constants, but they are inserted into the function at the time
            ///  #Evaluate is called instead of when the function is parsed. This function can
            ///  be called at any time, and it will take effect in the next call to #Evaluate.
            ///  This function will delete all of the parameters, and replace all of them with only
            ///  the ones in the map argument.
            LIB_UTILITIES_EXPORT void SetParameters(std::map<std::string, double> const& params);

            ///  This function behaves in the same way as #SetParameters, but it only adds one
            ///  parameter and it does not delete the others. If the parameter existed previously,
            ///  it will be overridden and replaced with the new value. If it did not exist previously,
            ///  it will be added to the current parameters.
            LIB_UTILITIES_EXPORT void SetParameter(std::string const& name, double value);

            ///  If a parameter with the specified name exists, it returns the double value that the
            ///  parameter stores. If the parameter doesn't exist, it throws an exception.
            LIB_UTILITIES_EXPORT double GetParameter(std::string const& name);



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
            LIB_UTILITIES_EXPORT double Evaluate0(const int AnalyticExpression_id);

            ///  Evaluation method for expressions depending on 4 variables (+parameters).
            LIB_UTILITIES_EXPORT double Evaluate4(
                        const int AnalyticExpression_id,
                        const double,
                        const double,
                        const double,
                        const double);

            ///  Evaluation method for expressions depending on unspecified number of variables.
            ///  This suitable for expressions depending on more than 4 variables or for the dynamic
            ///  setting some variables as parameters (there is currently no interface method
            ///  for removing a variable from parameter map though).
            LIB_UTILITIES_EXPORT double EvaluateAtPoint(
                        const int AnalyticExpression_id,
                        std::vector<double> point);

            ///  Vectorized evaluation method for expressions depending on 4 variables.
            LIB_UTILITIES_EXPORT Array<OneD, NekDouble> Evaluate4Array(
                        const int expression_id,
                        const Array<OneD, const NekDouble>&,
                        const Array<OneD, const NekDouble>&,
                        const Array<OneD, const NekDouble>&,
                        const Array<OneD, const NekDouble>&);

            ///  Vectorized evaluation method for expressions depending on unspecified
            ///  number of variables.
            LIB_UTILITIES_EXPORT Array<OneD, NekDouble> EvaluateAtPoints(
                    const int expression_id,
                    const std::vector<Array<OneD, const NekDouble> > points);

        private:

            // ======================================================
            //  Private parsing and partial evaluation method
            // ======================================================

            ///  This method prepares the execution stack (an ordered sequence of
            ///  operators that perform the evaluation) for the parsed evaluation tree.
            ///
            ///  In order to do this, it unrolls binary tree representing
            ///  the recursive evaluation into an ordered sequence of commands.
            ///  That ordered sequence of commands is equivalent to bottom-up
            ///  walk up the evaluation tree, but this allows not to form tree explicitly.
            ///
            ///  This approach requires to introduce explicitly an execution state
            ///  (memory) shared by commands in the evaluation sequence: recursively
            ///  dependent commands need to pass data between each other. Such state
            ///  for the recursive evaluation is passed via return values of a recursive
            ///  evaluation function --- which is bad if one wants to implement vectorized
            ///  evaluator.
            ///
            ///  On the other hand, to run through a sequential container of
            ///  functors is faster than to walk the tree and at each node to check
            ///  the node type.
            ///
            ///  \input  root   - iterator generated by boost::spirit;
            ///          stack  - initially empty sequential container of evaluation steps;
            ///          varMap - maps variable names to their ids;
            ///          stateIndex - an index in state[] array where an evaluation
            ///                     step corresponding to the current tree node
            ///                     is allowed to write.
            ///  \output an std::pair<bool, double> which encodes fully pre-evaluated
            ///          double value as pair <true, value> if all sub-tree down the
            ///          current node evaluates to constant, or flags the opposite
            ///          via pair <false,0>.
            LIB_UTILITIES_EXPORT PrecomputedValue PrepareExecutionAsYouParse(
                        const ParsedTreeIterator& root,
                        ExecutionStack& stack,
                        VariableMap &varMap,
                        int stateIndex);


            // ======================================================
            //  Boost::spirit related data structures
            // ======================================================


            /** This is a parser for spirit that parses the CONSTANT values. The default
                constants are those that are in math.h without the M_ prefix and they are
                initialized in the AnalyticExpressionEvaluator constructor. **/

            boost_spirit::symbols<double> m_constantsParser;


            /** This is the class that is used as the grammar parser for the spirit engine. **/
            class AnalyticExpression : public boost_spirit::grammar<AnalyticExpression>
            {
            private:
                const boost_spirit::symbols<double>* constants_p;

                /** Variables is a customized parser that will match the variables that the function
                    depends on (the first argument of #DefineFunction). **/
                struct variables : boost_spirit::symbols<double*>
                {
                    variables(std::vector<std::string> const& vars)
                    {
                        for (std::vector<std::string>::const_iterator it = vars.begin(); it != vars.end(); it++)
                            add(it->c_str(), 0);
                    }
                } variables_p;

            public:
                /** These constants are used to determine what parser was used to parse what value,
                    which allows for type identification when analyzing the parsed AST. **/
                static const int constantID     = 1;
                static const int numberID       = 2;
                static const int variableID     = 3;
                static const int parameterID    = 4;
                static const int functionID     = 5;
                static const int factorID       = 6;
                static const int operatorID     = 7;

                AnalyticExpression(const boost_spirit::symbols<double>* constants, const std::vector<std::string>& variables) :
                        boost_spirit::grammar<AnalyticExpression>(), constants_p(constants), variables_p(variables) {}

                template <typename ScannerT>
                struct definition
                {
                    /** This function specifies the grammar of the MathAnalyticExpression parser. **/
                    definition(AnalyticExpression const& self);

                    /** This holds the double value that is parsed by spirit so it can be stored in the AST. **/
                    double ParsedDouble;

                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<constantID> >     constant;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<numberID> >       number;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<variableID> >     variable;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<parameterID> >    parameter;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<functionID> >     function;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<factorID> >       factor;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     exponential;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     mult_div;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     add_sub;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     lt_gt;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     equality;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     logical_and;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     logical_or;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     expression;
                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     op;

                    boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >     const&
                        start() const { return expression; }
                };
            }; // class AnalyticExpression



            // ======================================================
            //  Pre-processed expressions
            // ======================================================

            ///  These vector and map store pre-processed evaluation sequences
            ///  for the analytic expressions. Each ExecutionStack is an ordered
            ///  container of steps of sequential execution process which
            ///  evaluates an analytic expression.

            ExpressionMap                m_parsedMapExprToExecStackId;
            std::vector<ExecutionStack>  m_executionStack;

            ///  Keeping map of variables individually per each analytic expression
            ///  allows correctly handling expressions which depend on different
            ///  number of variables.

            std::vector<VariableMap>     m_stackVariableMap;

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
            ///  stage when the analytic expression is parsed. This associates an integer
            ///  id with a parameter name in its string form. On evaluation stage the id
            ///  and a std::vector constant lookup time make evaluation faster compared
            ///  to permanent std::map<std::string, double> lookup.

            ParameterMap  m_parameterMapNameToId;
            ConstantMap   m_constantMapNameToId;
            VariableMap   m_expressionVariableMap;

            std::vector<double>  m_parameter;
            std::vector<double>  m_constant;
            std::vector<double>  m_variable;


            ///  This vector stores the execution state (memory) used by the
            ///  sequential execution process.
            std::vector<double>  m_state;

            int    m_state_size;

            ///  Vector of state sizes per each
            std::vector<int>     m_state_sizes;


            // ======================================================
            //  A map of (external) mathematical functions
            // ======================================================

            FunctionNameMap            m_functionMapNameToInstanceType;
            std::map<int, OneArgFunc>  m_function;

            // ======================================================
            //  Internal representation of evaluation step
            // ======================================================

            ///  Factory method which makes code little less messy
            template<typename StepType>
            EvaluationStep* makeStep(int dest, int src_left = 0, int src_right = 0)
            {
                return ( new StepType ( m_state,m_constant,m_parameter,m_variable,dest,src_left,src_right ) );
            }

            enum EvaluationStepType
            {
                    E_ABS,  
                    E_ASIN, 
                    E_ACOS, 
                    E_ATAN, 
                    E_ATAN2,
                    E_CEIL, 
                    E_COS,  
                    E_COSH, 
                    E_EXP,  
                    E_FABS, 
                    E_FLOOR,
                    E_FMOD, 
                    E_LOG,  
                    E_LOG10,
                    E_POW,  
                    E_SIN,  
                    E_SINH, 
                    E_SQRT, 
                    E_TAN,  
                    E_TANH
            };


            ///  Short names to minimise the infractructural code mess in defining functors below.
            typedef std::vector<double>& vr;
            typedef const std::vector<double>& cvr;

            ///  Function objects (functors)
            struct EvaluationStep
            {
                ///  references to arrays holding the common state
                vr       state;
                cvr      consts;
                cvr      params;
                cvr      vars;

                ///  indices in the above arrays uniquely defining
                ///  actual command arguments
                int      storeIdx;
                int      argIdx1;
                int      argIdx2;

                EvaluationStep(int i, int l, int r, vr s, cvr c, cvr p, cvr v):
                    storeIdx(i), argIdx1(l), argIdx2(r), state(s), consts(c), params(p), vars(v) {};

                ///  declaring this guy pure virtual shorten virtual table. It saves some execution time.
                virtual void operator()() = 0;
            };
            struct StoreConst: public EvaluationStep
            {
                StoreConst(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = consts[argIdx1]; }
            };
            struct StoreVar: public EvaluationStep
            {
                StoreVar(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = vars[argIdx1]; } 
            };
            struct StorePrm: public EvaluationStep
            {
                StorePrm(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = params[argIdx1]; } 
            };
            struct EvalSum: public EvaluationStep
            {
                EvalSum(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = state[argIdx1] + state[argIdx2]; }
            };
            struct EvalSub: public EvaluationStep
            {
                EvalSub(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = state[argIdx1] - state[argIdx2]; }
            };
            struct EvalMul: public EvaluationStep
            {
                EvalMul(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = state[argIdx1] * state[argIdx2]; }
            };
            struct EvalDiv: public EvaluationStep
            {
                EvalDiv(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = state[argIdx1] / state[argIdx2]; }
            };
            struct EvalNeg: public EvaluationStep
            {
                EvalNeg(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = - state[argIdx1]; }
            };
            struct EvalLogicalEqual: public EvaluationStep
            {
                EvalLogicalEqual(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = ( state[argIdx1] == state[argIdx2] ); }
            };
            struct EvalLogicalLeq: public EvaluationStep
            {
                EvalLogicalLeq(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = ( state[argIdx1] <= state[argIdx2] ); }
            };
            struct EvalLogicalLess: public EvaluationStep
            {
                EvalLogicalLess(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = ( state[argIdx1] < state[argIdx2] ); }
            };
            struct EvalLogicalGeq: public EvaluationStep
            {
                EvalLogicalGeq(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = ( state[argIdx1] >= state[argIdx2] ); }
            };
            struct EvalLogicalGreater: public EvaluationStep
            {
                EvalLogicalGreater(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = ( state[argIdx1] > state[argIdx2] ); }
            };
            struct EvalAbs: public EvaluationStep
            {
                EvalAbs(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::abs( state[argIdx1] ); }
            };
            struct EvalAsin: public EvaluationStep
            {
                EvalAsin(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::asin( state[argIdx1] ); }
            };
            struct EvalAcos: public EvaluationStep
            {
                EvalAcos(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::acos( state[argIdx1] ); }
            };
            struct EvalAtan: public EvaluationStep
            {
                EvalAtan(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::atan( state[argIdx1] ); }
            };
            struct EvalCeil: public EvaluationStep
            {
                EvalCeil(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::ceil( state[argIdx1] ); }
            };
            struct EvalCos: public EvaluationStep
            {
                EvalCos(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::cos( state[argIdx1] ); }
            };
            struct EvalCosh: public EvaluationStep
            {
                EvalCosh(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::cosh( state[argIdx1] ); }
            };
            struct EvalExp: public EvaluationStep
            {
                EvalExp(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::exp( state[argIdx1] ); }
            };
            struct EvalFabs: public EvaluationStep
            {
                EvalFabs(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::fabs( state[argIdx1] ); }
            };
            struct EvalFloor: public EvaluationStep
            {
                EvalFloor(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::floor( state[argIdx1] ); }
            };
            struct EvalFmod: public EvaluationStep
            {
                EvalFmod(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::fmod( state[argIdx1], state[argIdx2] ); }
            };
            struct EvalLog: public EvaluationStep
            {
                EvalLog(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::log( state[argIdx1] ); }
            };
            struct EvalLog10: public EvaluationStep
            {
                EvalLog10(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::log10( state[argIdx1] ); }
            };
            struct EvalPow: public EvaluationStep
            {
                EvalPow(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::pow( state[argIdx1], state[argIdx2] ); }
            };
            struct EvalSin: public EvaluationStep
            {
                EvalSin(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::sin( state[argIdx1] ); }
            };
            struct EvalSinh: public EvaluationStep
            {
                EvalSinh(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::sinh( state[argIdx1] ); }
            };
            struct EvalSqrt: public EvaluationStep
            {
                EvalSqrt(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::sqrt( state[argIdx1] ); }
            };
            struct EvalTan: public EvaluationStep
            {
                EvalTan(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::tan( state[argIdx1] ); }
            };
            struct EvalTanh: public EvaluationStep
            {
                EvalTanh(vr s, cvr c, cvr p, cvr v, int i, int l, int r): EvaluationStep(i,l,r,s,c,p,v) {}
                virtual void operator()() { state[storeIdx] = std::tanh( state[argIdx1] ); }
            };

        };
    };
};

#endif // _ANALYTIC_EXPRESSION_EVALUATOR_HPP
