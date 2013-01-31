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
#include <LibUtilities/BasicUtils/Timer.h>

#include <boost/version.hpp>
#include <boost/random/mersenne_twister.hpp>  // for mt19937
#include <boost/random/variate_generator.hpp>  // for variate_generator
#include <boost/random/normal_distribution.hpp>


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

#include <string>
#include <vector>
#include <map>
#if defined(__INTEL_COMPILER)
#include <mathimf.h>
#else
#include <cmath>
#endif

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
            typedef std::pair<bool, NekDouble>       PrecomputedValue;
            typedef NekDouble (*OneArgFunc)(NekDouble);

            typedef boost_spirit::tree_parse_info<
                            std::string::const_iterator,
                            boost_spirit::node_val_data_factory<NekDouble>
                    > ParsedTreeInfo;
            typedef boost_spirit::tree_match<
                            std::string::const_iterator,
                            boost_spirit::node_val_data_factory<NekDouble>
                    >::tree_iterator ParsedTreeIterator;

            typedef std::vector<const Array<OneD, const NekDouble>* > VariableArray;

            typedef boost::mt19937  RandomGeneratorType;


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


            LIB_UTILITIES_EXPORT void SetRandomSeed(unsigned int seed = 123u);

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
            LIB_UTILITIES_EXPORT void AddConstants(std::map<std::string, NekDouble> const& constants);

            ///  This function behaves in the same way as #AddConstants, but it only adds one
            ///    constant at a time. If the constant existed previously, an exception will be thrown
            ///    stating the fact. If it did not exist previously, it will be added to the global
            ///    constants and will be used the next time #DefineFunction is called.
            LIB_UTILITIES_EXPORT int AddConstant(std::string const& name, NekDouble value);

            ///  If a constant with the specified name exists, it returns the NekDouble value that the
            ///  constant stores. If the constant doesn't exist, it throws an exception.
            LIB_UTILITIES_EXPORT NekDouble GetConstant(std::string const& name);

            ///  Parameters are like constants, but they are inserted into the function at the time
            ///  #Evaluate is called instead of when the function is parsed. This function can
            ///  be called at any time, and it will take effect in the next call to #Evaluate.
            ///  This function will delete all of the parameters, and replace all of them with only
            ///  the ones in the map argument.
            LIB_UTILITIES_EXPORT void SetParameters(std::map<std::string, NekDouble> const& params);

            ///  This function behaves in the same way as #SetParameters, but it only adds one
            ///  parameter and it does not delete the others. If the parameter existed previously,
            ///  it will be overridden and replaced with the new value. If it did not exist previously,
            ///  it will be added to the current parameters.
            LIB_UTILITIES_EXPORT void SetParameter(std::string const& name, NekDouble value);

            ///  If a parameter with the specified name exists, it returns the NekDouble value that the
            ///  parameter stores. If the parameter doesn't exist, it throws an exception.
            LIB_UTILITIES_EXPORT NekDouble GetParameter(std::string const& name);

            ///  Returns the total time spent in evaluation procedures, seconds.
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
            LIB_UTILITIES_EXPORT void EvaluateAtPoints(
                    const int expression_id,
                    const std::vector<Array<OneD, const NekDouble> > points,
                    Array<OneD, NekDouble>& result);

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
            ///  \output an std::pair<bool, NekDouble> which encodes fully pre-evaluated
            ///          NekDouble value as pair <true, value> if all sub-tree down the
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

            boost_spirit::symbols<NekDouble> m_constantsParser;


            /** This is the class that is used as the grammar parser for the spirit engine. **/
            class AnalyticExpression : public boost_spirit::grammar<AnalyticExpression>
            {
            private:
                const boost_spirit::symbols<NekDouble>* constants_p;

                /** Variables is a customized parser that will match the variables that the function
                    depends on (the first argument of #DefineFunction). **/
                struct variables : boost_spirit::symbols<NekDouble*>
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

                AnalyticExpression(const boost_spirit::symbols<NekDouble>* constants, const std::vector<std::string>& variables) :
                        boost_spirit::grammar<AnalyticExpression>(), constants_p(constants), variables_p(variables) {}

                template <typename ScannerT>
                struct definition
                {
                    /** This function specifies the grammar of the MathAnalyticExpression parser. **/
                    definition(AnalyticExpression const& self);

                    /** This holds the NekDouble value that is parsed by spirit so it can be stored in the AST. **/
                    NekDouble ParsedDouble;

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
            ///  to permanent std::map<std::string, NekDouble> lookup.

            ParameterMap  m_parameterMapNameToId;
            ConstantMap   m_constantMapNameToId;
            VariableMap   m_expressionVariableMap;

            std::vector<NekDouble>  m_parameter;
            std::vector<NekDouble>  m_constant;
            std::vector<NekDouble>  m_variable;


            ///  This vector stores the execution state (memory) used by the
            ///  sequential execution process.
            std::vector<NekDouble>  m_state;

            ///  Vector of state sizes per each
            std::vector<int>     m_state_sizes;

            ///  This counter is used by PrepareExecutionAsYouParse for finding
            ///  the minimal state size necessary for evaluation of function parsed.
            int                  m_state_size;


            ///  Timer and sum of evaluation times
            Timer        m_timer;
            NekDouble       m_total_eval_time;


            RandomGeneratorType  m_generator;
            // boost::variate_generator<RandomGeneratorType&, boost::normal_distribution<> >
            //                     m_normal;



            // ======================================================
            //  A map of (external) mathematical functions
            // ======================================================


            FunctionNameMap            m_functionMapNameToInstanceType;
            std::map<int, OneArgFunc>  m_function;


            // ======================================================
            //  Internal representation of evaluation step
            // ======================================================

            ///  Short names to minimise the infractructural code mess in defining functors below.
            typedef std::vector<NekDouble>& vr;
            typedef const std::vector<NekDouble>& cvr;
            typedef const int  ci;
            typedef RandomGeneratorType& rgt;

            ///  Factory method which makes code little less messy
            template<typename StepType>
            EvaluationStep* makeStep(ci dest, ci src_left = 0, ci src_right = 0)
            {
                return ( new StepType ( m_generator,m_state,m_constant,m_parameter,m_variable,dest,src_left,src_right ) );
            }

            enum EvaluationStepType
            {
                    E_ABS,    E_ASIN,  E_ACOS,  E_ATAN,  E_ATAN2,
                    E_CEIL,   E_COS,   E_COSH,  E_EXP,   E_FABS,
                    E_FLOOR,  E_LOG,   E_LOG10, E_POW,   E_SIN,
                    E_SINH,   E_SQRT,  E_TAN,   E_TANH,  E_SIGN,
                    E_AWGN
            };



            ///  Function objects (functors)
            struct EvaluationStep
            {
                ///  reference to random number generator
                rgt     rng;

                ///  references to arrays holding the common state
                vr      state;
                cvr     consts;
                cvr     params;
                cvr     vars;

                ///  indices in the above arrays uniquely defining actual command arguments
                ci      storeIdx;
                ci      argIdx1;
                ci      argIdx2;

                EvaluationStep(rgt rn, ci i, ci l, ci r, vr s, cvr c, cvr p, cvr v):
                    rng(rn), state(s), consts(c), params(p), vars(v), storeIdx(i), argIdx1(l), argIdx2(r) {};

                virtual ~EvaluationStep() {}

                ///  declaring this guy pure virtual shortens virtual table. It saves some execution time.
                virtual void run_many(ci n) = 0;
                virtual void run_once() = 0;
            };
            struct CopyState: public EvaluationStep
            {
                CopyState(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = state[argIdx1]; }
                virtual void run_once() { state[storeIdx] = state[argIdx1]; }
            };
            struct StoreConst: public EvaluationStep
            {
                StoreConst(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = consts[argIdx1]; }
                virtual void run_once() { state[storeIdx] = consts[argIdx1]; }
            };
            struct StoreVar: public EvaluationStep
            {
                StoreVar(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = vars[argIdx1*n+i]; }
                virtual void run_once() { state[storeIdx] = vars[argIdx1]; }
            };
            struct StorePrm: public EvaluationStep
            {
                StorePrm(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = params[argIdx1]; }
                virtual void run_once() { state[storeIdx] = params[argIdx1]; }
            };
            struct EvalSum: public EvaluationStep
            {
                EvalSum(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = state[argIdx1*n+i] + state[argIdx2*n+i]; }
                virtual void run_once() { state[storeIdx] = state[argIdx1] + state[argIdx2]; }
            };
            struct EvalSub: public EvaluationStep
            {
                EvalSub(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = state[argIdx1*n+i] - state[argIdx2*n+i]; }
                virtual void run_once() { state[storeIdx] = state[argIdx1] - state[argIdx2]; }
            };
            struct EvalMul: public EvaluationStep
            {
                EvalMul(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = state[argIdx1*n+i] * state[argIdx2*n+i]; }
                virtual void run_once() { state[storeIdx] = state[argIdx1] * state[argIdx2]; }
            };
            struct EvalDiv: public EvaluationStep
            {
                EvalDiv(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = state[argIdx1*n+i] / state[argIdx2*n+i]; }
                virtual void run_once() { state[storeIdx] = state[argIdx1] / state[argIdx2]; }
            };
            struct EvalPow: public EvaluationStep
            {
                EvalPow(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::pow( state[argIdx1*n+i], state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::pow( state[argIdx1], state[argIdx2] ); }
            };
            struct EvalNeg: public EvaluationStep
            {
                EvalNeg(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = - state[argIdx1*n+i]; }
                virtual void run_once() { state[storeIdx] = - state[argIdx1]; }
            };
            struct EvalLogicalEqual: public EvaluationStep
            {
                EvalLogicalEqual(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ( state[argIdx1*n+i] == state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = ( state[argIdx1] == state[argIdx2] ); }
            };
            struct EvalLogicalLeq: public EvaluationStep
            {
                EvalLogicalLeq(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ( state[argIdx1*n+i] <= state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = ( state[argIdx1] <= state[argIdx2] ); }
            };
            struct EvalLogicalLess: public EvaluationStep
            {
                EvalLogicalLess(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ( state[argIdx1*n+i] < state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = ( state[argIdx1] < state[argIdx2] ); }
            };
            struct EvalLogicalGeq: public EvaluationStep
            {
                EvalLogicalGeq(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ( state[argIdx1*n+i] >= state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = ( state[argIdx1] >= state[argIdx2] ); }
            };
            struct EvalLogicalGreater: public EvaluationStep
            {
                EvalLogicalGreater(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ( state[argIdx1*n+i] > state[argIdx2*n+i] ); }
                virtual void run_once() { state[storeIdx] = ( state[argIdx1] > state[argIdx2] ); }
            };
            struct EvalAbs: public EvaluationStep
            {
                EvalAbs(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::abs( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::abs( state[argIdx1] ); }
            };
            struct EvalSign: public EvaluationStep
            {
                EvalSign(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = ((state[argIdx1*n+i] > 0.0) - (state[argIdx1*n+i] < 0.0)); }
                virtual void run_once() { state[storeIdx] = ((state[argIdx1] > 0.0) - (state[argIdx1] < 0.0)); }
            };
            struct EvalAsin: public EvaluationStep
            {
                EvalAsin(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::asin( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::asin( state[argIdx1] ); }
            };
            struct EvalAcos: public EvaluationStep
            {
                EvalAcos(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::acos( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::acos( state[argIdx1] ); }
            };
            struct EvalAtan: public EvaluationStep
            {
                EvalAtan(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::atan( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::atan( state[argIdx1] ); }
            };
            struct EvalCeil: public EvaluationStep
            {
                EvalCeil(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::ceil( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::ceil( state[argIdx1] ); }
            };
            struct EvalCos: public EvaluationStep
            {
                EvalCos(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::cos( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::cos( state[argIdx1] ); }
            };
            struct EvalCosh: public EvaluationStep
            {
                EvalCosh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::cosh( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::cosh( state[argIdx1] ); }
            };
            struct EvalExp: public EvaluationStep
            {
                EvalExp(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::exp( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::exp( state[argIdx1] ); }
            };
            struct EvalFabs: public EvaluationStep
            {
                EvalFabs(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::fabs( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::fabs( state[argIdx1] ); }
            };
            struct EvalFloor: public EvaluationStep
            {
                EvalFloor(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::floor( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::floor( state[argIdx1] ); }
            };
            struct EvalLog: public EvaluationStep
            {
                EvalLog(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::log( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::log( state[argIdx1] ); }
            };
            struct EvalLog10: public EvaluationStep
            {
                EvalLog10(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::log10( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::log10( state[argIdx1] ); }
            };
            struct EvalSin: public EvaluationStep
            {
                EvalSin(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::sin( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::sin( state[argIdx1] ); }
            };
            struct EvalSinh: public EvaluationStep
            {
                EvalSinh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::sinh( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::sinh( state[argIdx1] ); }
            };
            struct EvalSqrt: public EvaluationStep
            {
                EvalSqrt(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::sqrt( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::sqrt( state[argIdx1] ); }
            };
            struct EvalTan: public EvaluationStep
            {
                EvalTan(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::tan( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::tan( state[argIdx1] ); }
            };
            struct EvalTanh: public EvaluationStep
            {
                EvalTanh(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n) { for(int i=0;i<n;i++) state[storeIdx*n+i] = std::tanh( state[argIdx1*n+i] ); }
                virtual void run_once() { state[storeIdx] = std::tanh( state[argIdx1] ); }
            };
            struct EvalAWGN: public EvaluationStep
            {
                EvalAWGN(rgt rn, vr s, cvr c, cvr p, cvr v, ci i, ci l, ci r): EvaluationStep(rn,i,l,r,s,c,p,v) {}
                virtual void run_many(ci n)
                {
                    // assuming the argument to AWGN does not depend on spatial variables =>
                    boost::variate_generator<RandomGeneratorType&, boost::normal_distribution<> >
                                       _normal(rng, boost::normal_distribution<>(0, state[storeIdx*n]) );
                    for(int i=0;i<n;i++) { state[storeIdx*n+i] = _normal(); }
                }
                virtual void run_once()
                {
                    boost::variate_generator<RandomGeneratorType&, boost::normal_distribution<> >
                                       _normal(rng, boost::normal_distribution<>(0, state[storeIdx]) );
                    state[storeIdx] = _normal();
                }
            };

        };
    };
};

#endif // _ANALYTIC_EXPRESSION_EVALUATOR_HPP
