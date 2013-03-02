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


#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/random/detail/seed.hpp>

#ifdef _MSC_VER
#include <boost/preprocessor/cat.hpp>
#endif //MSC_VER

#ifdef _MSC_VER
#define NEKTAR_MATH_NAME(x) BOOST_PP_CAT(_, x)
#else
#define NEKTAR_MATH_NAME(x) x
#endif

#if( BOOST_VERSION / 100 % 1000 >= 36 )
using namespace boost::spirit::classic;
#else
using namespace boost::spirit;
#endif

// trying to avoid incompatibility between standart <algorithm> header and
// windows.h header which defines max and min macros.
using std::min;

namespace Nektar
{
    namespace LibUtilities
    {

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
            func(PFD1 p) : func1(p), size(1) {};
            func(PFD2 p) : func2(p), size(2) {};
            func(PFD3 p) : func3(p), size(3) {};
            func(PFD4 p) : func4(p), size(4) {};

            union	// Pointer to a function 
            {
                PFD1 func1;
                PFD2 func2;
                PFD3 func3;
                PFD4 func4;
            };
            size_t size;
        };

        // signum function
        NekDouble sign(NekDouble arg)
        {
            return (arg > 0.0) - (arg < 0.0);
        }

        // Additive white Gaussian noise function.
        // Arg: sigma of the zero-mean gaussian distribution
        // Attention: this function is not actually used for
        // evaluation purposes.
        NekDouble awgn(NekDouble sigma)
        {
            AnalyticExpressionEvaluator::RandomGeneratorType rng;
            boost::variate_generator<
                    AnalyticExpressionEvaluator::RandomGeneratorType&,
                    boost::normal_distribution<>
                >  _normal(rng, boost::normal_distribution<>(0, sigma) );
            return _normal();
        }


        /** This struct creates a parser that matches the function
        definitions from math.h. All of the functions accept one
        of more NekDoubles as arguments and returns a NekDouble. **/
        static struct functions : symbols<func>
        {
            functions()
            {
                // Add all of the functions from math.h 
                add
                    ("abs",		std::abs)
                    ("asin",	asin)
                    ("acos",	acos)
                    ("atan",	atan)
                    ("ceil",	ceil)
                    ("cos",		cos)
                    ("cosh",	cosh)
                    ("exp",		exp)
                    ("fabs",	fabs)	// same as abs
                    ("floor",	floor)
                    ("log",		log)
                    ("log10",	log10)
                    ("sin",		sin)
                    ("sinh",	sinh)
                    ("sqrt",	sqrt)
                    ("tan",		tan)
                    ("tanh",	tanh)
                    // and few more custom functions
                    ("sign",	sign)
                    ("awgn",	awgn)
                    ;
            }
        } functions_p;


        /** This function specifies the grammar of the MathAnalyticExpression parser. **/
        template <typename ScannerT>
        AnalyticExpressionEvaluator::AnalyticExpression::definition<ScannerT>::definition(AnalyticExpression const& self)
        {
            expression	=	logical_or;

            logical_or	=	logical_and >> *(  (root_node_d[str_p("||")] >> logical_and) );

            logical_and	=	equality >> *(  (root_node_d[str_p("&&")] >> equality) );

            equality	=	lt_gt >> *(  (root_node_d[str_p("==")] >> lt_gt) );

            lt_gt		=	add_sub >>
                *(  (root_node_d[str_p("<=")] >> add_sub)
                | (root_node_d[str_p(">=")] >> add_sub)
                | (root_node_d[ch_p('<')] >> add_sub)
                | (root_node_d[ch_p('>')] >> add_sub)
                );

            add_sub		=	mult_div >>
                *(  (root_node_d[ch_p('+')] >> mult_div)
                  | (root_node_d[ch_p('-')] >> mult_div)
                );

            mult_div	=	exponential >>
                *(  (root_node_d[ch_p('*')] >> exponential)
                  | (root_node_d[ch_p('/')] >> exponential)
                );

            exponential	=	factor >>
                *(	(root_node_d[ch_p('^')] >> factor)	);

            factor		=	number
                |	function
                |	variable
                |	constant
                |	parameter
                |	inner_node_d[ch_p('(') >> expression >> ch_p(')')]
            |	(root_node_d[ch_p('-')] >> factor);

            parameter	=	leaf_node_d[ lexeme_d[
                (alpha_p | '_' | '$') >> *(alnum_p | '_' | '$') 
            ] ] >> op;

            function	=	root_node_d[functions_p] >>
                infix_node_d[inner_node_d[ch_p('(') >> expression >> *(',' >> expression) >> ch_p(')')]];

            variable	=	leaf_node_d[ lexeme_d[
                self.variables_p
            ] ] >> op;

            number		=	leaf_node_d[ lexeme_d[
                real_p
            ] ] >> op;

            constant	=	leaf_node_d[ lexeme_d[
                *self.constants_p
            ] ] >> op;

            op = eps_p( end_p | "||" | "&&" | "==" | "<=" | ">=" | '<' | '>' | '+' | '-' | '*' | '/' | '^' | ')' );
        }

        // =========================================================================
        //      AnalyticExpressionEvaluator constructor and setting up methods
        // =========================================================================

        // \brief Initializes the evaluator. Call DefineFunction(...) next.
        AnalyticExpressionEvaluator::AnalyticExpressionEvaluator():
                m_timer(),
                m_total_eval_time(0)
        {
            m_state_size = 1;

            AddConstant("MEANINGLESS", 0.0);
            AddConstant("E",           2.71828182845904523536);     // Natural logarithm
            AddConstant("LOG2E",       1.4426950408889634074);      // log_2 e
            AddConstant("LOG10E",      0.43429448190325182765);     // log_10 e
            AddConstant("LN2",         0.69314718055994530942);     // log_e 2
            AddConstant("LN10",        2.30258509299404568402);     // log_e 10
            AddConstant("PI",          3.14159265358979323846);     // pi
            AddConstant("PI_2",        1.57079632679489661923);     // pi/2
            AddConstant("PI_4",        0.78539816339744830962);     // pi/4
            AddConstant("1_PI",        0.31830988618379067154);     // 1/pi
            AddConstant("2_PI",        0.63661977236758134308);     // 2/pi
            AddConstant("2_SQRTPI",    1.12837916709551257390);     // 2/sqrt(pi)
            AddConstant("SQRT2",       1.41421356237309504880);     // sqrt(2)
            AddConstant("SQRT1_2",     0.70710678118654752440);     // 1/sqrt(2)
            AddConstant("GAMMA",       0.57721566490153286060);     // Euler
            AddConstant("DEG",         57.2957795130823208768);     // deg/radian
            AddConstant("PHI",         1.61803398874989484820);     // golden ratio

            m_functionMapNameToInstanceType["abs"]   =  E_ABS;
            m_functionMapNameToInstanceType["asin"]  =  E_ASIN;
            m_functionMapNameToInstanceType["acos"]  =  E_ACOS;
            m_functionMapNameToInstanceType["atan"]  =  E_ATAN;
            m_functionMapNameToInstanceType["ceil"]  =  E_CEIL;
            m_functionMapNameToInstanceType["cos"]   =  E_COS;
            m_functionMapNameToInstanceType["cosh"]  =  E_COSH;
            m_functionMapNameToInstanceType["exp"]   =  E_EXP;
            m_functionMapNameToInstanceType["fabs"]  =  E_FABS;
            m_functionMapNameToInstanceType["floor"] =  E_FLOOR;
            m_functionMapNameToInstanceType["log"]   =  E_LOG;
            m_functionMapNameToInstanceType["log10"] =  E_LOG10;
            m_functionMapNameToInstanceType["sin"]   =  E_SIN;
            m_functionMapNameToInstanceType["sinh"]  =  E_SINH;
            m_functionMapNameToInstanceType["sqrt"]  =  E_SQRT;
            m_functionMapNameToInstanceType["tan"]   =  E_TAN;
            m_functionMapNameToInstanceType["tanh"]  =  E_TANH;
            m_functionMapNameToInstanceType["sign"]  =  E_SIGN;
            m_functionMapNameToInstanceType["awgn"]  =  E_AWGN;

            m_function[ E_ABS  ] = std::abs;
            m_function[ E_ASIN ] = asin;
            m_function[ E_ACOS ] = acos;
            m_function[ E_ATAN ] = atan;
            m_function[ E_CEIL ] = ceil;
            m_function[ E_COS  ] = cos;
            m_function[ E_COSH ] = cosh;
            m_function[ E_EXP  ] = exp;
            m_function[ E_FABS ] = fabs;
            m_function[ E_FLOOR] = floor;
            m_function[ E_LOG  ] = log;
            m_function[ E_LOG10] = log10;
            m_function[ E_SIN  ] = sin;
            m_function[ E_SINH ] = sinh;
            m_function[ E_SQRT ] = sqrt;
            m_function[ E_TAN  ] = tan;
            m_function[ E_TANH ] = tanh;
            m_function[ E_SIGN ] = sign;
            // there is no entry to m_function that correspond to awgn function.
            // this is made in purpose. This function need not be pre-evaluated once!
        }


        AnalyticExpressionEvaluator::~AnalyticExpressionEvaluator(void)
        {
            for (std::vector<ExecutionStack>::iterator it_es = m_executionStack.begin(); it_es != m_executionStack.end(); ++it_es)
            {
                for (std::vector<EvaluationStep*>::iterator it = (*it_es).begin(); it != (*it_es).end(); ++it)
                {
                    delete *it;
                }
                (*it_es).clear();
            }
            m_executionStack.clear();
        }


        void AnalyticExpressionEvaluator::SetRandomSeed(unsigned int seed)
        {
            m_generator.seed(seed);
        }


        void AnalyticExpressionEvaluator::AddConstants(std::map<std::string, NekDouble> const& constants)
        {
            for (std::map<std::string, NekDouble>::const_iterator it = constants.begin(); it != constants.end(); ++it)
            {
                AddConstant(it->first, it->second);
            }
        }

        int AnalyticExpressionEvaluator::AddConstant(std::string const& name, NekDouble value)
        {
            ConstantMap::const_iterator it = m_constantMapNameToId.find(name);
            if (it == m_constantMapNameToId.end())
            {
                // we are trying to avoid duplicating entries in m_constantParser and m_constants
                m_constantsParser.add(name.c_str(), value);
                int index = m_constant.size();
                m_constantMapNameToId[name] = index;
                m_constant.push_back(value);
                return index;
            }
            else
            {
	      if(m_constant[it->second] != value)
		  {
		    std::string errormsg("Attempt to add numerically different constants under the same name: ");
		    errormsg += name; 
		    std::cout << errormsg << std::endl;
		  }
	    //ASSERTL1(m_constant[it->second] == value, "Attempt to add numerically different constants under the same name: " + name);
            }
            return it->second;
        }

        NekDouble AnalyticExpressionEvaluator::GetConstant(std::string const& name)
        {
            NekDouble* value = find(m_constantsParser, name.c_str());

            ASSERTL1(value != NULL, "Constant variable not found: " + name);

            return *value;
        }

        void AnalyticExpressionEvaluator::SetParameters(std::map<std::string, NekDouble> const& params)
        {
            for (std::map<std::string, NekDouble>::const_iterator it = params.begin(); it != params.end(); it++)
            {
                SetParameter(it->first, it->second);
            }
        }

        void AnalyticExpressionEvaluator::SetParameter(std::string const& name, NekDouble value)
        {
            ParameterMap::const_iterator it = m_parameterMapNameToId.find(name);
            if (it == m_parameterMapNameToId.end())
            {
                m_parameterMapNameToId[name] = m_parameter.size();
                m_parameter.push_back(value);
            }
            else
            {
                // if parameter is known, change its value
                m_parameter[ it->second ] = value;
            }
        }


        NekDouble AnalyticExpressionEvaluator::GetParameter(std::string const& name)
        {
            ParameterMap::const_iterator it = m_parameterMapNameToId.find(name);

            ASSERTL1(it != m_parameterMapNameToId.end(), "Parameter not found: " + name);

            return m_parameter[ it->second ];
        }


        NekDouble AnalyticExpressionEvaluator::GetTime() const
        {
            return m_total_eval_time;
        }


        // ======================================================
        //  Public evaluate methods
        // ======================================================

        int AnalyticExpressionEvaluator::DefineFunction(const std::string& vlist, const std::string& expr)
        {
            // Find the previous parsing.
            ExpressionMap::const_iterator it = m_parsedMapExprToExecStackId.find(expr);
            if (it != m_parsedMapExprToExecStackId.end())
            {
                // if this function is already defined, don't do anything but
                // return its ID.
                return it->second;
            }

            // ----------------------------------------------
            // Prepare an iterator that allows to walk along
            // the string representing an analytic expression in the order
            // that respects its recursive structure (thanks to boost::spirit).
            // ----------------------------------------------

            // Parse the vlist input and separate the variables into ordered entries
            // in a vector<string> object. These need to be ordered because this is
            // the order the variables will get assigned to in the Map when Evaluate(...)
            // is called.
            std::vector<std::string> variableNames;
            parse((char*) vlist.c_str(), ( *space_p >>
                       *(
                                +(+graph_p)[push_back_a(variableNames)]
                                    >> +space_p
                             )
                        )
            );
            // Set up our grammar
            AnalyticExpression myGrammar(&m_constantsParser, variableNames);

            // Do the actual parsing with boost::spirit and alert the user if there was an error with an exception.
            ParsedTreeInfo   parseInfo = ast_parse<
                                                node_val_data_factory<NekDouble>,
                                                std::string::const_iterator,
                                                AnalyticExpression,
                                                space_parser
                                             >
                                             (expr.begin(), expr.end(), myGrammar, space_p);

            ASSERTL1(parseInfo.full != false, "Unable to fully parse function. Stopped just before: "
                                         + std::string(parseInfo.stop, parseInfo.stop + 15));

            // ----------------------------------------------
            // Data parsed, start setting up internal data structures.
            // ----------------------------------------------

            ExecutionStack  stack;
            VariableMap     variableMap;

            int stackId = m_executionStack.size();
            m_state_size = 1;

            // register all variables declared in the expression
            for (int i = 0; i < variableNames.size(); i++)
            {
                variableMap[variableNames[i]] = i;
            }

            // then prepare an execution stack.
            // this method also calculates a length of internal
            // state storage (m_state_size) for this function.
            PrecomputedValue v = PrepareExecutionAsYouParse(parseInfo.trees.begin(), stack, variableMap, 0);

            // constant expression, fully evaluated
            if (true == v.first)
            {
                ASSERTL1(stack.size() == 0, "Constant expression yeilds non-empty execution stack. Bug in PrepareExecutionAsYouParse()");

                int const_index = AddConstant(std::string("EXPRESSION_") + boost::lexical_cast<std::string>(stackId), v.second);
                stack.push_back ( makeStep<StoreConst>( 0, const_index ) );
            }

            m_parsedMapExprToExecStackId[expr] = stackId;

            // the execution stack and its corresponding variable index map are
            // two parallel std::vectors that share their ids. This split helps
            // to achieve some performance improvement.
            m_executionStack.push_back(stack);
            m_stackVariableMap.push_back(variableMap);
            m_state_sizes.push_back(m_state_size);
            return stackId;
        }


        NekDouble AnalyticExpressionEvaluator::Evaluate(const int expression_id)
        {
            m_timer.Start();

            ASSERTL1(m_executionStack.size() > expression_id, "unknown analytic expression, it must first be defined with DefineFunction(...)");

            ExecutionStack &stack = m_executionStack[expression_id];

            m_state.resize(m_state_sizes[expression_id]);
            for (int i = 0; i < stack.size(); i++)
            {
                (*stack[i]).run_once();
            }

            m_timer.Stop();
            m_total_eval_time += m_timer.TimePerTest(1);

            return m_state[0];
        }

        NekDouble AnalyticExpressionEvaluator::Evaluate(
                const int expression_id,
                const NekDouble x,
                const NekDouble y,
                const NekDouble z,
                const NekDouble t)
        {
            m_timer.Start();

            ASSERTL1(m_executionStack.size() > expression_id, "unknown analytic expression, it must first be defined with DefineFunction(...)");

            ExecutionStack &stack = m_executionStack[expression_id];

            // initialise internal vector of variable values
            m_state.resize(m_state_sizes[expression_id]);

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

        NekDouble AnalyticExpressionEvaluator::EvaluateAtPoint(const int expression_id, const std::vector<NekDouble> point)
        {
            m_timer.Start();

            ASSERTL1(m_executionStack.size() > expression_id, "unknown analytic expression, it must first be defined with DefineFunction(...)");

            ExecutionStack&  stack    = m_executionStack[expression_id];
            VariableMap&  variableMap = m_stackVariableMap[expression_id];

            ASSERTL1(point.size() == variableMap.size(), "The number of variables used to define this expression should match the point dimensionality.");

            // initialise internal vector of variable values
            m_state.resize(m_state_sizes[expression_id]);
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


       void AnalyticExpressionEvaluator::Evaluate(
                    const int expression_id,
                    const Array<OneD, const NekDouble>& x,
                    const Array<OneD, const NekDouble>& y,
                    const Array<OneD, const NekDouble>& z,
                    const Array<OneD, const NekDouble>& t,
                    Array<OneD, NekDouble>& result)
        {
            m_timer.Start();

            const int num_points = x.num_elements();
            ASSERTL1(m_executionStack.size() > expression_id, "unknown analytic expression, it must first be defined with DefineFunction(...)");
            ASSERTL1(result.num_elements() >= num_points, "destination array must have enough capacity to store expression values at each given point");

            ExecutionStack &stack = m_executionStack[expression_id];

            /// If number of points tends to 10^6, one may end up
            /// with up to ~0.5Gb data allocated for m_state only.
            /// Lets split the work into cache-sized chunks.
            /// Ahtung, magic constant!
            const int max_chunk_size = 1024;

            /// please don't remove brackets around std::min, it screws up windows compilation
            const int chunk_size = (std::min)(max_chunk_size, num_points);
            if (m_state.size() < chunk_size * m_state_sizes[expression_id] )
            {
                m_state.resize( m_state_sizes[expression_id] * chunk_size, 0.0 );
            }
            if (m_variable.size() < 4 * chunk_size )
            {
                m_variable.resize( 4 * chunk_size, 0.0);
            }
            if (result.num_elements() < num_points)
            {
                result = Array<OneD, NekDouble>(num_points, 0.0);
            }

            int offset = 0;
            int work_left = num_points;
            while(work_left > 0)
            {
                const int this_chunk_size = (std::min)(work_left, 1024);
                for (int i = 0; i < this_chunk_size; i++)
                {
                    m_variable[i+this_chunk_size*0] = x[offset + i];
                    m_variable[i+this_chunk_size*1] = y[offset + i];
                    m_variable[i+this_chunk_size*2] = z[offset + i];
                    m_variable[i+this_chunk_size*3] = t[offset + i];
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
                offset    += this_chunk_size;
            }
            m_timer.Stop();
            m_total_eval_time += m_timer.TimePerTest(1);
        }


        void AnalyticExpressionEvaluator::EvaluateAtPoints(
                    const int expression_id,
                    const std::vector<Array<OneD, const NekDouble> > points,
                    Array<OneD, NekDouble>& result)
        {
            m_timer.Start();

            /// \todo test this function properly/update as the method above

            ASSERTL1(m_executionStack.size() > expression_id, "unknown analytic expression, it must first be defined with DefineFunction(...)");

            ExecutionStack&  stack    = m_executionStack[expression_id];
            VariableMap&  variableMap = m_stackVariableMap[expression_id];

            const int num = points[0].num_elements();
            m_state.resize(m_state_sizes[expression_id]*num);

            // assuming all points have same # of coordinates
            m_variable.resize(4*num,0.0);

            for (int i = 0; i < points.size(); i++)
            {
                for (VariableMap::const_iterator it = variableMap.begin(); it != variableMap.end(); ++it)
                {
                    m_variable[it->second] = points[i][it->second];
                }
            }
            for (int j = 0; j < stack.size(); j++)
            {
                (*stack[j]).run_many(num);
            }
            for (int i = 0; i < num; ++i)
            {
                result[i] = m_state[i];
            }

            m_timer.Stop();
            m_total_eval_time += m_timer.TimePerTest(1);
        }



        AnalyticExpressionEvaluator::PrecomputedValue AnalyticExpressionEvaluator::PrepareExecutionAsYouParse(
                    const ParsedTreeIterator& location,
                    ExecutionStack& stack,
                    VariableMap& variableMap,
                    int stateIndex)
        {
            std::string valueStr(location->value.begin(), location->value.end());
            boost::algorithm::trim(valueStr);

            const parser_id parserID  = location->value.id();
#if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
            const int num_children    = location->children.size();
#endif

            if (parserID == AnalyticExpression::constantID)
            {
                ASSERTL1(num_children == 0, "Illegal children under constant node: " + valueStr);

                ConstantMap::const_iterator it = m_constantMapNameToId.find(valueStr);
                ASSERTL1(it != m_constantMapNameToId.end(), "Cannot find the value for the specified constant: " + valueStr);

                return std::make_pair(true, m_constant[it->second]);
            }
            else if (parserID == AnalyticExpression::numberID)
            {
                ASSERTL1(num_children == 0, "Illegal children under number node: " + valueStr);
                return std::make_pair(true, boost::lexical_cast<NekDouble>(valueStr.c_str()) );
            }
            else if (parserID == AnalyticExpression::variableID)
            {
                ASSERTL1(num_children == 0, "Illegal children under variable node: " + valueStr);

                VariableMap::const_iterator it = variableMap.find(valueStr);
                ASSERTL1(it != variableMap.end(), "Unknown variable parsed: " + valueStr);

                // Variables are not defined at the time of this parse.
                stack.push_back ( makeStep<StoreVar>( stateIndex, it->second ) );
                return std::make_pair(false, 0);
            }
            else if (parserID == AnalyticExpression::parameterID)
            {
                ASSERTL1(num_children == 0, "Illegal children under parameter node: " + valueStr);

                ParameterMap::const_iterator it = m_parameterMapNameToId.find(valueStr);
                ASSERTL1(it != m_parameterMapNameToId.end(), "Unknown parameter parsed: " + valueStr);

                // Parameters may change in between of evalutions.
                stack.push_back ( makeStep<StorePrm>( stateIndex, it->second ) );
                return std::make_pair(false, 0);
            }
            else if (parserID == AnalyticExpression::functionID)
            {
                FunctionNameMap::const_iterator it = m_functionMapNameToInstanceType.find(valueStr);
                ASSERTL1(it != m_functionMapNameToInstanceType.end(), "Invalid function specified: " + valueStr);
                ASSERTL1(num_children == 1, "Function " + valueStr + " would like to have too few or too many arguments. This is not implemented yet");

                PrecomputedValue v = PrepareExecutionAsYouParse(location->children.begin(), stack, variableMap, stateIndex);

                // additive white gaussian noise function
                if (it->second == E_AWGN)
                {
                    int const_index = AddConstant(std::string("SUB_EXPR_") + boost::lexical_cast<std::string>(m_constant.size()), v.second);
                    stack.push_back ( makeStep<StoreConst>( stateIndex, const_index ) );
                    stack.push_back ( makeStep<EvalAWGN>( stateIndex, stateIndex ) );
                    return std::make_pair(false,0);
                }

                // if precomputed value is valid, return function(value).
                if (true == v.first)
                {
                    return std::make_pair( true, m_function[it->second](v.second) );
                }


                // if somewhere down the parse tree there is a variable or parameter, set up an
                // evaluation sequence.
                switch (it->second)
                {
                    case E_ABS:
                        stack.push_back ( makeStep<EvalAbs>( stateIndex, stateIndex) );
                        return std::make_pair(false,0);
                    case E_ASIN:
                        stack.push_back ( makeStep<EvalAsin>( stateIndex, stateIndex) );
                        return std::make_pair(false,0);
                    case E_ACOS:
                        stack.push_back ( makeStep<EvalAcos>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_ATAN:
                        stack.push_back ( makeStep<EvalAtan>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_CEIL:
                        stack.push_back ( makeStep<EvalCeil>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_COS:
                        stack.push_back ( makeStep<EvalCos>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_COSH:
                        stack.push_back ( makeStep<EvalCosh>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_EXP: 
                        stack.push_back ( makeStep<EvalExp>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_FABS:
                        stack.push_back ( makeStep<EvalFabs>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_FLOOR:
                        stack.push_back ( makeStep<EvalFloor>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_LOG: 
                        stack.push_back ( makeStep<EvalLog>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_LOG10:
                        stack.push_back ( makeStep<EvalLog10>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_SIN: 
                        stack.push_back ( makeStep<EvalSin>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_SINH:
                        stack.push_back ( makeStep<EvalSinh>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_SQRT:
                        stack.push_back ( makeStep<EvalSqrt>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_TAN: 
                        stack.push_back ( makeStep<EvalTan>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_TANH:
                        stack.push_back ( makeStep<EvalTanh>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    case E_SIGN:
                        stack.push_back ( makeStep<EvalSign>( stateIndex, stateIndex ) );
                        return std::make_pair(false,0);
                    default:
                        ASSERTL0(false, "Evaluation of " + valueStr + " is not implemented yet");
                }
                return std::make_pair(false,0);
            }
            else if (parserID == AnalyticExpression::factorID)
            {
                ASSERTL1(*valueStr.begin() == '-', "Illegal factor - it can only be '-' and it was: " + valueStr);
                ASSERTL1(num_children == 1, "Illegal number of children under factor node: " + valueStr);

                PrecomputedValue v = PrepareExecutionAsYouParse(location->children.begin(), stack, variableMap, stateIndex);

                // if precomputed value is valid, process it further.
                if (true == v.first)
                {
                    return std::make_pair( true, - v.second );
                }
                stack.push_back (makeStep<EvalNeg>( stateIndex, stateIndex) );
                return std::make_pair(false,0);
            }
            else if (parserID == AnalyticExpression::operatorID)
            {
                ASSERTL1(num_children == 2, "Too few or too many arguments for mathematical operator: " + valueStr);
                PrecomputedValue left  = PrepareExecutionAsYouParse(location->children.begin()+0, stack, variableMap, stateIndex);
                PrecomputedValue right = PrepareExecutionAsYouParse(location->children.begin()+1, stack, variableMap, stateIndex+1);
                m_state_size++;

                // if both precomputed values are valid, process them further.
                if ((true == left.first) && (true == right.first))
                {
                    switch(*valueStr.begin())
                    {
                    case '+':
                        return std::make_pair( true, left.second + right.second );
                    case '-':
                        return std::make_pair( true, left.second - right.second );
                    case '*':
                        return std::make_pair( true, left.second * right.second );
                    case '/':
                        return std::make_pair( true, left.second / right.second );
                    case '^':
                        return std::make_pair( true, std::pow(left.second, right.second) );
                    case '=':
                        return std::make_pair( true, left.second == right.second );
                    case '<':
                        if (*(valueStr.end()-1) == '=')
                        {
                            return std::make_pair( true, left.second <= right.second );
                        }
                        else
                        {
                            return std::make_pair( true, left.second < right.second );
                        }
                        return std::make_pair(false,0);
                    case '>':
                        if (*(valueStr.end()-1) == '=')
                        {
                            return std::make_pair( true, left.second >= right.second );
                        }
                        else
                        {
                            return std::make_pair( true, left.second > right.second );
                        }
                        return std::make_pair(false,0);
                    default:
                        ASSERTL0(false, "Invalid operator encountered: " + valueStr);
                    }
                    return std::make_pair(false,0);
                }

                // either operator argument is not fully evaluated
                // add pre-evaluated value to the contaner of constants
                if (true == left.first)
                {
                    int const_index = AddConstant(std::string("SUB_EXPR_") + boost::lexical_cast<std::string>(m_constant.size()), left.second);
                    stack.push_back ( makeStep<StoreConst>( stateIndex, const_index ) );
                }
                if (true == right.first)
                {
                    int const_index = AddConstant(std::string("SUB_EXPR_") + boost::lexical_cast<std::string>(m_constant.size()), right.second);
                    stack.push_back ( makeStep<StoreConst>( stateIndex+1, const_index ) );
                }


                switch(*valueStr.begin())
                {
                case '+':
                    stack.push_back (makeStep<EvalSum>( stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '-':
                    stack.push_back (makeStep<EvalSub> (stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '*':
                    stack.push_back (makeStep<EvalMul>( stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '/':
                    stack.push_back (makeStep<EvalDiv>( stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '^':
                    stack.push_back (makeStep<EvalPow>( stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '=':
                    stack.push_back (makeStep<EvalLogicalEqual>( stateIndex, stateIndex, stateIndex+1 ) );
                    return std::make_pair(false,0);
                case '<':
                    if (*(valueStr.end()-1) == '=')
                    {
                        stack.push_back (makeStep<EvalLogicalLeq>( stateIndex, stateIndex, stateIndex+1 ) );
                    }
                    else
                    {
                        stack.push_back (makeStep<EvalLogicalLess>( stateIndex, stateIndex, stateIndex+1 ) );
                    }
                    return std::make_pair(false,0);

                case '>':
                    if (*(valueStr.end()-1) == '=')
                    {
                        stack.push_back (makeStep<EvalLogicalGeq>( stateIndex, stateIndex, stateIndex+1 ) );
                    }
                    else
                    {
                        stack.push_back (makeStep<EvalLogicalGreater>( stateIndex, stateIndex, stateIndex+1 ) );
                    }
                    return std::make_pair(false,0);

                default:
                    ASSERTL0(false, "Invalid operator encountered: " + valueStr);
                }
                return std::make_pair(false,0);
            }
            ASSERTL0(false, "Illegal expression encountered: " + valueStr);
            return std::make_pair(false,0);
        }

    };
};
