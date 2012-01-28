// Expression Evaluator
//
// File:		ExpressionEvaluator.h
// Description:	Implementation of an expression parser.
// Author:		Michael DeLisi
// Date:		July 27, 2007

#ifndef _EXPRESSION_EVALUATOR_H
#define _EXPRESSION_EVALUATOR_H

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <boost/version.hpp>

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
#include <math.h>
#include <errno.h>
#include <stdlib.h>

// Microsoft fix for math functions
#ifndef _MSC_VER 
#define _hypot hypot
#define _jn jn
#define _j0 j0
#define _j1 j1
#define _yn yn
#define _y0 y0
#define _y1 y1
#endif

namespace Nektar
{
    namespace LibUtilities
    {
        class ExpressionEvaluator
        {
        public:
			/** Initializes the evaluator to a state where it is ready to accept input
			    from the #DefineFunction function. **/
	        LIB_UTILITIES_EXPORT ExpressionEvaluator(void);

	        /** Destroys the evaluator object and deletes all of the saves ASTs. **/
	        LIB_UTILITIES_EXPORT ~ExpressionEvaluator(void);

	        /** This function allows one to define a function to evaluate. The first argument (vlist)
	            is a list of variables (separated by spaces) that the second argument (function)
	            depends on. For example, if function = "x + y", then vlist should most likely be
			    "x y", unless you are defining x or y as parameters with #SetParameters.
	        
	            If the function definition is changed, all of the current values such as parsed AST and
	            parameters are saved if you call DefineFunction with the same "function" string. Also
	            note that if you have previously set a function to use a certain "vlist", and then
	            you call DefineFunction again with the same function string but a different vlist string,
	            it will use the old vlist string from the first declaration. **/
	        LIB_UTILITIES_EXPORT int DefineFunction(const std::string& vlist, const std::string& function);

	        /** Constants are evaluated and inserted into the function at the time it is parsed
			    when calling the #DefineFunction function. After parsing, if a constant is
	            changed, it will not be reflected in the function when Evaluate is called. This
				also means that if a function with an unknown constant is added, and then the
				constant is added, the function will not see the added constant and through an
				exception. This function will add all of the constants in the map argument to
				the global internal constants. If a constant was already loaded previously, it will
				throw an exception stating which constants in the map had this issue. It will add
				all of the constants it can and output the constants it couldn't add in the string
				exception. **/
	        LIB_UTILITIES_EXPORT void AddConstants(std::map<std::string, double> const& constants);

			/** This function behaves in the same way as #AddConstants, but it only adds one
				constant at a time. If the constant existed previously, an exception will be thrown
				stating the fact. If it did not exist previously, it will be added to the global
				constants and will be used the next time #DefineFunction is called. **/
	        LIB_UTILITIES_EXPORT void AddConstant(std::string const& name, double value);

	        /** If a constant with the specified name exists, it returns the double value that the
	            constant stores. If the constant doesn't exist, it throws an exception. **/
	        LIB_UTILITIES_EXPORT double GetConstant(std::string const& name);

	        /** Parameters are like constants, but they are inserted into the function at the time
				#Evaluate is called instead of when the function is parsed. This function can
				be called at any time, and it will take effect in the next call to #Evaluate.
				This function will delete all of the parameters, and replace all of them with only
				the ones in the map argument. **/
	        LIB_UTILITIES_EXPORT void SetParameters(std::map<std::string, double> const& params);

			/** This function behaves in the same way as #SetParameters, but it only adds one
				parameter and it does not delete the others. If the parameter existed previously,
				it will be overridden and replaced with the new value. If it did not exist previously,
				it will be added to the current parameters. **/
	        LIB_UTILITIES_EXPORT void SetParameter(std::string const& name, double value);

	        /** If a parameter with the specified name exists, it returns the double value that the
				parameter stores. If the parameter doesn't exist, it throws an exception. **/
	        LIB_UTILITIES_EXPORT double GetParameter(std::string const& name);

            /** This function evaluates at one point the expression defined via
                #DefineFunction. If #DefineFunction was called multiple times for
                different expressions, this function evaluates the latest expression
                defined via #DefineFunction.
                The arguments to the function are the values that were defined in the
                first argument of #DefineFunction in the same order as they were
                listed. It results the result of the evaluation as a double. If a function
                is not currently defined, behavior may be unpredictable since I do not
                know how many arguments were passed in, so it will probably not be cleaned
                up correctly. **/
            LIB_UTILITIES_EXPORT double Evaluate(double start = 0, ...);

            /** This function works the same as above but takes additional argument ---
                the string form of evaluated expression --- in order to choose which
                expression (among others defined via #DefineFunction) to evaluate.
                \attention This version is slower than the version below due to a
                logarithmic time hash map lookups with std::string keys.  **/
            LIB_UTILITIES_EXPORT double Evaluate(std::string const& function, double start = 0, ...);

            /** This function takes an unique ID of previously defined expression. It performs
                exactly the same operation as the other two functions above.  **/
            LIB_UTILITIES_EXPORT double Evaluate(const int expression_id, double start = 0, ...);


            /** This function evaluates the function defined from #DefineFunction. It accepts
                a "vector<double> const*" argument for each of the variables given in the first
                argument (vlist) of the #DefineFunction function. The vectors for each variable
                must be in the same order as they appear in the vlist argument. This function will
                calculate the function from #DefineFunction with values from each of the vectors
                and then store the result in another vector. If the vectors aren't the same length,
                an exception will be thrown.
                For example: vlist = "x y", Function: "x+y", arg1={1,2}, arg2={5,10}, out={6,12} **/
            LIB_UTILITIES_EXPORT std::vector<double> Evaluate(std::vector<double> const* start, ...);

        private:

	        /** This is the class that is used as the grammar parser for the spirit engine. **/
	        class MathExpression : public boost_spirit::grammar<MathExpression>
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
		        static const int constantID		= 1;
		        static const int numberID		= 2;
		        static const int variableID		= 3;
		        static const int parameterID	= 4;
		        static const int functionID		= 5;
		        static const int factorID		= 6;
		        static const int operatorID		= 7;

		        MathExpression(boost_spirit::symbols<double> const* constants, std::vector<std::string> const& variables) :
					        boost_spirit::grammar<MathExpression>(), constants_p(constants), variables_p(variables) {}

		        template <typename ScannerT>
		        struct definition
		        {
			        /** This function specifies the grammar of the MathExpression parser. **/
			        definition(MathExpression const& self);

			        /** This holds the double value that is parsed by spirit so it can be stored in the AST. **/
			        double ParsedDouble;

			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<constantID> >		constant;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<numberID> >		number;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<variableID> >		variable;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<parameterID> >		parameter;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<functionID> >		function;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<factorID> >		factor;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		exponential;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		mult_div;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		add_sub;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		lt_gt;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		equality;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		logical_and;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		logical_or;
			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		expression;

			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> >		op;

			        boost_spirit::rule<ScannerT, boost_spirit::parser_context<>, boost_spirit::parser_tag<operatorID> > const&
				        start() const { return expression; }
		        };
	        };

	        /** This structure stores the data that was parsed in a tree. It uses
				TermType a lot to determine what kind of node it is so the appropriate
				action can be performed. **/
	        struct Node
	        {
		        enum TermType { DOUBLE, VARIABLE, PARAMETER, CONSTANT, FUNCTION,
						        ADDITION, SUBTRACTION, MULTIPLICATION, DIVISION,
						        EXPONENTIATION, LT, GT, LTEQ, GTEQ, EQUAL,
						        LOGICAL_AND, LOGICAL_OR };
		        enum AST_Mode { DEFAULT, SAVE_PARAMETERS, USE_SAVED_PARAMETERS };

		        union
		        {
			        double (*Function1)(double);
			        double (*Function2)(double, double);
			        double (*Function3)(double, double, double);
			        double (*Function4)(double, double, double, double);
		        };
		        double DoubleValue;
		        std::string* StringValue;
		        TermType DataType;
		        std::vector<Node*>* children;
		        double* PointerToVariable;	// pointer to var in map for quicker access

		        Node()
		        {
			        Function1 = NULL;
			        DoubleValue = 0;
			        StringValue = NULL;
			        DataType = DOUBLE;
			        children = NULL;
			        PointerToVariable = NULL;
		        }
		        ~Node()
		        {
			        if (StringValue != NULL) delete StringValue;
			        if (children != NULL)
			        {
				        for (std::vector<Node*>::iterator it = children->begin(); it != children->end(); it++)
					        delete *it;
				        delete children;
			        }
		        }
	        };

            /** This class stores pointers to the value for each map variable for the
                variables defined in #DefineFunction. Therefore, instead of having
                to always search the map for the variable name, you just call GetNext()
                since it is always in the same order (the order in which you specified
                the variables in #DefineFunction. This makes evaluating the
                function MUCH faster. **/
            class NextMapVariable
            {
            private:
                std::vector<double*> addresses;
                std::vector<double*>::size_type index;

            public:
                NextMapVariable(std::vector<double*> const& addr) : addresses(addr), index(0) { }
                ~NextMapVariable() { }
                double* GetNext() { return addresses[index++ % addresses.size()]; }
                void ResetIndex() { index = 0; }
            };

            /** This structure holds all of the data needed to evaluate an expression.
                Therefore, if another function is defined, this structure can be saved
                and reverted to later to parse data from a previous function definition. **/
            struct ParsedAST
            {
                /** This is a map that looks like <var_name, var_value>. It is set when
                    #Evaluate is called and used when evaluating the expression in
                    eval_expression. **/
                std::map<std::string, double>* VariableMap;

                /** This is a vector that holds the names of the variables specified in
                    the first parameter of #DefineFunction (vlist) in the order in
                    which they appear. It is used to fill in the VariableMap appropriately
                    so the variables are in the correct order. **/
                std::vector<std::string>* VariableVector;

                /** This is an object that stores the pointer to the VariableMap value.
                    Therefore, instead of having to always search the map for the variable
                    name, you just call NextMapVar->GetNext() since it is always in the
                    same order (the order in which you specified the variables in
                    #DefineFunction. This makes evaluating the function MUCH faster. **/
                NextMapVariable* NextMapVar;

                /** This is the number of variables that were specified in the vlist parameter
                    of #DefineFunction. This is the same number as VariableVector->size(). **/
                std::string::size_type NumberVariables;

                /** This stores the simplified AST that is created with #CreateAST in the
                    #DefineFunction function. **/
                Node* AST;

                /** This is the mode that will be used when evaluating the AST as it pertains to
                    the lookup of defined parameters. There are three options: DEFAULT,
                    SAVE_PARAMETERS, and USE_SAVED_PARAMETERS. These determine if the parameter
                    lookups should be cached in the node, so if the parameters don't change they
                    don't need to be looked up in the map each time. **/
                Node::AST_Mode ASTMode;

                /** This holds the string that was used for the function definition and variables
                    list. They are used to compare against the current state in #DefineFunction
                    to determine if a change is needed. **/
                std::string FunctionString;
                std::string VariableListString;

                ParsedAST(std::string functionStr, std::string varListStr)
                    : FunctionString(functionStr), VariableListString(varListStr)
                {
                    VariableMap = new std::map<std::string, double>;
                    VariableVector = new std::vector<std::string>;
                    NextMapVar = NULL;
                    NumberVariables = 0;
                    AST = NULL;
                    ASTMode = Node::SAVE_PARAMETERS;
                }

                ~ParsedAST()
                {
                    if (VariableMap != NULL) delete VariableMap;
                    if (VariableVector != NULL) delete VariableVector;
                    if (NextMapVar != NULL) delete NextMapVar;
                    if (AST != NULL) delete AST;
                }
            };

            /** This is the currently active ParsedAST that is being used. It can be
                changed with defining a new function in #DefineFunction. **/
            ParsedAST* ParsedData;

            /** These two containers hold processed string expressions defined
                via the second argument of #DefineFunction. Each processed expression
                gets its unique ID (to be stored in ParsedMapExprToId
                and unique AST (to be stored in ParsedASTs). ParsedMapExprToId is here for
                quick lookups of string expressions while ParsedASTs is split from the map
                in order to get constant time access to the AST via its ID. **/
            std::map<std::string, int> ParsedMapExprToId;
            std::vector<ParsedAST*> ParsedASTs;

	        /** This is a map that looks like <parameter_name, parameter_value>. It is set in
				the #SetParameters function and used for parameter lookup during function
				evaluation. Therefore, using this instead of constants (which are evaluated
				with #DefineFunction) is slower. **/
	        std::map<std::string, double>* ParametersMap;

	        /** This is a parser for spirit that parses the CONSTANT values. The default
				constants are those that are in math.h without the M_ prefix and they are
				initialized in the ExpressionEvaluator constructor. **/
	        boost_spirit::symbols<double>* constants_p;

			/** This function evaluates the AST created from #CreateAST and returns
				the result as a double. It will throw an exception if it encounters a
				parameter that isn't defined. If the "mode" argument is DEFAULT, it will
				look up the parameters from the ParametersMap. If "mode" is SAVE_PARAMETERS,
				it will get the current value from the map and then save it in the DoubleValue
				field. This can then be used with the USE_SAVED_PARAMETERS mode which will use
				the DoubleValue field for the parameter instead of the map. This can be used to
				speed up the vector #Evaluate function some since the parameters don't
				change there. **/
	        double EvaluateExpression(Node* const n, Node::AST_Mode mode);

	        /** This function walks the AST that is created with the spirit parser and creates a
				simplified AST that only holds the information required to parse the expression.
				It also performs any simplifications possible without having the final variable
				and parameter values. For example, it will simplifiy sin(5)*10+y to 9.04...+y so
				the calculation doesn't need to be done for every evaluation. It also performs the
				checks to make sure everything is in the correct range so these don't need to be
				performed at evaluation either. **/
	        Node* CreateAST(boost_spirit::tree_match<std::string::const_iterator,
				            boost_spirit::node_val_data_factory<double> >::tree_iterator const &i);

        };
    };
};
#endif
