// Expression Evaluator
//
// File:		ExpressionEvaluator.hpp
// Description:	Implementation of an expression parser.
// Author:		Michael DeLisi
// Date:		July 25, 2007

#ifndef _EXPRESSION_EVALUATOR_HPP
#define _EXPRESSION_EVALUATOR_HPP

#include <boost/spirit/core.hpp>
#include <boost/spirit/tree/ast.hpp>
#include <boost/spirit/symbols/symbols.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>

#include <iostream>
#include <stdarg.h>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <errno.h>
#include <stdlib.h>

class ExpressionEvaluator
{
public:
	// Initializes the evaluator to a state where it is ready to accept input with
	// the DefineFunction(...) function.
	ExpressionEvaluator(void);

	 // Destroys the evaluator object and deletes all of the saves ASTs.
	~ExpressionEvaluator(void);

	// This function allows one to define a function to evaluate. The first argument (vlist)
	// is a list of variables (separated by spaces) that the second argument (function)
	// depends on. For example, if function = "x + y", then vlist should most likely be
	// "x y", unless you are defining x or y as parameters with SetParameters(...).
	//
	// If the function definition is changed, all of the current values such as parsed AST and
	// parameters are saved if you call DefineFunction with the same "function" string. Also
	// note that if you have previously set a function to use a certain "vlist", and then
	// you call DefineFunction again with the same function string but a different vlist string,
	// it will use the old vlist string from the first declaration.
	void DefineFunction(char const *vlist, char const *function);

	// This function accepts parameters in the form of <parameter_name, double_value>.
	// It is called after DefineFunction(...) and before Evaluate(...). This function
	// will replace all of the parameters with whatever is in the input stl map. It will
	// not add to what was originally set.
	void SetParameters(std::map<std::string, double> const& params);

	 
	// This function evaluates the function defined from DefineFunction(...).
	// The arguments to the function are the values that were defined in the
	// first argument of DefineFunction(...) in the same order as they were
	// listed. It results the result of the evaluation as a double. If a function
	// is not currently defined, behavior may be unpredictable since I do not
	// know how many arguments were passed in, so it will probably not be cleaned
	// up correctly.
	double Evaluate(double start = 0, ...);

	 
	// This function evaluates the function defined from DefineFunction(...).
	// This function accepts one vector<double> argument that has the values
	// for the variables that were listed in the first argument (vlist) of the
	// DefineFunction(...) function. They must be in the vector in the same
	// order as they appeared in the vlist.
	double Evaluate(std::vector<double> const& vars);

private:
	// This is the class that is used as the grammar parser for the spirit engine.
	class MathExpression : public boost::spirit::grammar<MathExpression>
	{
	private:
		// Variables is a customized parser that will match the variables that the function
		// depends on (the first argument of DefineFunction(...)).
		struct variables : boost::spirit::symbols<double*>
		{
			variables(std::vector<std::string> const* vars)
			{
				for (std::vector<std::string>::const_iterator it = vars->begin(); it != vars->end(); it++)
					add(it->c_str(), 0);
			}
		} variables_p;

	public:
		// These constants are used to determine what parser was used to parse what value,
		// which allows for type identification when analyzing the parsed AST.
		static const int constantID		= 1;
		static const int numberID		= 2;
		static const int variableID		= 3;
		static const int parameterID	= 4;
		static const int functionID		= 5;
		static const int factorID		= 6;
		static const int expressionID	= 7;

		MathExpression(std::vector<std::string> const* vars) : boost::spirit::grammar<MathExpression>(), variables_p(vars) {}

		template <typename ScannerT>
		struct definition
		{
			// This function specifies the grammar of the MathExpression parser.
			definition(MathExpression const& self);

			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<constantID> >		constant;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<numberID> >		number;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<variableID> >		variable;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<parameterID> >		parameter;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<functionID> >		function;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<factorID> >		factor;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<expressionID> >	exponential;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<expressionID> >	multdiv;
			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<expressionID> >	expression;

			boost::spirit::rule<ScannerT, boost::spirit::parser_context<>, boost::spirit::parser_tag<expressionID> > const&
				start() const { return expression; }
		};
	};

	// This structure stores the data that was parsed in a tree. It uses
	// TermType a lot to determine what kind of node it is so the appropriate
	// action can be performed.
	struct Node
	{
		enum TermType { DOUBLE, VARIABLE, PARAMETER, CONSTANT, FUNCTION, ADDITION, SUBTRACTION, MULTIPLICATION, DIVISION, EXPONENTIATION };

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

		Node()
		{
			Function1 = NULL;
			Function2 = NULL;
			Function3 = NULL;
			Function4 = NULL;
			DoubleValue = 0;
			StringValue = NULL;
			DataType = DOUBLE;
			children = NULL;
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

	
	// This structure holds all of the data needed to evaluate an expression.
	// Therefore, if another function is defined, this structure can be saved
	// and reverted to later to parse data from a previous function definition.
	struct ParsedAST
	{
		
		// This is a map that looks like <var_name, var_value>. It is set when
		// Evaluate(...) is called and used when evaluating the expression in
		// eval_expression(...).
		std::map<std::string, double>* VariableMap;

		
		// This is a vector that holds the names of the variables specified in
		// the first parameter of DefineFunction(...) (vlist) in the order in
		// which they appear. It is used to fill in the VariableMap appropriately
		// so the variables are in the correct order.
		std::vector<std::string>* VariableVector;

		
		// This is a map that looks like <parameter_name, parameter_value>. It is
		// set in the SetParameters(...) function and used when evaluating the
		// expression in eval_expression(...).
		std::map<std::string, double>* ParametersMap;

		// This is the number of variables that were specified in the vlist parameter
		// of DefineFunction(...). This is the same number as VariableVector->size().
		std::string::size_type NumberVariables;

		// This stores the simplified AST that is created with CreateAST(...) in the
		// DefineFunction(...) function.
		Node* AST;

		ParsedAST()
		{
			VariableMap = new std::map<std::string, double>;
			VariableVector = new std::vector<std::string>;
			ParametersMap = NULL;
			NumberVariables = 0;
			AST = NULL;
		}

		~ParsedAST()
		{
			if (VariableMap != NULL) delete VariableMap;
			if (VariableVector != NULL) delete VariableVector;
			if (ParametersMap != NULL) delete ParametersMap;
			if (AST != NULL) delete AST;
		}
	};

	// This is the currently active ParsedAST that is being used. It can be
	// changed with defining a new function in DefineFunction(...).
	ParsedAST* ParsedData;

	// This is a map of <string, ParsedAST*>. The string key is the second argument
	// of DefineFunction(...) that is used to find the ParsedAST* again. The found
	// AST is set to be default be storing it in ParsedData.
	std::map<std::string, ParsedAST*>* ParsedMap;

	// This function evaluates the AST created from CreateAST(...) and returns
	// the result as a double. It will throw an exception if it encounters a
	// parameter that isn't defined.
	double ExpressionEvaluator::EvaluateExpression(Node* const n);

	
	// This function walks the AST that is created with the spirit parser and creates a
	// simplified AST that only holds the information required to parse the expression.
	// It also performs any simplifications possible without having the final variable
	// and parameter values. For example, it will simplifiy sin(5)*10+y to 9.04...+y so
	// the calculation doesn't need to be done for every evaluation. It also performs the
	// checks to make sure everything is in the correct range so these don't need to be
	// performed at evaluation either.
	Node* ExpressionEvaluator::CreateAST(boost::spirit::tree_match<const char*>::tree_iterator const& n);

};

#endif
