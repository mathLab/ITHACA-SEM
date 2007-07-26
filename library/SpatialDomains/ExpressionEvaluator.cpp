// Expression Evaluator
//
// File:		ExpressionEvaluator.cpp
// Description:	Parses char* math expressions with Spirit.
// Author:		Michael DeLisi
// Date:		July 25, 2007

#include "ExpressionEvaluator.h"

using namespace std;
using namespace boost::spirit;

// ---------------------------------------------------------------------
//                       Math Function Declarations
// ---------------------------------------------------------------------

// This function is called each time a math function is parsed and executed. It checks
// the errno to see if there was an error that occured in the math.h function.
static void CheckMathOperationForErrors(string const& functionName)
{
	if (errno == EDOM)
		throw "Argument of " + functionName + " exceeds the range of the function.";
	else if (errno == ERANGE)
		throw "The result from " + functionName + " overflowed the double type.";

	errno = 0;
}

// Bessel function of the first kind. This is used since the math.h function
// requires the first argument to be an integer.
static double Jn (double i, double x)
{
	return jn((int) i, x);
}

// Bessel function of the second kind. This is used since the math.h function
// requires the first argument to be an integer.
static double Yn (double i, double x)
{
	return yn((int) i, x);
}

// ---------------------------------------------------------------------
//         Structures Used to Create Customized Spirit Parsers
// ---------------------------------------------------------------------

// This struct is a parser for spirit that parses the CONSTANT values.
// The constants are those that are in math.h without the M_ prefix.
static struct constants : symbols<double>
{
    constants()
    {
		add
			// Constants from math.h
			("E",			2.71828182845904523536)		// Natural logarithm 
			("LOG2E",		1.4426950408889634074)		// log_2 e 
			("LOG10E",		0.43429448190325182765)		// log_10 e 
			("LN2",			0.69314718055994530942)		// log_e 2 
			("LN10",		2.30258509299404568402)		// log_e 10 
			("PI",			3.14159265358979323846)		// pi 
			("PI_2",		1.57079632679489661923)		// pi/2 
			("PI_4",		0.78539816339744830962)		// pi/4 
			("1_PI",		0.31830988618379067154)		// 1/pi 
			("2_PI",		0.63661977236758134308)		// 2/pi 
			("2_SQRTPI",	1.12837916709551257390)		// 2/sqrt(pi) 
			("SQRT2",		1.41421356237309504880)		// sqrt(2) 
			("SQRT1_2",		0.70710678118654752440)		// 1/sqrt(2) 

			// Additional constants from the old parser
			("GAMMA",		0.57721566490153286060)		// Euler 
			("DEG",			57.29577951308232087680)	// deg/radian 
			("PHI",			1.61803398874989484820)		// golden ratio 
		;
    }
} constants_p;

// This is a helper function that stores a pointer to the math.h functions. A
// pointer to a function that returns a double cannot be used because it results
// in multiple overloaded function definitions matching. Therefore, I explicitly
// specify the arguments and just a union to store the correct version.
typedef double (*PFD)();
typedef double (*PFD1)(double);
typedef double (*PFD2)(double, double);
typedef double (*PFD3)(double, double, double);
typedef double (*PFD4)(double, double, double, double);
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

// This struct creates a parser that matches the function definitions from
// math.h. All of the functions accept one of more doubles as arguments and
// returns a double.
 
static struct functions : symbols<func>
{
	functions()
	{
		// Add all of the functions from math.h 
		add
			("abs",		abs)
			("asin",	asin)
			("atan",	atan)
			("atan2",	atan2)
			("ceil",	ceil)
			("cos",		cos)
			("cosh",	cosh)
			("exp",		exp)
			("fabs",	fabs)	// same as abs
			("floor",	floor)
			("fmod",	fmod)
			("log",		log)
			("log10",	log10)
			("pow",		pow)
			("sin",		sin)
			("sinh",	sinh)
			("sqrt",	sqrt)
			("tan",		tan)
			("tanh",	tanh)
			("hypot",	hypot)
			("j0",		j0)
			("j1",		j1)
			("jn",		Jn)
			("y0",		y0)
			("y1",		y1)
			("yn",		Yn)
		
			// These won't work because they require pointer or integer input. 
			//("frexp",	frexp)
			//("ldexp",	ldexp)
			//("modf",	modf)

			// These say they are defined, but they aren't in Visual Studio. 
			//("erf",		erf)
			//("erfc",	erfc)
			//("gamma",	gamma)
			//("lgamma",	lgamma)
			//("isnan",	isnan)
			//("acosh",	acosh)
			//("asinh",	asinh)
			//("atanh",	atanh)
			//("cbrt",	cbrt)
			//("expm1",	expm1)
			//("ilogb",	ilogb)
			//("log1p",	log1p)
			//("logb",	logb)
			//("nextafter",	nextafter)
			//("remainder",	remainder)
			//("rint",	rint)
			//("scalb",	scalb)

			// These are custom functions that were in manager.c that aren't in math.h. 
			//("rand",	Rand)		// random number (input the magnitude)
			//("bump",	Bump)
			//("Single",	Single)
			//("rad",		Radius)		// rad = sqrt(x^2 + y^2)
			//("ang",		Angle)		// ang = atan2(x,y)
			//("step",	Step)		// step(x,a) = 0 (if x < a) else 1
			//("step2",	Step2)		// step2(x,a) = 0 (if x <= a) else 1
			//("shock",	Shock)		// shock(x,a,b) = a (if x < 0), (a+b)/2 (if x==0) or b (if x > 0)
			//("jacobi",	Jacobi)
		;
	}
} functions_p;

// ---------------------------------------------------------------------
//             MathExpression definitions for Spirit Parser
// --------------------------------------------------------------------- 

// This function specifies the grammar of the MathExpression parser.
template <typename ScannerT>
ExpressionEvaluator::MathExpression::definition<ScannerT>::definition(MathExpression const& self)
{
	expression	=	multdiv >>
					*(  (root_node_d[ch_p('+')] >> multdiv)
					  | (root_node_d[ch_p('-')] >> multdiv)
					);

	multdiv		=	exponential >>
					*(  (root_node_d[ch_p('*')] >> exponential)
					  | (root_node_d[ch_p('/')] >> exponential)
					);

	exponential	=	factor >>
					*(	(root_node_d[ch_p('^')] >> factor)	);

	factor		=	constant
				|	number
				|	variable
				|	function
				|	parameter
				|	inner_node_d[ch_p('(') >> expression >> ch_p(')')]
				|	(root_node_d[ch_p('-')] >> factor);

	parameter	=	leaf_node_d[ lexeme_d[
						+alpha_p
					] ];

	function	=	root_node_d[functions_p] >> infix_node_d[inner_node_d[ch_p('(') >> expression >> *(',' >> expression) >> ch_p(')')]];

	variable	=	leaf_node_d[ lexeme_d[
						self.variables_p
					] ];

	number		=	leaf_node_d[ lexeme_d[
						real_p
					] ];

	constant	=	leaf_node_d[ lexeme_d[
						constants_p
					] ];
}

// ---------------------------------------------------------------------
//      ExpressionEvaluator Definitions that Wraps the Parser Class
// ---------------------------------------------------------------------

// Initializes the evaluator to a state where it is ready to accept input with
// the DefineFunction(...) function.
ExpressionEvaluator::ExpressionEvaluator(void)
{
	ParsedMap = new map<string, ParsedAST*>;
	ParsedData = NULL;
}

// Destroys the evaluator object and deletes all of the saves ASTs.
ExpressionEvaluator::~ExpressionEvaluator(void)
{
	if (ParsedMap != NULL)
	{
		for (map<string, ParsedAST*>::iterator it = ParsedMap->begin(); it != ParsedMap->end(); it++)
			if (it->second != NULL)
				delete it->second;
		delete ParsedMap;
	}
}


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
void ExpressionEvaluator::DefineFunction(char const *vlist, char const *function)
{
	// Find the previous parsing, or create a new one if it doesn't exist.
	string functionStr(function);
	map<string, ParsedAST*>::iterator it = ParsedMap->find(functionStr);
	if (it == ParsedMap->end())
	{
		ParsedData = new ParsedAST;
	}
	else
	{
		ParsedData = it->second;
		return;
	}

	// Parse the vlist input and separate the variables into ordered entries
	// in a vector<string> object. These need to be ordered because this is
	// the order the variables will get assigned to in the Map when Evaluate(...)
	// is called.
	ParsedData->VariableVector->clear();
	parse(vlist, ( *space_p >>
				   *(
						+(+graph_p)[push_back_a(*(ParsedData->VariableVector))]
							>> +space_p
					 )
				 )
		 );
	ParsedData->NumberVariables = ParsedData->VariableVector->size();
	for (vector<string>::iterator it = ParsedData->VariableVector->begin();
			it != ParsedData->VariableVector->end(); it++)
		(*ParsedData->VariableMap)[*it] = 0;

	// Do the actual parsing with spirit and alert the user if there was an error with an exception.
	MathExpression myGrammar(ParsedData->VariableVector);
	tree_parse_info<> parseInfo = ast_parse(function, myGrammar, space_p);
	if (parseInfo.full == false)
	{
		delete ParsedData;
		throw string("Unable to fully parse function. Stopped just before: ") + string(parseInfo.stop).substr(0, 15);
		return;
	}

	// This creates a simplified AST that only holds the information required to
	// parse the expression. It also performs any simplifications possible without
	// having the final variable and parameter values.
	ParsedData->AST = CreateAST(parseInfo.trees.begin());

	// Store the AST so this function can be Evaluated again very quickly.
	(*ParsedMap)[functionStr] = ParsedData;
}

// This function accepts parameters in the form of <parameter_name, double_value>.
// It is called after DefineFunction(...) and before Evaluate(...). This function
// will replace all of the parameters with whatever is in the input stl map. It will
// not add to what was originally set.
void ExpressionEvaluator::SetParameters(map<string, double> const& params)
{
	if (ParsedData == NULL)
	{
		throw string("Unable to set parameters because a function must first be defined with DefineFunction(...).");
		return;
	}

	if (ParsedData->ParametersMap == NULL)
		ParsedData->ParametersMap = new map<string, double>;
	else
		ParsedData->ParametersMap->clear();

	for (map<string, double>::const_iterator it = params.begin(); it != params.end(); it++)
		(*ParsedData->ParametersMap)[it->first] = it->second;
}

// This function evaluates the function defined from DefineFunction(...).
// The arguments to the function are the values that were defined in the
// first argument of DefineFunction(...) in the same order as they were
// listed. It results the result of the evaluation as a double. If a function
// is not currently defined, behavior may be unpredictable since I do not
// know how many arguments were passed in, so it will probably not be cleaned
// up correctly.
double ExpressionEvaluator::Evaluate(double start, ...)
{
	va_list ap;
	
	va_start(ap, start);
	if (ParsedData == NULL)
	{
		va_end(ap);
		throw string("Unable to evaluate because a function must first be defined with DefineFunction(...).");
		return -1;
	}
	else
	{
		ParsedData->VariableMap->clear();
		if (ParsedData->NumberVariables > 0)
		{
			(*ParsedData->VariableMap)[(*ParsedData->VariableVector)[0]] = start;
			for (string::size_type i = 1 ; i < ParsedData->NumberVariables; i++)
				(*ParsedData->VariableMap)[(*ParsedData->VariableVector)[i]] = va_arg(ap, double);
		}
	}
	va_end(ap);

	return EvaluateExpression(ParsedData->AST);
}

// This function evaluates the function defined from DefineFunction(...).
// This function accepts one vector<double> argument that has the values
// for the variables that were listed in the first argument (vlist) of the
// DefineFunction(...) function. They must be in the vector in the same
// order as they appeared in the vlist.
double ExpressionEvaluator::Evaluate(vector<double> const& vars)
{
	if (ParsedData == NULL)
	{
		throw string("Unable to evaluate because a function must first be defined with DefineFunction(...).");
		return -1;
	}

	if (ParsedData->NumberVariables != vars.size())
	{
		throw string("Illegal number of variables in the Evaluate input vector.");
		return -1;
	}

	ParsedData->VariableMap->clear();
	if (ParsedData->NumberVariables > 0)
		for (string::size_type i = 0 ; i < ParsedData->NumberVariables; i++)
			(*ParsedData->VariableMap)[(*ParsedData->VariableVector)[i]] = vars[i];

	return EvaluateExpression(ParsedData->AST);
}

// This function walks the AST that is created with the spirit parser and creates a
// simplified AST that only holds the information required to parse the expression.
// It also performs any simplifications possible without having the final variable
// and parameter values. For example, it will simplifiy sin(5)*10+y to 9.04...+y so
// the calculation doesn't need to be done for every evaluation. It also performs the
// checks to make sure everything is in the correct range so these don't need to be
// performed at evaluation either.
ExpressionEvaluator::Node* ExpressionEvaluator::CreateAST(tree_match<const char*>::tree_iterator const& i)
{
	const parser_id parserID = i->value.id();
	if (parserID == MathExpression::constantID)
	{
		string constantName(i->value.begin(), i->value.end());

		if (i->children.size() != 0)
		{
			throw "Illegal children under constant node: " + constantName;
			return NULL;
		}

		double* value = find(constants_p, constantName.c_str() );
		if (value == NULL)
		{
			throw "Cannot find the value for the specified constant: " + constantName;
			return NULL;
		}

		Node* n = new Node;
		n->DataType = Node::CONSTANT;
		n->DoubleValue = *value;
		return n;
	}
	else if (parserID == MathExpression::numberID)
	{
		string doubleString(i->value.begin(), i->value.end());
		if (i->children.size() != 0)
		{
			throw "Illegal children under number node: " + doubleString;
			return NULL;
		}

		Node* n = new Node;
		n->DataType = Node::DOUBLE;

		// Parse the number and check for over/underflows.
		errno = 0;
		n->DoubleValue = strtod(doubleString.c_str(), NULL);
		if (errno == ERANGE)
		{
			errno = 0;
			if (n->DoubleValue == HUGE_VAL)
			{
				delete n;
				throw "Double value in function resulted in overflow: " + doubleString;
				return NULL;
			}
			else if (n->DoubleValue == 0)
			{
				delete n;
				throw "Double value in function resulted in underflow: " + doubleString;
				return NULL;
			}
		}

		return n;
	}
	else if (parserID == MathExpression::variableID)
	{
		if (i->children.size() != 0)
		{
			throw "Illegal children under variable node: " + string(i->value.begin(), i->value.end());
			return NULL;
		}

		Node* n = new Node;
		n->DataType = Node::VARIABLE;
		n->StringValue = new string(i->value.begin(), i->value.end());
		return n;
	}
	else if (parserID == MathExpression::parameterID)
	{
		if (i->children.size() != 0)
		{
			throw "Illegal children under parameter node: " + string(i->value.begin(), i->value.end());
			return NULL;
		}

		Node* n = new Node;
		n->DataType = Node::PARAMETER;
		n->StringValue = new string(i->value.begin(), i->value.end());
		return n;
	}
	else if (parserID == MathExpression::functionID)
	{
		string fname(i->value.begin(), i->value.end());
		func* funcptr = find(functions_p, fname.c_str());
		if (funcptr == NULL)
		{
			throw "Invalid function specified: " + fname;
			return NULL;
		}
		if (funcptr->size != i->children.size())
		{
			throw "Illegal number or arguments for math function: " + fname;
			return NULL;
		}

		bool allDoubles = true;
		vector<Node*>* arguments = new vector<Node*>;
		for (tree_match<char const*>::tree_iterator it = i->children.begin(); it != i->children.end(); it++)
		{
			Node* node = CreateAST(it);
			if (node->DataType != Node::DOUBLE)
				allDoubles = false;
			arguments->push_back(node);
		}

		// If allDoubles == true, we can evaluate the expression since there are no unknowns.
		if (allDoubles == true)
		{
			Node* node = new Node;
			node->DataType = Node::DOUBLE;
			errno = 0;

			switch (i->children.size())
			{
			case 1:
				node->DoubleValue = (*funcptr->func1)( (*arguments)[0]->DoubleValue );
				break;
			case 2:
				node->DoubleValue = (*funcptr->func2)( (*arguments)[0]->DoubleValue, (*arguments)[1]->DoubleValue );
				break;
			case 3:
				node->DoubleValue = (*funcptr->func3)( (*arguments)[0]->DoubleValue, (*arguments)[1]->DoubleValue, (*arguments)[2]->DoubleValue );
				break;
			case 4:
				node->DoubleValue = (*funcptr->func4)( (*arguments)[0]->DoubleValue, (*arguments)[1]->DoubleValue, (*arguments)[2]->DoubleValue, (*arguments)[3]->DoubleValue );
				break;
			}

			// Delete the arguments vector.
			for (vector<Node*>::iterator it = arguments->begin(); it != arguments->end(); it++)
				delete *it;
			delete arguments;

			CheckMathOperationForErrors(fname);

			return node;
		}
		else
		{
			Node* node = new Node;
			node->DataType = Node::FUNCTION;
			node->StringValue = new string(fname);
			switch (i->children.size())
			{
			case 1:
				node->Function1 = funcptr->func1;
				break;
			case 2:
				node->Function2 = funcptr->func2;
				break;
			case 3:
				node->Function3 = funcptr->func3;
				break;
			case 4:
				node->Function4 = funcptr->func4;
				break;
			}
			node->children = arguments;

			return node;
		}
	}
	else if (parserID == MathExpression::factorID)
	{
		if (*i->value.begin() != '-')
		{
			throw "Illegal factor - it can only be '-' and it was: " + string(i->value.begin(), i->value.end());
			return NULL;
		}
		if (i->children.size() != 1)
		{
			throw "Illegal number of children under factor node: " + string(i->value.begin(), i->value.end());
			return NULL;
		}

		Node* value = CreateAST(i->children.begin());
		if (value->DataType == Node::DOUBLE)
		{
			value->DoubleValue = -1 * value->DoubleValue;
			return value;
		}
		else
		{
			Node* zero = new Node;
			zero->DataType = Node::DOUBLE;
			zero->DoubleValue = 0.;

			Node* subtraction = new Node;
			subtraction->DataType = Node::SUBTRACTION;
			subtraction->children = new vector<Node*>;
			subtraction->children->push_back(zero);
			subtraction->children->push_back(value);

			return subtraction;
		}
	}
	else if (parserID == MathExpression::expressionID)
	{
		if (i->children.size() != 2)
		{
			throw "Too many arguments for mathematical operator: " + string(i->value.begin(), i->value.end());
			return NULL;
		}

		Node* left = CreateAST(i->children.begin()+0);
		Node* right = CreateAST(i->children.begin()+1);

		// If they are both of double type, we can simplify to a number.
		if (left->DataType == Node::DOUBLE && right->DataType == Node::DOUBLE)
		{
			Node* node = new Node;
			node->DataType = Node::DOUBLE;
			switch(*i->value.begin())
			{
			case '+':
				node->DoubleValue = left->DoubleValue + right->DoubleValue;
				break;
			case '-':
				node->DoubleValue = left->DoubleValue - right->DoubleValue;
				break;
			case '*':
				node->DoubleValue = left->DoubleValue * right->DoubleValue;
				break;
			case '/':
				node->DoubleValue = left->DoubleValue / right->DoubleValue;
				break;
			case '^':
				node->DoubleValue = pow(left->DoubleValue, right->DoubleValue);
				break;
			default:
				delete left;
				delete right;
				delete node;
				throw "Invalid operator encountered: " + string(i->value.begin(), i->value.end());
				return NULL;
			}

			delete left, right;
			return node;
		}
		else
		{
			Node* node = new Node;
			switch(*i->value.begin())
			{
			case '+':
				node->DataType = Node::ADDITION;
				break;
			case '-':
				node->DataType = Node::SUBTRACTION;
				break;
			case '*':
				node->DataType = Node::MULTIPLICATION;
				break;
			case '/':
				node->DataType = Node::DIVISION;
				break;
			case '^':
				node->DataType = Node::EXPONENTIATION;
				break;
			default:
				delete left;
				delete right;
				delete node;
				throw "Invalid operator encountered: " + string(i->value.begin(), i->value.end());
				return NULL;
			}

			node->children = new vector<Node*>;
			node->children->push_back(left);
			node->children->push_back(right);
			return node;
		}
	}

	throw "Illegal expression encountered: " + string(i->value.begin(), i->value.end());
	return NULL;
}


// This function evaluates the AST created from CreateAST(...) and returns
// the result as a double. It will throw an exception if it encounters a
// parameter that isn't defined.
double ExpressionEvaluator::EvaluateExpression(Node* const n)
{
	switch (n->DataType)
	{
	case Node::CONSTANT:
	case Node::DOUBLE:
		return n->DoubleValue;
	case Node::VARIABLE:
		return (*ParsedData->VariableMap)[*n->StringValue];
	case Node::PARAMETER:
		{
			map<string, double>::iterator it;
			if (ParsedData->ParametersMap == NULL ||
					(it = ParsedData->ParametersMap->find(*n->StringValue)) == ParsedData->ParametersMap->end())
			{
				throw "Illegal parameter specified: " + *n->StringValue;
				return -1;
			}
			return it->second;
		}
	case Node::FUNCTION:
		{
			errno = 0;
			double result = 0;
			switch (n->children->size())
			{
			case 1:
				result = (*n->Function1)( EvaluateExpression((*n->children)[0]) );
				break;
			case 2:
				result = (*n->Function2)( EvaluateExpression((*n->children)[0]), EvaluateExpression((*n->children)[1]) );
				break;
			case 3:
				result = (*n->Function3)( EvaluateExpression((*n->children)[0]), EvaluateExpression((*n->children)[1]), EvaluateExpression((*n->children)[2]) );
				break;
			case 4:
				result = (*n->Function4)( EvaluateExpression((*n->children)[0]), EvaluateExpression((*n->children)[1]), EvaluateExpression((*n->children)[2]), EvaluateExpression((*n->children)[3]) );
				break;
			}
			CheckMathOperationForErrors(*n->StringValue);
			return result;
		}
	case Node::ADDITION:
		return EvaluateExpression((*n->children)[0]) + EvaluateExpression((*n->children)[1]);
	case Node::SUBTRACTION:
		return EvaluateExpression((*n->children)[0]) - EvaluateExpression((*n->children)[1]);
	case Node::MULTIPLICATION:
		return EvaluateExpression((*n->children)[0]) * EvaluateExpression((*n->children)[1]);
	case Node::DIVISION:
		return EvaluateExpression((*n->children)[0]) / EvaluateExpression((*n->children)[1]);
	case Node::EXPONENTIATION:
		return pow(EvaluateExpression((*n->children)[0]), EvaluateExpression((*n->children)[1]));
	}

	throw string("Illegal expression was encountered.");
	return -1;
}

#if 0

// ---------------------------------------------------------------------
//             Math Functions that are Used when Parsing Input
// ---------------------------------------------------------------------

static double Pow(double x, double y)
{
	const double yn = floor(y + .5);
	double px = 1.;

	if (yn >= 0 && yn == y)	// Do it inline if y is an integer power 
	{
		register int n = (int) yn;
		while (n--) 
			px *= x;
	}
	else  
		px = pow(x,y);

	return px;
}

static double Rand(double x)
{
	return x * drand();
}

static double Bump(double x)
{
	if(x >= 0. && x < .125)
		return -1;
	if(x >= 0.125 && x < .25)
		return 0.;
	if(x >= 0.25 && x < .375)
		return 1.;
	if(x >= 0.375 && x <= .5)
		return 0.;

	return -9999.;
}

static double Radius(double x, double y)
{
	if (x != 0. || y != 0.)
		return sqrt (x*x + y*y);
	else
		return 0.;
}

static double Angle(double x, double y)
{
	double theta = 0.;

	if ((x != 0.) || (y != 0.))
		theta = atan2 (y,x);

	return theta;
}

// Heaviside step function H(x-a) =1 if x >= 0 else =0 
static double Step(double x, double a)
{
	double H = 1.0;
	if (x < a)
		H = 0.0;

	return H;
}

// Heaviside step function H(x-a) =1 if x >= 0 else =0 
static double Step2(double x, double a)
{
	double H = 1.0;
	if (x <= a)
		H = 0.0;

	return H;
}

static double Shock(double x, double a, double b)
{
	if(x == 0)
		return 0.5*(a+b);
	if(x > 0)
		return b;
	if(x < 0)
		return a;
	return 0;
}

// -----------------------------------------------------------------
//   jacobi() - jacobi polynomials 
//   
//   Get a vector 'poly' of values of the n_th order Jacobi polynomial
//   P^(alpha,beta)_n(z) alpha > -1, beta > -1 at the z
// ----------------------------------------------------------------- 
static double Jacobi(double z, double n, double alpha, double beta)
{
	register int k;
	double one = 1.0;
	double a1,a2,a3,a4;
	double two = 2.0, apb = alpha + beta;
	double poly, polyn1,polyn2;

	polyn2 = one;
	polyn1 = 0.5*(alpha - beta + (alpha + beta + 2)*z);

	for(k = 2; k <= n; ++k)
	{
		a1 =  two*k*(k + apb)*(two*k + apb - two);
		a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
		a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
		a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);

		a2 /= a1;
		a3 /= a1;
		a4 /= a1;

		poly   = (a2 + a3*z)*polyn1 - a4*polyn2;
		polyn2 = polyn1;
		polyn1 = poly;
	}

	return poly;
}

static double Single(double x, double y)
{
#if 1
	double gamma = 64.0*64.0;
	double tmp;

	if (y>=3.0 && y<=4.0)
	{
		tmp = (y-3.)*(y-4.)*(y-3.)*(y-4.);
		if(x>=1.0 && x<=2.)
			return gamma*(x-1.0)*(x-2.)*(x-1.0)*(x-2.);
		if(x>=3.0 && x<=4.0)
			return gamma*(x-3.0)*(x-4.0)*(x-3.0)*(x-4.0);
		if(x>=5. && x<=6.) 
			return gamma*(x-5.0)*(x-6.)*(x-5.0)*(x-6.);
	}
#else
	double gamma = 64.0*64.0;
	double xa,xb,ya,yb;

	if(x>1. && x<2. && y>1. && y<2.)
	{
		xa = 1.;    xb = 2.;
		ya = 1.;    yb = 2.;
	}
	else if(x>3. && x<4. && y>3. && y<4.)
	{
		xa = 3.;    xb = 4.;
		ya = 3.;    yb = 4.;
	}
	else if(x>5. && x<6. && y>5. && y<6.)
	{
		xa = 5.;    xb = 6.;
		ya = 5.;    yb = 6.;
	}
	else
		return 0.;

	return gamma*(x-xa)*(x-xa)*(x-xb)*(x-xb)*(y-ya)*(y-ya)*(y-yb)*(y-yb);
#endif
}

#endif
