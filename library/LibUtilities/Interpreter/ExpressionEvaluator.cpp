// Expression Evaluator
//
// File:		ExpressionEvaluator.cpp
// Description:	Parses string math expressions with Spirit.
// Author:		Michael DeLisi
// Date:		July 27, 2007
#include "ExpressionEvaluator.h"
#include <stdexcept>

#ifdef _MSC_VER
#include <boost/preprocessor/cat.hpp>  
#endif //MSC_VER

using namespace std;

#if( BOOST_VERSION / 100 % 1000 >= 36 )
using namespace boost::spirit::classic;
#else
using namespace boost::spirit;
#endif

#ifdef _MSC_VER
#define NEKTAR_MATH_NAME(x) BOOST_PP_CAT(_, x)
#else
#define NEKTAR_MATH_NAME(x) x
#endif

namespace Nektar
{
    namespace LibUtilities
    {

        // ---------------------------------------------------------------------
        //                       Math Function Declarations
        // ---------------------------------------------------------------------

        /** This function is called each time a math function is parsed and executed. It checks
			the errno to see if there was an error that occured in the math.h function. **/
        static void CheckMathOperationForErrors(string const& functionName)
        {
            if (errno == EDOM)
                throw std::runtime_error("Argument of " + functionName + " exceeds the range of the function.");
            else if (errno == ERANGE)
                throw std::runtime_error("The result from " + functionName + " overflowed the double type.");

            errno = 0;
        }

        /** Bessel function of the first kind. This is used since the math.h function
			requires the first argument to be an integer. **/
        static double Jn (double i, double x)
        {
            return NEKTAR_MATH_NAME(jn) ((int) i, x);
        }

        /** Bessel function of the second kind. This is used since the math.h function
			requires the first argument to be an integer. **/
        static double Yn (double i, double x)
        {
            return NEKTAR_MATH_NAME(yn) ((int) i, x);
        }
		
        // ---------------------------------------------------------------------
        //        Definition of the Customized Function Parser for Spirit
        // ---------------------------------------------------------------------

        /** This is a helper function that stores a pointer to the
            math.h functions. A pointer to a function that returns a
            double cannot be used because it results in multiple
            overloaded function definitions matching. Therefore, I
            explicitly specify the arguments and just a union to store
            the correct version. **/

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

        /** This struct creates a parser that matches the function
 	    definitions from math.h. All of the functions accept one
 	    of more doubles as arguments and returns a double. **/
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
                    ("hypot",	NEKTAR_MATH_NAME(hypot))
                    ("j0",      NEKTAR_MATH_NAME(j0))
                    ("j1",		NEKTAR_MATH_NAME(j1))
                    ("jn",		Jn)
                    ("y0",		NEKTAR_MATH_NAME(y0))
                    ("y1",		NEKTAR_MATH_NAME(y1))
                    ("yn",		Yn)

                    // These won't work because they require pointer or integer input. 
                    //("frexp",	frexp)
                    //("ldexp",	ldexp)
                    //("modf",	modf)

                    // These say they are defined in math.h, but they won't compile in Visual Studio. 
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
                    // All of these functions are commented out at the bottom of this file, so to add
                    // one just move it up to the top where with the other math functions and uncomment.
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

        /** This function specifies the grammar of the MathExpression parser. **/
        template <typename ScannerT>
        ExpressionEvaluator::MathExpression::definition<ScannerT>::definition(MathExpression const& self)
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
                //access_node_d[real_p[assign_a(ParsedDouble)]][StoreParsedDoubleInAST(ParsedDouble)]
				real_p
            ] ] >> op;

            constant	=	leaf_node_d[ lexeme_d[
                *self.constants_p
            ] ] >> op;

            op = eps_p( end_p | "||" | "&&" | "==" | "<=" | ">=" | '<' | '>' | '+' | '-' | '*' | '/' | '^' | ')' );
        }

        // ---------------------------------------------------------------------
        //      ExpressionEvaluator Definitions that Wraps the Parser Class
        // ---------------------------------------------------------------------

        // \brief Initializes the evaluator. Call DefineFunction(...) next.
        ExpressionEvaluator::ExpressionEvaluator(void)
        {
            //ParsedMap = new <ParsedAST*>;
            ParsedData = NULL;
            ParametersMap = new map<string, double>;

            // Populate constants with default values.
            constants_p = new symbols<double>;
            constants_p->add
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

        /** \brief Destroys the evaluator object and deletes all of the saves ASTs. **/
        ExpressionEvaluator::~ExpressionEvaluator(void)
        {
            if (!ParsedASTs.empty())
            {
                for (vector<ParsedAST*>::iterator it = ParsedASTs.begin(); it != ParsedASTs.end(); it++)
                        delete *it;
            }
            if (ParametersMap != NULL) delete ParametersMap;
            if (constants_p != NULL) delete constants_p;
        }


        /** \brief Initializes evaluator with list of non-constant/parameter variables and the
            function to evaluate. **/
        int ExpressionEvaluator::DefineFunction(const std::string& vlist, const std::string& function)
        {

            // Find the previous parsing.
            map<string, int>::const_iterator it = ParsedMapExprToId.find(function);
            if (it != ParsedMapExprToId.end())
            {
                // if this function is already defined, don't do anything but
                // return its ID.
                ParsedData = ParsedASTs[it->second];
                ParsedData->ASTMode = Node::SAVE_PARAMETERS;
                return it->second;
            }

            // Create a new parsing if a previous one doesn't exist.
            string functionStr(function.c_str());
            string vlistStr(vlist.c_str());
            ParsedData = new ParsedAST(functionStr, vlistStr);

            // Parse the vlist input and separate the variables into ordered entries
            // in a vector<string> object. These need to be ordered because this is
            // the order the variables will get assigned to in the Map when Evaluate(...)
            // is called.
            vector<string> VariableVector;
            parse((char*) vlistStr.c_str(), ( *space_p >>
                       *(
                                +(+graph_p)[push_back_a(VariableVector)]
                                    >> +space_p
                             )
                        )
             );
            ParsedData->NumberVariables = VariableVector.size();
            vector<double*> mapAddresses;
            for (vector<string>::iterator it = VariableVector.begin(); it != VariableVector.end(); it++)
            {
                (*ParsedData->VariableMap)[*it] = 0;
                mapAddresses.push_back(&(*ParsedData->VariableMap)[*it]);
            }
            ParsedData->NextMapVar = new NextMapVariable(mapAddresses);

            // Do the actual parsing with spirit and alert the user if there was an error with an exception.
            MathExpression myGrammar(constants_p, VariableVector);
            //const char* last = functionStr.c_str(); while (*last != '\0') last++;
            //tree_parse_info<const char*, node_val_data_factory<double> > parseInfo = ast_parse<node_val_data_factory<double> >((const char*) functionStr.c_str(), (const char*) last, myGrammar, space_p);
            tree_parse_info<string::const_iterator, node_val_data_factory<double> > parseInfo = ast_parse<node_val_data_factory<double>, string::const_iterator, MathExpression, space_parser>(functionStr.begin(), functionStr.end(), myGrammar, space_p);
            if (parseInfo.full == false)
            {
                delete ParsedData;
                throw std::runtime_error("Unable to fully parse function. Stopped just before: " + string(parseInfo.stop, parseInfo.stop + 15));
                return -1;
            }

            // This creates a simplified AST that only holds the information required to
            // parse the expression. It also performs any simplifications possible without
            // having the final variable and parameter values.
            ParsedData->AST = CreateAST(parseInfo.trees.begin());

            // Store the AST so this function can be Evaluated again very quickly.
            //(*ParsedMap)[functionStr] = ParsedData;
            ParsedASTs.push_back(ParsedData);
            ParsedMapExprToId[functionStr] = ( ParsedASTs.size() - 1 );

            return ParsedASTs.size() - 1;
        }


        /** \brief Adds constants that can be used in all functions defined in the future. **/
        void ExpressionEvaluator::AddConstants(std::map<std::string, double> const& constants)
        {
            vector<string> AlreadyAdded;
            for (map<string, double>::const_iterator it = constants.begin(); it != constants.end(); it++)
            {
                if (find(*constants_p, it->first.c_str()) != NULL)
                    AlreadyAdded.push_back(it->first);
                else
                    constants_p->add(it->first.c_str(), it->second);
            }

            if (AlreadyAdded.size() > 0)
            {
                string exception("The following constant(s) were not added because they already existed: ");
                for (vector<string>::iterator it = AlreadyAdded.begin(); it != AlreadyAdded.end(); it++)
                {
                    if (it == AlreadyAdded.begin())
                        exception += *it;
                    else
                        exception += ", " + *it;
                }

                throw std::runtime_error(exception + ".");
            }
        }

		/** \brief Adds a constant that can be used in all functions defined in the future. **/
		void ExpressionEvaluator::AddConstant(std::string const& name, double value)
        {
            if (find(*constants_p, name.c_str()))
                throw std::runtime_error("Cannot add specified constant because it is already added: " + name);
            else
                constants_p->add(name.c_str(), value);
        }

        /** \brief Returns the constant's value. **/
		double ExpressionEvaluator::GetConstant(std::string const& name)
        {
            double* value = find(*constants_p, name.c_str());
            if (value == NULL)
            {
                throw std::runtime_error("Constant variable not found: " + name);
                return -1;
            }

            return *value;
        }

        /** \brief Sets parameters, which will change all previously defined
			functions, and functions defined in the future. **/
		void ExpressionEvaluator::SetParameters(std::map<std::string, double> const& params)
        {
            if (ParsedData != NULL)
                ParsedData->ASTMode = Node::SAVE_PARAMETERS;

            if (ParametersMap == NULL)
                ParametersMap = new map<string, double>;
            else
                ParametersMap->clear();

            for (map<string, double>::const_iterator it = params.begin(); it != params.end(); it++)
                (*ParametersMap)[it->first] = it->second;
        }

        /** \brief Sets s parameter, which will change all previously defined
			functions, and functions defined in the future. **/
		void ExpressionEvaluator::SetParameter(std::string const& name, double value)
        {
            if (ParsedData != NULL)
                ParsedData->ASTMode = Node::SAVE_PARAMETERS;

            if (ParametersMap == NULL)
                ParametersMap = new map<string, double>;

            (*ParametersMap)[name] = value;
        }


		/** \brief Returns the value of a parameter. **/
		double ExpressionEvaluator::GetParameter(std::string const& name)
        {
            map<string, double>::iterator it;
            if (ParametersMap == NULL || (it = ParametersMap->find(name)) == ParametersMap->end())
            {
                throw "Parameter not found: " + name;
                return -1;
            }

            return it->second;
        }


		/** \brief Evaluates a function and returns the result. **/
        double ExpressionEvaluator::Evaluate(double start, ...)
        {
            va_list ap;

            va_start(ap, start);
            if (ParsedData == NULL)
            {
                va_end(ap);
                throw std::runtime_error("Unable to evaluate because a function must first be defined with DefineFunction(...).");
                return -1;
            }
            else
            {
                if (ParsedData->NumberVariables > 0)
                {
                    ParsedData->NextMapVar->ResetIndex();
                    *(ParsedData->NextMapVar->GetNext()) = start;
                    for (string::size_type i = 1 ; i < ParsedData->NumberVariables; i++)
                        *(ParsedData->NextMapVar->GetNext()) = va_arg(ap, double);
                }
            }
            va_end(ap);

            if (ParsedData->ASTMode == Node::SAVE_PARAMETERS)
            {
                ParsedData->ASTMode = Node::USE_SAVED_PARAMETERS;
                return EvaluateExpression(ParsedData->AST, Node::SAVE_PARAMETERS);
            }
            else
            {
                return EvaluateExpression(ParsedData->AST, ParsedData->ASTMode);
            }
        }


        /** \brief Evaluates a function and returns the result. **/
        double ExpressionEvaluator::Evaluate(std::string const& function, double start, ...)
        {
            va_list ap;

            va_start(ap, start);

            // Find the previous parsing.
            if (ParsedData == NULL)
            {
                va_end(ap);
                throw std::runtime_error("Unable to evaluate because a function must first be defined with DefineFunction(...).");
                return -1;
            }
            else
            {
                map<string, int>::const_iterator it = ParsedMapExprToId.find(function);
                if (it == ParsedMapExprToId.end())
                {
                    va_end(ap);
                    throw std::runtime_error("Unable to evaluate because a function must first be defined with DefineFunction(...).");
                    return -1;
                }

                ParsedData = ParsedASTs[it->second];
                ParsedData->ASTMode = Node::SAVE_PARAMETERS;
            }

            if (ParsedData->NumberVariables > 0)
            {
                ParsedData->NextMapVar->ResetIndex();
                *(ParsedData->NextMapVar->GetNext()) = start;
                for (string::size_type i = 1 ; i < ParsedData->NumberVariables; i++)
                    *(ParsedData->NextMapVar->GetNext()) = va_arg(ap, double);
            }

            va_end(ap);

            if (ParsedData->ASTMode == Node::SAVE_PARAMETERS)
            {
                ParsedData->ASTMode = Node::USE_SAVED_PARAMETERS;
                return EvaluateExpression(ParsedData->AST, Node::SAVE_PARAMETERS);
            }
            else
            {
                return EvaluateExpression(ParsedData->AST, ParsedData->ASTMode);
            }
        }

        /** \brief Evaluates a function and returns the result. **/
        double ExpressionEvaluator::Evaluate(const int expression_id, double start, ...)
        {
            va_list ap;

            va_start(ap, start);

            if (expression_id >= ParsedASTs.size())
            {
                va_end(ap);
                throw std::runtime_error("Unable to evaluate because a function must first be defined with DefineFunction(...).");
                return -1;
            }

            ParsedData = ParsedASTs[expression_id];
            ParsedData->ASTMode = Node::SAVE_PARAMETERS;


            if (ParsedData->NumberVariables > 0)
            {
                ParsedData->NextMapVar->ResetIndex();
                *(ParsedData->NextMapVar->GetNext()) = start;
                for (string::size_type i = 1 ; i < ParsedData->NumberVariables; i++)
                    *(ParsedData->NextMapVar->GetNext()) = va_arg(ap, double);
            }

            va_end(ap);

            if (ParsedData->ASTMode == Node::SAVE_PARAMETERS)
            {
                ParsedData->ASTMode = Node::USE_SAVED_PARAMETERS;
                return EvaluateExpression(ParsedData->AST, Node::SAVE_PARAMETERS);
            }
            else
            {
                return EvaluateExpression(ParsedData->AST, ParsedData->ASTMode);
            }
        }


        /** \brief Evaluates a function with different values, and returns the result in a vector. **/
        vector<double> ExpressionEvaluator::Evaluate(vector<double> const* start, ...)
        {
            va_list ap;
            vector<vector<double> const*> vectors;
            vector<double> rtn;

            va_start(ap, start);
            if (ParsedData == NULL)
            {
                va_end(ap);
                throw std::runtime_error("Unable to evaluate because a function must first be defined with DefineFunction(...).");
                return rtn;
            }
            else
            {
                vectors.push_back(start);
                for (string::size_type i = 1 ; i < ParsedData->NumberVariables; i++)
                    vectors.push_back(va_arg(ap, vector<double> const*));
                va_end(ap);
            }

            vector<vector<double> const*>::size_type size = vectors.size();
            if (size != ParsedData->NumberVariables)
            {
                throw std::runtime_error("The input vectors must all be the same length as the number of function variables.");
                return rtn;
            }

            vector<double>::size_type entries = vectors[0]->size();
            for (vector<vector<double> const*>::iterator it = vectors.begin()+1; it != vectors.end(); it++)
            {
                if ((*it)->size() != entries)
                {
                    throw std::runtime_error("The input vectors must all be the same length.");
                    return rtn;
                }
            }

            // Do the work of constructing the return array. For the first run through, save the parameters
            // in the DoubleValue variable for each node, and then use that value all of the other run
            // throughs since it won't change.
            ParsedData->NextMapVar->ResetIndex();
            vector<double>::size_type i;
            vector<vector<double> const*>::size_type j;
            if (entries > 0)
            {
                for (j = 0; j < size; j++)
                    *(ParsedData->NextMapVar->GetNext()) = (*vectors[j])[0];
                rtn.push_back( EvaluateExpression(ParsedData->AST, Node::SAVE_PARAMETERS) );
                ParsedData->ASTMode = Node::USE_SAVED_PARAMETERS;
            }
            for (i = 1; i < entries; i++)
            {
                for (j = 0; j < size; j++)
                    *(ParsedData->NextMapVar->GetNext()) = (*vectors[j])[i];
                rtn.push_back( EvaluateExpression(ParsedData->AST, ParsedData->ASTMode) );
            }

            return rtn;
        }

        /** Gcc 4.0.1 on the mac did not compile when isspace was passed as the third 
			parameter to std::remove_if below (it takes an int as a parameter 
			instead of a char).  This wrapper method gets around this problem. **/
        static int ConvertIsSpaceForGcc(char c) { return isspace(c); }
        
        /** \brief Creates a simplified AST from the Spirit parse. **/
		ExpressionEvaluator::Node* ExpressionEvaluator::CreateAST(boost_spirit::tree_match<std::string::const_iterator, boost_spirit::node_val_data_factory<double> >::tree_iterator const &i)
		{
			string valueStr(i->value.begin(), i->value.end());
			valueStr.erase(std::remove_if(valueStr.begin(), valueStr.end(), ConvertIsSpaceForGcc), valueStr.end());	// trim
			const parser_id parserID = i->value.id();
			if (parserID == MathExpression::constantID)
			{
				if (i->children.size() != 0)
				{
					throw std::runtime_error("Illegal children under constant node: " + valueStr);
					return NULL;
				}

				double* value = find(*constants_p, valueStr.c_str() );
				if (value == NULL)
				{
					throw std::runtime_error("Cannot find the value for the specified constant: " + valueStr);
					return NULL;
				}

				Node* n = new Node;
				n->DataType = Node::CONSTANT;
				n->DoubleValue = *value;
				return n;
			}
			else if (parserID == MathExpression::numberID)
			{
				if (i->children.size() != 0)
				{
					throw std::runtime_error("Illegal children under number node: " + valueStr);
					return NULL;
				}

				Node* n = new Node;
				n->DataType = Node::DOUBLE;
				//n->DoubleValue = i->value.value();
				n->DoubleValue = atof(valueStr.c_str());
				return n;
			}
			else if (parserID == MathExpression::variableID)
			{
				if (i->children.size() != 0)
				{
					throw std::runtime_error("Illegal children under variable node: " + valueStr);
					return NULL;
				}

				map<string, double>::iterator it;
				if (ParsedData->VariableMap == NULL ||
						(it = ParsedData->VariableMap->find(valueStr)) == ParsedData->VariableMap->end())
				{
					throw std::runtime_error("Unknown variable parsed: " + valueStr);
					return NULL;
				}

				Node* n = new Node;
				n->DataType = Node::VARIABLE;
				n->StringValue = new string(valueStr);
				n->PointerToVariable = &it->second;
				return n;
			}
			else if (parserID == MathExpression::parameterID)
			{
				if (i->children.size() != 0)
				{
					throw std::runtime_error("Illegal children under parameter node: " + valueStr);
					return NULL;
				}

				Node* n = new Node;
				n->DataType = Node::PARAMETER;
				n->StringValue = new string(valueStr);
				return n;
			}
			else if (parserID == MathExpression::functionID)
			{
				func* funcptr = find(functions_p, valueStr.c_str());
				if (funcptr == NULL)
				{
					throw std::runtime_error("Invalid function specified: " + valueStr);
					return NULL;
				}
				if (funcptr->size != i->children.size())
				{
					throw std::runtime_error("Illegal number or arguments for math function: " + valueStr);
					return NULL;
				}

				bool allDoubles = true;
				vector<Node*>* arguments = new vector<Node*>;
				for (tree_match<string::const_iterator, node_val_data_factory<double> >::tree_iterator it = i->children.begin(); it != i->children.end(); it++)
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

					CheckMathOperationForErrors(valueStr);

					return node;
				}
				else
				{
					Node* node = new Node;
					node->DataType = Node::FUNCTION;
					node->StringValue = new string(valueStr);
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
				if (*valueStr.begin() != '-')
				{
					throw std::runtime_error("Illegal factor - it can only be '-' and it was: " + valueStr);
					return NULL;
				}
				if (i->children.size() != 1)
				{
					throw std::runtime_error("Illegal number of children under factor node: " + valueStr);
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
			else if (parserID == MathExpression::operatorID)
			{
				if (i->children.size() != 2)
				{
					throw std::runtime_error("Too many arguments for mathematical operator: " + valueStr);
					return NULL;
				}

				Node* left = CreateAST(i->children.begin()+0);
				Node* right = CreateAST(i->children.begin()+1);

				// If they are both of double type, we can simplify to a number.
				if (left->DataType == Node::DOUBLE && right->DataType == Node::DOUBLE)
				{
					Node* node = new Node;
					node->DataType = Node::DOUBLE;
					switch(*valueStr.begin())
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
					case '<':
						if (*(valueStr.end()-1) == '=')
							node->DoubleValue = left->DoubleValue <= right->DoubleValue;
						else
							node->DoubleValue = left->DoubleValue < right->DoubleValue;
						break;
					case '>':
						if (*(valueStr.end()-1) == '=')
							node->DoubleValue = left->DoubleValue >= right->DoubleValue;
						else
							node->DoubleValue = left->DoubleValue > right->DoubleValue;
						break;
					case '=':
						node->DoubleValue = left->DoubleValue == right->DoubleValue;
						break;
					case '&':
						node->DoubleValue = left->DoubleValue && right->DoubleValue;
						break;
					case '|':
						node->DoubleValue = left->DoubleValue || right->DoubleValue;
						break;

					default:
						delete left;
						delete right;
						delete node;
						throw "Invalid operator encountered: " + valueStr;
						return NULL;
					}

					delete left, right;
					return node;
				}
				else
				{
					Node* node = new Node;
					switch(*valueStr.begin())
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
					case '<':
						if (*(valueStr.end()-1) == '=')
							node->DataType = Node::LTEQ;
						else
							node->DataType = Node::LT;
						break;
					case '>':
						if (*(valueStr.end()-1) == '=')
							node->DataType = Node::GTEQ;
						else
							node->DataType = Node::GT;
						break;
					case '=':
						node->DataType = Node::EQUAL;
						break;
					case '&':
						node->DataType = Node::LOGICAL_AND;
						break;
					case '|':
						node->DataType = Node::LOGICAL_OR;
						break;
					default:
						delete left;
						delete right;
						delete node;
						throw "Invalid operator encountered: " + valueStr;
						return NULL;
					}

					node->children = new vector<Node*>;
					node->children->push_back(left);
					node->children->push_back(right);
					return node;
				}
			}

			throw "Illegal expression encountered: " + valueStr;
			return NULL;
		}


		/** \brief Evaluates the AST, returning a double result. **/
        double ExpressionEvaluator::EvaluateExpression(Node* const n, Node::AST_Mode mode)
        {
            switch (n->DataType)
            {
            case Node::CONSTANT:
            case Node::DOUBLE:
                return n->DoubleValue;
            case Node::VARIABLE:
                return *n->PointerToVariable;	// this is a pointer to the entry in the map
            case Node::PARAMETER:
                if (mode == Node::USE_SAVED_PARAMETERS)
                {
                    return n->DoubleValue;
                }
                else
                {
                    map<string, double>::iterator it;
                    if (ParametersMap == NULL || (it = ParametersMap->find(*n->StringValue)) == ParametersMap->end())
                    {
                        throw std::runtime_error("Illegal parameter specified: " + *n->StringValue);
                        return -1;
                    }
                    if (mode == Node::SAVE_PARAMETERS) n->DoubleValue = it->second;
                    return it->second;
                }
            case Node::FUNCTION:
                {
                    errno = 0;
                    double result = 0;
                    switch (n->children->size())
                    {
                    case 1:
                        result = (*n->Function1)( EvaluateExpression((*n->children)[0], mode) );
                        break;
                    case 2:
                        result = (*n->Function2)( EvaluateExpression((*n->children)[0], mode), EvaluateExpression((*n->children)[1], mode) );
                        break;
                    case 3:
                        result = (*n->Function3)( EvaluateExpression((*n->children)[0], mode), EvaluateExpression((*n->children)[1], mode), EvaluateExpression((*n->children)[2], mode) );
                        break;
                    case 4:
                        result = (*n->Function4)( EvaluateExpression((*n->children)[0], mode), EvaluateExpression((*n->children)[1], mode), EvaluateExpression((*n->children)[2], mode), EvaluateExpression((*n->children)[3], mode) );
                        break;
                    }
                    CheckMathOperationForErrors(*n->StringValue);
                    return result;
                }
            case Node::ADDITION:
                return EvaluateExpression((*n->children)[0], mode) + EvaluateExpression((*n->children)[1], mode);
            case Node::SUBTRACTION:
                return EvaluateExpression((*n->children)[0], mode) - EvaluateExpression((*n->children)[1], mode);
            case Node::MULTIPLICATION:
                return EvaluateExpression((*n->children)[0], mode) * EvaluateExpression((*n->children)[1], mode);
            case Node::DIVISION:
                return EvaluateExpression((*n->children)[0], mode) / EvaluateExpression((*n->children)[1], mode);
            case Node::EXPONENTIATION:
                return pow(EvaluateExpression((*n->children)[0], mode), EvaluateExpression((*n->children)[1], mode));
            case Node::LT:
                return EvaluateExpression((*n->children)[0], mode) < EvaluateExpression((*n->children)[1], mode);
            case Node::LTEQ:
                return EvaluateExpression((*n->children)[0], mode) <= EvaluateExpression((*n->children)[1], mode);
            case Node::GT:
                return EvaluateExpression((*n->children)[0], mode) > EvaluateExpression((*n->children)[1], mode);
            case Node::GTEQ:
                return EvaluateExpression((*n->children)[0], mode) >= EvaluateExpression((*n->children)[1], mode);
            case Node::EQUAL:
                return EvaluateExpression((*n->children)[0], mode) == EvaluateExpression((*n->children)[1], mode);
            case Node::LOGICAL_AND:
                return EvaluateExpression((*n->children)[0], mode) && EvaluateExpression((*n->children)[1], mode);
            case Node::LOGICAL_OR:
                return EvaluateExpression((*n->children)[0], mode) || EvaluateExpression((*n->children)[1], mode);
            }

            throw std::runtime_error("Illegal expression was encountered.");
            return -1;
        }

#if 0

        // ---------------------------------------------------------------------
        //             Math Functions that are Used when Parsing Input
        // ---------------------------------------------------------------------

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

#endif

    };
};
