/**
 * \page pageCodeStandard Coding Standards Guide
 * The purpose of this page is to detail the coding standards of the project
 * which all contributers are requested to follow.
 *
 * This page describes the coding style standard for C++. A coding style
 * standard defines the visual layout of source code.  Presenting source code
 * in a uniform fashion facilitates the use of code by different developers.
 * In addition, following a standard prevents certain types of coding errors.
 *
 * All of the items below, unless otherwise noted, are guidelines. They are
 * recommendations about how to lay out a given block of code. Use common
 * sense and provide comments to describe any deviation from the standard.
 * Sometimes, violating a guideline may actually improve readability.
 *
 * If you are working with code that does not follow the standard, bring the
 * code up-to-date or follow the existing style.  Don’t mix styles.
 *
 * \section secCodeStandardLayout Code Layout
 * The aim here is to maximise readability on all platforms and editors.
 * - Code width of 80 characters maximum - hard-wrap longer lines.
 * - Use sensible wrapping for long statements in a way which maximises
 *   readability.
 * - Do not put multiple statements on the same line.
 * - Do not declare multiple variables on the same line.
 * - Provide a default value on all variable declarations.
 * - Enclose every program block (if, else, for, while, etc) in braces, even if
 *   empty or just a single line.
 * - Opening braces ({) should be on their own line.
 * - Braces at same indentation as preceeding statement.
 * - One class per .cpp and .h file only, unless nested.
 * - Define member functions in the .cpp file in the same order as defined in
 *   the .h file.
 * - Templated classes defined and implemented in a single .hpp file.
 * - Do not put inline functions in the header file unless the function is
 *   trivial (e.g. accessor, empty destructor), or profiling explicitly suggests
 *   to.
 * - Inline functions should be declared within the class declaration but
 *   defined outside the class declaration at the bottom of the header file.
 *   <b>Note:</b><tt>virtual</tt> and <tt>inline</tt> are mutually exclusive.
 *   Virtual functions should therefore be implemented in the .cpp file.
 *
 * \section White Space
 * Adding an appropriate amount of white space enhances readability.
 * Too much white space, on the other hand, detracts from that readability.
 * - Indent using a four-space tab.  Consistent tab spacing is necessary to
 *   maintain formatting.  Note that this means when a tab is pressed, four
 *   physical spaces are inserted into the source instead.
 * - Put a blank line at the end of a public/protected/private block.
 * - Put a blank line at the end of every file.
 * - Put a space after every keyword (if, while, for, etc.).
 * - Put a space after every comma, unless the comma is at the end of the line.
 * - Do not put a space before the opening parenthesis of an argument list to
 *   a function.
 * - Declare pointers and references with the * or & symbol next to the
 *   declarator, not the type; e.g., Object *object.  Do not put multiple
 *   variables in the same declaration.
 * - Place a space on both sides of a binary operator.
 * - Do not use a space to separate a unary operator from its operand.
 * - Place open and close braces on their own line.  No executable statements
 *   should appear on the line with the brace, but comments are allowed.
 *   Indent opening braces at the same level as the statement above and indent
 *   the closing brace at the same level as the corresponding opening brace.
 * - Indent all statements following an open brace by one tab. Developer Studio
 *   puts any specifier terminated with a colon at the same indentation level
 *   as the enclosing brace.  Examples of such specifiers include case
 *   statements, access specifiers (public, private, protected), and goto
 *   labels.  This is not acceptable and should be manually corrected so that
 *   all statements appearing within a block and delineated by braces are
 *   indented.
 * - Break a line into multiple lines when it becomes too long to read.  Use
 *   at least two tabs to start the new line, so it does not look like the
 *   start of a block.
 * - Follow C++ style comments with one space. It is also preferable to
 *   consider any text that follows C++ style comments as a sentence and to
 *   begin this text with a capital letter.  This helps to distinguish the
 *   line from a continuation of a previous line; i.e., // This is my comment.
 * - As a general rule, don’t keep commented out source code in the final
 *   baselined product.  Such code leads the reader to believe there was
 *   uncertainty in the code as it currently exists.
 * - Place the # of a preprocessor directive at column one.  An exception is
 *   the use of nested #ifdefs where the bodies only contain other preprocessor
 *   directives.  Add tabs to enhance readability.
 * @begin{code}
 * void foo()
 * {
 *   for(int I = 0; I < 10; ++i)
 *   {
 * #ifdef BAR
 *     do_something();
 * #endif
 *     for_loop_code();
 *    }
 * }
 * @end{code}
 * - Use tabular white space if it enhances readability.
 * - Use only one return statement.  Structure the code so that only one return
 *   statement is necessary.
 *
 *
 * \section secCodeStandardNaming Naming Conventions
 * Keep variable and function names meaningful but concise.
 * - Begin variable names with lower-case letter.
 * - Begin function names and class names with upper-case letter.
 * - All function, class and variable names should be written in CamelCase, e.g.
 *   MyClass, DoFunction() or myVariableName.
 * - All preprocessor definitions written in UPPER_CASE with words separated by
 *   underscores, e.g. USE_SPECIFIC_FEATURE.
 * - All member variables prefixed with m_.
 * - All constants prefixed with a k.
 * - All function parameters prefixed with a p.
 * - All enumerations prefixed with an e.
 * - Do not use leading underscores.
 *
 *
 * \section secCodeStandardNamespaces Namespaces
 * The top-level namespace is "Nektar". All code should reside in this namespace
 * or a sub-space of this.
 * - Namespaces correspond to code structure.
 * - Namespaces should be kept to a minimum to simplify the interface to their
 *   contents.
 *
 *
 * \section secCodeStandardDocumentation Documentation
 * - Briefs for classes, functions and types in header files using /// notation.
 * - Full documentation with implementation using /** ... *\/ notation.
 * - Use \@ symbol for \@class, \@param, \@returns, etc for ease of identification.
 * - Any separate documentation pages not directly associated with a portion of
 * the code should be in a separate file in /docs/html/doxygen.
 */
