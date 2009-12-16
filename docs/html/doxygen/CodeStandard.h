/**
 * \page pageCodeStandard Coding Standards Guide
 * The purpose of this page is to detail the coding standards of the project
 * which all contributers are requested to follow.
 *
 * This document describes the coding style standard for C++. A coding style
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
 * Portions of this standard associated with the use of Doxygen for automatic
 * documentation generation or SourceSafe for configuration management may be
 * tailored to suit the particular needs of projects that do not use these
 * tools.  An attempt should be made, however, to adhere to the standard using
 * similar features from other tools.
 *
 * \section secCodeStandardLayout Code Layout
 * The aim here is to maximise readability on all platforms.
 * - Code width of 80 characters maximum.
 * - Use sensible wrapping for long statements in a way which maximises
 *   readability.
 * - Opening braces ({) should be on their own line.
 * - Braces at same indentation as preceeding statement.
 * - Inline functions should be declared within the class declaration but
 *   defined outside the class declaration at the bottom of the header file.
 *   <b>Note:</b><tt>virtual</tt> and <tt>inline</tt> are mutually exclusive in
 *   almost all circumstances.
 *
 * \section secCodeStandardNaming Naming Conventions
 * Keep variable and function names meaningful but concise. 
 * - All member variables prefixed with m and use CamelCase. e.g. mMemVar. 
 * - All function names using CamelCase
 *
 * \section secCodeStandardDocumentation Documentation
 * - Briefs for classes, functions and types in header files using /// notation.
 * - Full documentation with implementation using /** ... *\/ notation.
 * - Use @ symbol for @class, @param, @returns, etc for ease of identification.
 * - Any separate documentation pages not directly associated with a portion of
 * the code should be in a separate file in /docs/html/doxygen.
 */
