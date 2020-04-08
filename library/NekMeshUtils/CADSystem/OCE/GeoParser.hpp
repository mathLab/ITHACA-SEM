////////////////////////////////////////////////////////////////////////////////
//
//  File: GeoParser.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Simple parser for the .geo format based on boost::spirit.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <LibUtilities/Interpreter/Interpreter.h>

namespace Nektar
{
namespace NekMeshUtils
{

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

/**
 * @brief Simple AST for the Geo format which will be populated by the
 * #GeoParser class.
 */
namespace GeoAst
{

/**
 * @brief Evaluate an expression using the LibUtilities::Interpreter class.
 *
 * @param interp  Handle to the interpreter to use.
 * @param expr    String expression returned from boost::qi, handled as a vector
 *                of characters.
 *
 * @return The double-precision evaluation of the expression.
 */
double Eval(LibUtilities::Interpreter &interp,
            std::vector<char> &expr)
{
    std::string exprStr(expr.begin(), expr.end());
    int interpId = interp.DefineFunction("", exprStr);
    return interp.Evaluate(interpId);
}

/**
 * @brief Set a parameter inside the @p interp object. Used to store Gmsh
 * variables defined in the .geo file.
 *
 * @param interp  Handle to the interpreter to use.
 * @param expr    Variable name returned from boost::qi, handled as a vector
 *                of characters.
 * @param val     Value of the parameter.
 */
void SetParam(LibUtilities::Interpreter &interp,
              std::vector<char> &varname,
              double &val)
{
    std::string varStr(varname.begin(), varname.end());
    interp.SetParameter(varStr, val);
}

/**
 * @brief Print a warning for unknown geometry types.
 *
 * @param geomname  The interpreted geomtery name.
 */
void PrintWarning(std::vector<char> &geomname)
{
    std::string geom(geomname.begin(), geomname.end());
    std::cout << "Warning: ignoring unknown geometry entity"
              << geom << std::endl;
}

/**
 * @brief Wrapper for the Gmsh Point geometry type.
 */
struct Point
{
    /// Point ID.
    unsigned int id;
    /// Point x-coordinate.
    double x;
    /// Point y-coordinate.
    double y;
    /// Point z-coordinate.
    double z;
    /// Characteristic length at this point.
    double clen;
};

/**
 * @brief Wrapper for all other Geometry types, which are defined using a list
 * of integer IDs. For example, a line is defined from the IDs of two points.
 */
struct Geom
{
    /// ID of the geometry object.
    unsigned int id;
    /// List of IDs that define the geometry object.
    std::vector<int> ids;
};

/**
 * @brief Lightweight wrapper for a .geo file structure.
 */
struct GeoFile
{
    /// Vector of points.
    std::vector<Point> points;
    /// Vector of lines.
    std::vector<Geom> lines;
    /// Vector of splines.
    std::vector<Geom> splines;
    /// Vector of bsplines.
    std::vector<Geom> bsplines;
    /// Vector of circles
    std::vector<Geom> circles;
    /// Vector of ellipses.
    std::vector<Geom> ellipses;
    /// Vector of line loops.
    std::vector<Geom> lineLoops;
    /// Vector of plane surfaces.
    std::vector<Geom> planeSurfs;
    /// Vector of ruled surfaces (or 'surfaces' in the newer .geo formats).
    std::vector<Geom> ruledSurfs;
    /// Vector of surface loops.
    std::vector<Geom> surfLoops;
    /// Vector of volumes.
    std::vector<Geom> volumes;
};

}

/**
 * @brief A skipper that boost::qi can use to ignore all comments in the .geo
 * file. We support the use of single and multi-line comments.
 *
 * @tparam Iterator   Iterator type, e.g. std::string::const_iterator.
 */
template<typename Iterator>
struct CommentSkipper : public qi::grammar<Iterator>
{
    /**
     * @brief Constructor for a simple comment skipper grammar.
     */
    CommentSkipper() : CommentSkipper::base_type(skip)
    {
        // Inline comments.
        line = "//" >> *(qi::char_ - qi::eol) >> (qi::eol | qi::eoi);
        // Block comments.
        block = "/*" >> *(block | (qi::char_ - "*/")) > "*/";
        // Also skip all whitespace.
        skip = boost::spirit::ascii::space | line | block;
    }

    /// Grammar rules.
    qi::rule<Iterator> skip, block, line;
};

/**
 * @brief A lightweight parser for the .geo format.
 *
 * This class defines a grammar that parses a small subset of the .geo
 * format. Presently, it supports:
 *
 * - 0D: points;
 * - 1D: lines, splines, bsplines, circles, ellipses and line loops;
 * - 2D: surfaces and ruled surfaces;
 * - 3D: volumes
 *
 * In order to parse mathematical expressions, we use the
 * LibUtilities::Interpreter class, which can then be used to evaluate and store
 * .geo variables also.
 */
template <typename Iterator, typename Skipper = CommentSkipper<Iterator>>
struct GeoParser : qi::grammar<Iterator, GeoAst::GeoFile(), Skipper>
{
    using Interp = LibUtilities::Interpreter;

    GeoParser(Interp &interp) : GeoParser::base_type(geoFile), m_interp(interp)
    {
        using GeoAst::GeoFile;
        using GeoAst::Geom;
        using GeoAst::Point;
        using phx::bind;
        using phx::push_back;
        using qi::_1;
        using qi::_2;
        using qi::_val;

        // Defines the rules for arithmetic expressions. These are then parsed
        // by the Nektar++ interpreter since we are a bit lazy.
        expr = (*(
            qi::alnum | qi::char_('+') | qi::char_('-') | qi::char_('_') |
            qi::char_('*') | qi::char_('/') | qi::char_('.') | qi::char_('(') |
            qi::char_(')')))[_val = phx::bind(GeoAst::Eval, phx::ref(interp), _1)];

        // Defines the rule for a point.
        point = '(' >> qi::uint_[bind(&Point::id, _val) = _1] >> ')' >> '=' >>
                '{' >> (qi::double_ | expr)[bind(&Point::x, _val) = _1] >>
                ',' >> (qi::double_ | expr)[bind(&Point::y, _val) = _1] >>
                ',' >> (qi::double_ | expr)[bind(&Point::z, _val) = _1] >>
                ',' >> (qi::double_ | expr)[bind(&Point::clen, _val) = _1] >>
                '}';

        // Defines the rule for a geometry object. This approach parses
        // everything after the geometry name, which we use to make the grammar
        // a bit simpler.
        geom = "(" >> qi::uint_[bind(&Geom::id, _val) = _1] >> ")" >> "=" >>
               "{" >> (qi::int_[push_back(bind(&Geom::ids, _val), _1)] % ",") >>
               "}";

        // Define the rule for the geo file itself.
        geoFile = *(
            // Geometry types
            (("Point" >> point[push_back(bind(&GeoFile::points, _val), _1)]) |
             ("Line" >> geom[push_back(bind(&GeoFile::lines, _val), _1)]) |
             ("Spline" >> geom[push_back(bind(&GeoFile::splines, _val), _1)]) |
             ((qi::lit("BSpline") | "Bezier") >>
              geom[push_back(bind(&GeoFile::bsplines, _val), _1)]) |
             ("Circle" >> geom[push_back(bind(&GeoFile::circles, _val), _1)]) |
             ("Ellipse" >>
              geom[push_back(bind(&GeoFile::ellipses, _val), _1)]) |
             ((qi::lit("Line Loop") | "Curve Loop") >>
              geom[push_back(bind(&GeoFile::lineLoops, _val), _1)]) |
             ("Plane Surface" >>
              geom[push_back(bind(&GeoFile::planeSurfs, _val), _1)]) |
             ((qi::lit("Surface") | "Ruled Surface") >>
              geom[push_back(bind(&GeoFile::ruledSurfs, _val), _1)]) |
             ("Surface Loop" >>
              geom[push_back(bind(&GeoFile::surfLoops, _val), _1)]) |
             ("Volume" >> geom[push_back(bind(&GeoFile::volumes, _val), _1)]) |

             // Unknown geometry type
             ((*qi::alpha >> geom)[phx::bind(&GeoAst::PrintWarning, _1)]) |

             // Variables
             (*(qi::alnum | qi::char_('_') | qi::char_('.')) >> '=' >>
              expr)[phx::bind(&GeoAst::SetParam, phx::ref(interp), _1, _2)]) >>
            ";");
    }

    template <typename T> using rule = qi::rule<Iterator, T(), Skipper>;

    /// Rule for points.
    rule<GeoAst::Point> point;
    /// Rule for geometry objects.
    rule<GeoAst::Geom> geom;
    /// Rule for geo file.
    rule<GeoAst::GeoFile> geoFile;
    /// Rule for expressions
    rule<double> expr;
    /// Reference to an interpreter that is used to .
    Interp &m_interp;
};
}
}
