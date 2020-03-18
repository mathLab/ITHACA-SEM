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

namespace ast
{

double Eval(LibUtilities::Interpreter &interp,
            std::vector<char> &expr)
{
    std::string exprStr(expr.begin(), expr.end());
    int interpId = interp.DefineFunction("", exprStr);
    return interp.Evaluate(interpId);
}

void SetParam(LibUtilities::Interpreter &interp,
              std::vector<char> &varname,
              double &val)
{
    std::string varStr(varname.begin(), varname.end());
    interp.SetParameter(varStr, val);
}

void PrintWarning(std::vector<char> &geomname)
{
    std::string geom(geomname.begin(), geomname.end());
    std::cout << "Warning: ignoring unknown geometry entity"
              << geom << std::endl;
}

struct Point
{
    unsigned int id;
    double x;
    double y;
    double z;
    double clen;
};

struct Geom
{
    unsigned int id;
    std::vector<int> ids;
};

struct GeoFile
{
    std::vector<Point> points;
    std::vector<Geom> lines;
    std::vector<Geom> splines;
    std::vector<Geom> bsplines;
    std::vector<Geom> circles;
    std::vector<Geom> ellipses;
    std::vector<Geom> lineLoops;
    std::vector<Geom> planeSurfs;
    std::vector<Geom> ruledSurfs;
    std::vector<Geom> surfLoops;
    std::vector<Geom> volumes;
};

}

template<typename Iterator>
struct CommentSkipper : public qi::grammar<Iterator> {

    CommentSkipper() : CommentSkipper::base_type(skip)
    {
        line = "//" >> *(qi::char_ - qi::eol) >> (qi::eol | qi::eoi);
        block = "/*" >> *(block | qi::char_ - "*/") > "*/";
        skip = boost::spirit::ascii::space | line | block;
    }
    qi::rule<Iterator> skip, block, line;
};

template <typename Iterator, typename Skipper = CommentSkipper<Iterator>>
struct GeoParser : qi::grammar<Iterator, ast::GeoFile(), Skipper>
{
    using Interp = LibUtilities::Interpreter;

    GeoParser(Interp &interp) : GeoParser::base_type(geoFile), m_interp(interp)
    {
        using phx::push_back;
        using phx::bind;
        using qi::_1;
        using qi::_2;
        using qi::_val;
        using ast::Point;
        using ast::Geom;
        using ast::GeoFile;

        expr = (*(qi::alnum | qi::char_('+') | qi::char_('-') | qi::char_('_') |
                  qi::char_('*') | qi::char_('/') | qi::char_('.') |
                  qi::char_('(') | qi::char_(')')))[
            _val = phx::bind(ast::Eval, phx::ref(interp), _1)
            ];

        point = '('
            >> qi::uint_[bind(&Point::id, _val) = _1] >> ')'
            >> '=' >> '{'
            >> (qi::double_ | expr) [bind(&Point::x, _val) = _1] >> ','
            >> (qi::double_ | expr) [bind(&Point::y, _val) = _1] >> ','
            >> (qi::double_ | expr) [bind(&Point::z, _val) = _1] >> ','
            >> (qi::double_ | expr) [bind(&Point::clen, _val) = _1] >> '}';

        geom = "("
            >> qi::uint_[bind(&Geom::id, _val) = _1] >> ")" >> "=" >> "{"
            >> (qi::int_[push_back(bind(&Geom::ids, _val), _1)] % ",")
            >> "}";

        geoFile = *(
            (
                ("Point" >> point[push_back(bind(&GeoFile::points, _val), _1)]) |
                ("Line"  >> geom[push_back(bind(&GeoFile::lines, _val), _1)]) |
                ("Spline"  >> geom[push_back(bind(&GeoFile::splines, _val), _1)]) |
                (
                    (qi::lit("BSpline") | "Bezier")
                    >> geom[push_back(bind(&GeoFile::bsplines, _val), _1)]) |
                ("Circle"  >> geom[push_back(bind(&GeoFile::circles, _val), _1)]) |
                ("Ellipse"  >> geom[push_back(bind(&GeoFile::ellipses, _val), _1)]) |
                (
                    (qi::lit("Line Loop") | "Curve Loop")
                    >> geom[push_back(bind(&GeoFile::lineLoops, _val), _1)]) |
                ("Plane Surface"  >> geom[push_back(bind(&GeoFile::planeSurfs, _val), _1)]) |
                (
                    (qi::lit("Surface") | "Ruled Surface")
                    >> geom[push_back(bind(&GeoFile::ruledSurfs, _val), _1)]) |
                ("Surface Loop"  >> geom[push_back(bind(&GeoFile::surfLoops, _val), _1)]) |
                ("Volume"  >> geom[push_back(bind(&GeoFile::volumes, _val), _1)]) |
                ((*qi::alpha >> geom) [phx::bind(&ast::PrintWarning, _1)] ) |
                (*(qi::alnum | qi::char_('_') | qi::char_('.')) >> '=' >> expr)[
                    phx::bind(&ast::SetParam, phx::ref(interp), _1, _2)
                    ]
                )
            >> ";");
    }



    template<typename T>
    using rule = qi::rule<Iterator, T(), Skipper>;

    rule<ast::Point> point;
    rule<ast::Geom> geom;
    rule<ast::GeoFile> geoFile;
    rule<double> expr;

    Interp &m_interp;
};


}
}
