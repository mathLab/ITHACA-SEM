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

namespace Nektar
{
namespace NekMeshUtils
{

using namespace boost::spirit;
namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

namespace ast
{

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
    std::vector<unsigned int> ids;
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
        skip = ascii::space | line | block;
    }
    qi::rule<Iterator> skip, block, line;
};

template <typename Iterator, typename Skipper = CommentSkipper<Iterator>>
struct GeoParser : qi::grammar<Iterator, ast::GeoFile(), Skipper>
{
    GeoParser() : GeoParser::base_type(geoFile)
    {
        using phx::push_back;
        using phx::bind;

        using ast::Point;
        using ast::Geom;
        using ast::GeoFile;

        point = '(' >> qi::uint_[bind(&Point::id, _val) = _1] >> ')'
            >> '=' >> '{'
            >> qi::double_[bind(&Point::x, _val) = _1] >> ','
            >> qi::double_[bind(&Point::y, _val) = _1] >> ','
            >> qi::double_[bind(&Point::z, _val) = _1] >> ','
            >> qi::double_[bind(&Point::clen, _val) = _1] >> '}';

        geom = "(" >> qi::uint_[bind(&Geom::id, _val) = _1] >> ")"
            >> "=" >> "{"
            >> (qi::uint_[push_back(bind(&Geom::ids, _val), _1)] % ",")
            >> "}";

        geoFile = *(
            (
                ("Point" >> point[push_back(bind(&GeoFile::points, _val), _1)]) |
                ("Line"  >> geom[push_back(bind(&GeoFile::lines, _val), _1)]) |
                ("Spline"  >> geom[push_back(bind(&GeoFile::splines, _val), _1)]) |
                ("BSpline"  >> geom[push_back(bind(&GeoFile::bsplines, _val), _1)]) |
                ("Circle"  >> geom[push_back(bind(&GeoFile::circles, _val), _1)]) |
                ("Ellipse"  >> geom[push_back(bind(&GeoFile::ellipses, _val), _1)]) |
                ("Line Loop" >> geom[push_back(bind(&GeoFile::lineLoops, _val), _1)]) |
                ("Plane Surface"  >> geom[push_back(bind(&GeoFile::planeSurfs, _val), _1)]) |
                ("Ruled Surface"  >> geom[push_back(bind(&GeoFile::ruledSurfs, _val), _1)]) |
                ("Surface Loop"  >> geom[push_back(bind(&GeoFile::surfLoops, _val), _1)]) |
                ("Volume"  >> geom[push_back(bind(&GeoFile::volumes, _val), _1)])
                ) >> ";");
    }

    template<typename T>
    using rule = qi::rule<Iterator, T(), Skipper>;

    rule<ast::Point> point;
    rule<ast::Geom> geom;
    rule<ast::Geom> line;
    rule<ast::Geom> spline;
    rule<ast::Geom> bspline;
    rule<ast::Geom> circle;
    rule<ast::Geom> ellipse;
    rule<ast::Geom> lineLoop;
    rule<ast::Geom> planeSurf;
    rule<ast::Geom> ruledSurf;
    rule<ast::Geom> surfLoop;
    rule<ast::Geom> volume;
    rule<ast::GeoFile> geoFile;
};


}
}
