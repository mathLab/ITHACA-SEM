////////////////////////////////////////////////////////////////////////////////
//
//  File: Triangle.cpp
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
//  Description: class for triangle, originally the code of Jonathan Shewchuk
//               but heavily modified.
//               original file header below
//
////////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/*                                                                           */
/*      888888888        ,o,                          / 888                  */
/*         888    88o88o  "    o8888o  88o8888o o88888o 888  o88888o         */
/*         888    888    888       88b 888  888 888 888 888 d888  88b        */
/*         888    888    888  o88^o888 888  888 "88888" 888 8888oo888        */
/*         888    888    888 C888  888 888  888  /      888 q888             */
/*         888    888    888  "88o^888 888  888 Cb      888  "88oooo"        */
/*                                              "8oo8D                       */
/*                                                                           */
/*  A Two-Dimensional Quality Mesh Generator and Delaunay Triangulator.      */
/*  (triangle.c)                                                             */
/*                                                                           */
/*  Version 1.6                                                              */
/*  July 28, 2005                                                            */
/*                                                                           */
/*  Copyright 1993, 1995, 1997, 1998, 2002, 2005                             */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*  This program may be freely redistributed under the condition that the    */
/*    copyright notices (including this entire header and the copyright      */
/*    notice printed when the `-h' switch is selected) are not removed, and  */
/*    no compensation is received.  Private, research, and institutional     */
/*    use is free.  You may distribute modified versions of this code UNDER  */
/*    THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE   */
/*    SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE   */
/*    AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR    */
/*    NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution of this code as    */
/*    part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT  */
/*    WITH THE AUTHOR.  (If you are not directly supplying this code to a    */
/*    customer, and you are instead telling them how they can obtain it for  */
/*    free, then you are not required to make any arrangement with me.)      */
/*                                                                           */
/*  Hypertext instructions for Triangle are available on the Web at          */
/*                                                                           */
/*      http://www.cs.cmu.edu/~quake/triangle.html                           */
/*                                                                           */
/*  Disclaimer:  Neither I nor Carnegie Mellon warrant this code in any way  */
/*    whatsoever.  This code is provided "as-is".  Use at your own risk.     */
/*                                                                           */
/*  Some of the references listed below are marked with an asterisk.  [*]    */
/*    These references are available for downloading from the Web page       */
/*                                                                           */
/*      http://www.cs.cmu.edu/~quake/triangle.research.html                  */
/*                                                                           */
/*  Three papers discussing aspects of Triangle are available.  A short      */
/*    overview appears in "Triangle:  Engineering a 2D Quality Mesh          */
/*    Generator and Delaunay Triangulator," in Applied Computational         */
/*    Geometry:  Towards Geometric Engineering, Ming C. Lin and Dinesh       */
/*    Manocha, editors, Lecture Notes in Computer Science volume 1148,       */
/*    pages 203-222, Springer-Verlag, Berlin, May 1996 (from the First ACM   */
/*    Workshop on Applied Computational Geometry).  [*]                      */
/*                                                                           */
/*    The algorithms are discussed in the greatest detail in "Delaunay       */
/*    Refinement Algorithms for Triangular Mesh Generation," Computational   */
/*    Geometry:  Theory and Applications 22(1-3):21-74, May 2002.  [*]       */
/*                                                                           */
/*    More detail about the data structures may be found in my dissertation: */
/*    "Delaunay Refinement Mesh Generation," Ph.D. thesis, Technical Report  */
/*    CMU-CS-97-137, School of Computer Science, Carnegie Mellon University, */
/*    Pittsburgh, Pennsylvania, 18 May 1997.  [*]                            */
/*                                                                           */
/*  Triangle was created as part of the Quake Project in the School of       */
/*    Computer Science at Carnegie Mellon University.  For further           */
/*    information, see Hesheng Bao, Jacobo Bielak, Omar Ghattas, Loukas F.   */
/*    Kallivokas, David R. O'Hallaron, Jonathan R. Shewchuk, and Jifeng Xu,  */
/*    "Large-scale Simulation of Elastic Wave Propagation in Heterogeneous   */
/*    Media on Parallel Computers," Computer Methods in Applied Mechanics    */
/*    and Engineering 152(1-2):85-102, 22 January 1998.                      */
/*                                                                           */
/*  Triangle's Delaunay refinement algorithm for quality mesh generation is  */
/*    a hybrid of one due to Jim Ruppert, "A Delaunay Refinement Algorithm   */
/*    for Quality 2-Dimensional Mesh Generation," Journal of Algorithms      */
/*    18(3):548-585, May 1995 [*], and one due to L. Paul Chew, "Guaranteed- */
/*    Quality Mesh Generation for Curved Surfaces," Proceedings of the Ninth */
/*    Annual Symposium on Computational Geometry (San Diego, California),    */
/*    pages 274-280, Association for Computing Machinery, May 1993,          */
/*    http://portal.acm.org/citation.cfm?id=161150 .                         */
/*                                                                           */
/*  The Delaunay refinement algorithm has been modified so that it meshes    */
/*    domains with small input angles well, as described in Gary L. Miller,  */
/*    Steven E. Pav, and Noel J. Walkington, "When and Why Ruppert's         */
/*    Algorithm Works," Twelfth International Meshing Roundtable, pages      */
/*    91-102, Sandia National Laboratories, September 2003.  [*]             */
/*                                                                           */
/*  My implementation of the divide-and-conquer and incremental Delaunay     */
/*    triangulation algorithms follows closely the presentation of Guibas    */
/*    and Stolfi, even though I use a triangle-based data structure instead  */
/*    of their quad-edge data structure.  (In fact, I originally implemented */
/*    Triangle using the quad-edge data structure, but the switch to a       */
/*    triangle-based data structure sped Triangle by a factor of two.)  The  */
/*    mesh manipulation primitives and the two aforementioned Delaunay       */
/*    triangulation algorithms are described by Leonidas J. Guibas and Jorge */
/*    Stolfi, "Primitives for the Manipulation of General Subdivisions and   */
/*    the Computation of Voronoi Diagrams," ACM Transactions on Graphics     */
/*    4(2):74-123, April 1985, http://portal.acm.org/citation.cfm?id=282923 .*/
/*                                                                           */
/*  Their O(n log n) divide-and-conquer algorithm is adapted from Der-Tsai   */
/*    Lee and Bruce J. Schachter, "Two Algorithms for Constructing the       */
/*    Delaunay Triangulation," International Journal of Computer and         */
/*    Information Science 9(3):219-242, 1980.  Triangle's improvement of the */
/*    divide-and-conquer algorithm by alternating between vertical and       */
/*    horizontal cuts was introduced by Rex A. Dwyer, "A Faster Divide-and-  */
/*    Conquer Algorithm for Constructing Delaunay Triangulations,"           */
/*    Algorithmica 2(2):137-151, 1987.                                       */
/*                                                                           */
/*  The incremental insertion algorithm was first proposed by C. L. Lawson,  */
/*    "Software for C1 Surface Interpolation," in Mathematical Software III, */
/*    John R. Rice, editor, Academic Press, New York, pp. 161-194, 1977.     */
/*    For point location, I use the algorithm of Ernst P. Mucke, Isaac       */
/*    Saias, and Binhai Zhu, "Fast Randomized Point Location Without         */
/*    Preprocessing in Two- and Three-Dimensional Delaunay Triangulations,"  */
/*    Proceedings of the Twelfth Annual Symposium on Computational Geometry, */
/*    ACM, May 1996.  [*]  If I were to randomize the order of vertex        */
/*    insertion (I currently don't bother), their result combined with the   */
/*    result of Kenneth L. Clarkson and Peter W. Shor, "Applications of      */
/*    Random Sampling in Computational Geometry II," Discrete &              */
/*    Computational Geometry 4(1):387-421, 1989, would yield an expected     */
/*    O(n^{4/3}) bound on running time.                                      */
/*                                                                           */
/*  The O(n log n) sweepline Delaunay triangulation algorithm is taken from  */
/*    Steven Fortune, "A Sweepline Algorithm for Voronoi Diagrams",          */
/*    Algorithmica 2(2):153-174, 1987.  A random sample of edges on the      */
/*    boundary of the triangulation are maintained in a splay tree for the   */
/*    purpose of point location.  Splay trees are described by Daniel        */
/*    Dominic Sleator and Robert Endre Tarjan, "Self-Adjusting Binary Search */
/*    Trees," Journal of the ACM 32(3):652-686, July 1985,                   */
/*    http://portal.acm.org/citation.cfm?id=3835 .                           */
/*                                                                           */
/*  The algorithms for exact computation of the signs of determinants are    */
/*    described in Jonathan Richard Shewchuk, "Adaptive Precision Floating-  */
/*    Point Arithmetic and Fast Robust Geometric Predicates," Discrete &     */
/*    Computational Geometry 18(3):305-363, October 1997.  (Also available   */
/*    as Technical Report CMU-CS-96-140, School of Computer Science,         */
/*    Carnegie Mellon University, Pittsburgh, Pennsylvania, May 1996.)  [*]  */
/*    An abbreviated version appears as Jonathan Richard Shewchuk, "Robust   */
/*    Adaptive Floating-Point Geometric Predicates," Proceedings of the      */
/*    Twelfth Annual Symposium on Computational Geometry, ACM, May 1996. [*] */
/*    Many of the ideas for my exact arithmetic routines originate with      */
/*    Douglas M. Priest, "Algorithms for Arbitrary Precision Floating Point  */
/*    Arithmetic," Tenth Symposium on Computer Arithmetic, pp. 132-143, IEEE */
/*    Computer Society Press, 1991.  [*]  Many of the ideas for the correct  */
/*    evaluation of the signs of determinants are taken from Steven Fortune  */
/*    and Christopher J. Van Wyk, "Efficient Exact Arithmetic for Computa-   */
/*    tional Geometry," Proceedings of the Ninth Annual Symposium on         */
/*    Computational Geometry, ACM, pp. 163-172, May 1993, and from Steven    */
/*    Fortune, "Numerical Stability of Algorithms for 2D Delaunay Triangu-   */
/*    lations," International Journal of Computational Geometry & Applica-   */
/*    tions 5(1-2):193-213, March-June 1995.                                 */
/*                                                                           */
/*  The method of inserting new vertices off-center (not precisely at the    */
/*    circumcenter of every poor-quality triangle) is from Alper Ungor,      */
/*    "Off-centers:  A New Type of Steiner Points for Computing Size-Optimal */
/*    Quality-Guaranteed Delaunay Triangulations," Proceedings of LATIN      */
/*    2004 (Buenos Aires, Argentina), April 2004.                            */
/*                                                                           */
/*  For definitions of and results involving Delaunay triangulations,        */
/*    constrained and conforming versions thereof, and other aspects of      */
/*    triangular mesh generation, see the excellent survey by Marshall Bern  */
/*    and David Eppstein, "Mesh Generation and Optimal Triangulation," in    */
/*    Computing and Euclidean Geometry, Ding-Zhu Du and Frank Hwang,         */
/*    editors, World Scientific, Singapore, pp. 23-90, 1992.  [*]            */
/*                                                                           */
/*  The time for incrementally adding PSLG (planar straight line graph)      */
/*    segments to create a constrained Delaunay triangulation is probably    */
/*    O(t^2) per segment in the worst case and O(t) per segment in the       */
/*    common case, where t is the number of triangles that intersect the     */
/*    segment before it is inserted.  This doesn't count point location,     */
/*    which can be much more expensive.  I could improve this to O(d log d)  */
/*    time, but d is usually quite small, so it's not worth the bother.      */
/*    (This note does not apply when the -s switch is used, invoking a       */
/*    different method is used to insert segments.)                          */
/*                                                                           */
/*  The time for deleting a vertex from a Delaunay triangulation is O(d^2)   */
/*    in the worst case and O(d) in the common case, where d is the degree   */
/*    of the vertex being deleted.  I could improve this to O(d log d) time, */
/*    but d is usually quite small, so it's not worth the bother.            */
/*                                                                           */
/*  Ruppert's Delaunay refinement algorithm typically generates triangles    */
/*    at a linear rate (constant time per triangle) after the initial        */
/*    triangulation is formed.  There may be pathological cases where        */
/*    quadratic time is required, but these never arise in practice.         */
/*                                                                           */
/*  The geometric predicates (circumcenter calculations, segment             */
/*    intersection formulae, etc.) appear in my "Lecture Notes on Geometric  */
/*    Robustness" at http://www.cs.berkeley.edu/~jrs/mesh .                  */
/*                                                                           */
/*  If you make any improvements to this code, please please please let me   */
/*    know, so that I may obtain the improvements.  Even if you don't change */
/*    the code, I'd still love to hear what it's being used for.             */
/*                                                                           */
/*****************************************************************************/

#include <NekMeshUtils/Triangle/Triangle.h>

/*****************************************************************************/
/*                                                                           */
/*  Mesh manipulation primitives.  Each triangle contains three pointers to  */
/*  other triangles, with orientations.  Each pointer points not to the      */
/*  first byte of a triangle, but to one of the first three bytes of a       */
/*  triangle.  It is necessary to extract both the triangle itself and the   */
/*  orientation.  To save memory, I keep both pieces of information in one   */
/*  pointer.  To make this possible, I assume that all triangles are aligned */
/*  to four-byte boundaries.  The decode() routine below decodes a pointer,  */
/*  extracting an orientation (in the range 0 to 2) and a pointer to the     */
/*  beginning of a triangle.  The encode() routine compresses a pointer to a */
/*  triangle and an orientation into a single pointer.  My assumptions that  */
/*  triangles are four-byte-aligned and that the `unsigned long' type is     */
/*  long enough to hold a pointer are two of the few kludges in this program.*/
/*                                                                           */
/*  Subsegments are manipulated similarly.  A pointer to a subsegment        */
/*  carries both an address and an orientation in the range 0 to 1.          */
/*                                                                           */
/*  The other primitives take an oriented triangle or oriented subsegment,   */
/*  and return an oriented triangle or oriented subsegment or vertex; or     */
/*  they change the connections in the data structure.                       */
/*                                                                           */
/*  Below, triangles and subsegments are denoted by their vertices.  The     */
/*  triangle abc has origin (org) a, destination (dest) b, and apex (apex)   */
/*  c.  These vertices occur in counterclockwise order about the triangle.   */
/*  The handle abc may simultaneously denote vertex a, edge ab, and triangle */
/*  abc.                                                                     */
/*                                                                           */
/*  Similarly, the subsegment ab has origin (sorg) a and destination (sdest) */
/*  b.  If ab is thought to be directed upward (with b directly above a),    */
/*  then the handle ab is thought to grasp the right side of ab, and may     */
/*  simultaneously denote vertex a and edge ab.                              */
/*                                                                           */
/*  An asterisk (*) denotes a vertex whose identity is unknown.              */
/*                                                                           */
/*  Given this notation, a partial list of mesh manipulation primitives      */
/*  follows.                                                                 */
/*                                                                           */
/*                                                                           */
/*  For triangles:                                                           */
/*                                                                           */
/*  sym:  Find the abutting triangle; same edge.                             */
/*  sym(abc) -> ba*                                                          */
/*                                                                           */
/*  lnext:  Find the next edge (counterclockwise) of a triangle.             */
/*  lnext(abc) -> bca                                                        */
/*                                                                           */
/*  lprev:  Find the previous edge (clockwise) of a triangle.                */
/*  lprev(abc) -> cab                                                        */
/*                                                                           */
/*  onext:  Find the next edge counterclockwise with the same origin.        */
/*  onext(abc) -> ac*                                                        */
/*                                                                           */
/*  oprev:  Find the next edge clockwise with the same origin.               */
/*  oprev(abc) -> a*b                                                        */
/*                                                                           */
/*  dnext:  Find the next edge counterclockwise with the same destination.   */
/*  dnext(abc) -> *ba                                                        */
/*                                                                           */
/*  dprev:  Find the next edge clockwise with the same destination.          */
/*  dprev(abc) -> cb*                                                        */
/*                                                                           */
/*  rnext:  Find the next edge (counterclockwise) of the adjacent triangle.  */
/*  rnext(abc) -> *a*                                                        */
/*                                                                           */
/*  rprev:  Find the previous edge (clockwise) of the adjacent triangle.     */
/*  rprev(abc) -> b**                                                        */
/*                                                                           */
/*  org:  Origin          dest:  Destination          apex:  Apex            */
/*  org(abc) -> a         dest(abc) -> b              apex(abc) -> c         */
/*                                                                           */
/*  bond:  Bond two triangles together at the resepective handles.           */
/*  bond(abc, bad)                                                           */
/*                                                                           */
/*                                                                           */
/*  For subsegments:                                                         */
/*                                                                           */
/*  ssym:  Reverse the orientation of a subsegment.                          */
/*  ssym(ab) -> ba                                                           */
/*                                                                           */
/*  spivot:  Find adjoining subsegment with the same origin.                 */
/*  spivot(ab) -> a*                                                         */
/*                                                                           */
/*  snext:  Find next subsegment in sequence.                                */
/*  snext(ab) -> b*                                                          */
/*                                                                           */
/*  sorg:  Origin                      sdest:  Destination                   */
/*  sorg(ab) -> a                      sdest(ab) -> b                        */
/*                                                                           */
/*  sbond:  Bond two subsegments together at the respective origins.         */
/*  sbond(ab, ac)                                                            */
/*                                                                           */
/*                                                                           */
/*  For interacting tetrahedra and subfacets:                                */
/*                                                                           */
/*  tspivot:  Find a subsegment abutting a triangle.                         */
/*  tspivot(abc) -> ba                                                       */
/*                                                                           */
/*  stpivot:  Find a triangle abutting a subsegment.                         */
/*  stpivot(ab) -> ba*                                                       */
/*                                                                           */
/*  tsbond:  Bond a triangle to a subsegment.                                */
/*  tsbond(abc, ba)                                                          */
/*                                                                           */
/*****************************************************************************/

namespace Nektar
{
namespace NekMeshUtils
{

/********* User-defined triangle evaluation routine begins here      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  triunsuitable()   Determine if a triangle is unsuitable, and thus must   */
/*                    be further refined.                                    */
/*                                                                           */
/*  You may write your own procedure that decides whether or not a selected  */
/*  triangle is too big (and needs to be refined).  There are two ways to do */
/*  this.                                                                    */
/*                                                                           */
/*  (1)  Modify the procedure `triunsuitable' below, then recompile          */
/*  Triangle.                                                                */
/*                                                                           */
/*  (2)  Define the symbol EXTERNAL_TEST (either by adding the definition    */
/*  to this file, or by using the appropriate compiler switch).  This way,   */
/*  you can compile triangle.c separately from your test.  Write your own    */
/*  `triunsuitable' procedure in a separate C file (using the same prototype */
/*  as below).  Compile it and link the object code with triangle.o.         */
/*                                                                           */
/*  This procedure returns 1 if the triangle is too large and should be      */
/*  refined; 0 otherwise.                                                    */
/*                                                                           */
/*****************************************************************************/

#ifdef EXTERNAL_TEST

int triunsuitable();

#else /* not EXTERNAL_TEST */

int DelaunayTriangle::triunsuitable(vertex triorg, vertex tridest, vertex triapex, double area)
{
    double dxoa, dxda, dxod;
    double dyoa, dyda, dyod;
    double oalen, dalen, odlen;
    double maxlen;

    dxoa = triorg[0] - triapex[0];
    dyoa = triorg[1] - triapex[1];
    dxda = tridest[0] - triapex[0];
    dyda = tridest[1] - triapex[1];
    dxod = triorg[0] - tridest[0];
    dyod = triorg[1] - tridest[1];
    /* Find the squares of the lengths of the triangle's three edges. */
    oalen = dxoa * dxoa + dyoa * dyoa;
    dalen = dxda * dxda + dyda * dyda;
    odlen = dxod * dxod + dyod * dyod;
    /* Find the square of the length of the longest edge. */
    maxlen = (dalen > oalen) ? dalen : oalen;
    maxlen = (odlen > maxlen) ? odlen : maxlen;

    if (maxlen > 0.05 * (triorg[0] * triorg[0] + triorg[1] * triorg[1]) + 0.02)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

#endif /* not EXTERNAL_TEST */

/**                                                                         **/
/**                                                                         **/
/********* User-defined triangle evaluation routine ends here        *********/

/********* Memory allocation and program exit wrappers begin here    *********/
/**                                                                         **/
/**                                                                         **/

void DelaunayTriangle::triexit(int status)
{
    exit(status);
}

void *DelaunayTriangle::trimalloc(int size)
{
    void *memptr;

    memptr = (void *)malloc((unsigned int)size);
    if (memptr == (void *)NULL)
    {
        printf("Error:  Out of memory.\n");
        triexit(1);
    }
    return (memptr);
}

void DelaunayTriangle::trifree(void *memptr)
{
    free(memptr);
}

/**                                                                         **/
/**                                                                         **/
/********* Memory allocation and program exit wrappers end here      *********/

/********* User interaction routines begin here                      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  internalerror()   Ask the user to send me the defective product.  Exit.  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::internalerror()
{
    printf("  Please report this bug to jrs@cs.berkeley.edu\n");
    printf("  Include the message above, your input data set, and the exact\n");
    printf("    command line you used to run Triangle.\n");
    triexit(1);
}

/*****************************************************************************/
/*                                                                           */
/*  parsecommandline()   Read the command line, identify switches, and set   */
/*                       up options and file names.                          */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::parsecommandline(int argc, char **argv, struct behavior *b)
{
    int i, j, k;
    char workstring[2048];

    b->poly = b->quality = 0;
    b->usertest = 0;
    b->weighted = b->jettison = 0;
    b->nobisect                   = 0;
    b->minangle                   = 0.0;

    for (i = 0; i < argc; i++)
    {
        for (j = 0; argv[i][j] != '\0'; j++)
        {
            if (argv[i][j] == 'p')
            {
                b->poly = 1;
            }
            if (argv[i][j] == 'q')
            {
                b->quality = 1;
                if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                    (argv[i][j + 1] == '.'))
                {
                    k = 0;
                    while (
                        ((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                        (argv[i][j + 1] == '.'))
                    {
                        j++;
                        workstring[k] = argv[i][j];
                        k++;
                    }
                    workstring[k] = '\0';
                    b->minangle   = (double)strtod(workstring, (char **)NULL);
                }
                else
                {
                    b->minangle = 20.0;
                }
            }
            if (argv[i][j] == 'u')
            {
                b->quality  = 1;
                b->usertest = 1;
            }
            if (argv[i][j] == 'w')
            {
                b->weighted = 1;
            }
            if (argv[i][j] == 'W')
            {
                b->weighted = 2;
            }
            if (argv[i][j] == 'j')
            {
                b->jettison = 1;
            }
            if (argv[i][j] == 'Y')
            {
                b->nobisect++;
            }
        }
    }

    b->usesegments = b->poly || b->quality;
    b->goodangle   = cos(b->minangle * PI / 180.0);
    if (b->goodangle == 1.0)
    {
        b->offconstant = 0.0;
    }
    else
    {
        b->offconstant =
            0.475 * sqrt((1.0 + b->goodangle) / (1.0 - b->goodangle));
    }
    b->goodangle *= b->goodangle;

    /* Regular/weighted triangulations are incompatible with PSLGs */
    /*   and meshing.                                              */
    if (b->weighted && (b->poly || b->quality))
    {
        b->weighted = 0;
    }
}

/********* Memory management routines begin here                     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  poolzero()   Set all of a pool's fields to zero.                         */
/*                                                                           */
/*  This procedure should never be called on a pool that has any memory      */
/*  allocated to it, as that memory would leak.                              */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::poolzero(struct memorypool *pool)
{
    pool->firstblock       = (void **)NULL;
    pool->nowblock         = (void **)NULL;
    pool->nextitem         = (void *)NULL;
    pool->deaditemstack    = (void *)NULL;
    pool->pathblock        = (void **)NULL;
    pool->pathitem         = (void *)NULL;
    pool->alignbytes       = 0;
    pool->itembytes        = 0;
    pool->itemsperblock    = 0;
    pool->itemsfirstblock  = 0;
    pool->items            = 0;
    pool->maxitems         = 0;
    pool->unallocateditems = 0;
    pool->pathitemsleft    = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  poolrestart()   Deallocate all items in a pool.                          */
/*                                                                           */
/*  The pool is returned to its starting state, except that no memory is     */
/*  freed to the operating system.  Rather, the previously allocated blocks  */
/*  are ready to be reused.                                                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::poolrestart(struct memorypool *pool)
{
    unsigned long alignptr;

    pool->items    = 0;
    pool->maxitems = 0;

    /* Set the currently active block. */
    pool->nowblock = pool->firstblock;
    /* Find the first item in the pool.  Increment by the size of (void *). */
    alignptr = (unsigned long)(pool->nowblock + 1);
    /* Align the item on an `alignbytes'-byte boundary. */
    pool->nextitem = (void *)(alignptr + (unsigned long)pool->alignbytes -
                              (alignptr % (unsigned long)pool->alignbytes));
    /* There are lots of unallocated items left in this block. */
    pool->unallocateditems = pool->itemsfirstblock;
    /* The stack of deallocated items is empty. */
    pool->deaditemstack = (void *)NULL;
}

/*****************************************************************************/
/*                                                                           */
/*  poolinit()   Initialize a pool of memory for allocation of items.        */
/*                                                                           */
/*  This routine initializes the machinery for allocating items.  A `pool'   */
/*  is created whose records have size at least `bytecount'.  Items will be  */
/*  allocated in `itemcount'-item blocks.  Each item is assumed to be a      */
/*  collection of words, and either pointers or floating-point values are    */
/*  assumed to be the "primary" word type.  (The "primary" word type is used */
/*  to determine alignment of items.)  If `alignment' isn't zero, all items  */
/*  will be `alignment'-byte aligned in memory.  `alignment' must be either  */
/*  a multiple or a factor of the primary word size; powers of two are safe. */
/*  `alignment' is normally used to create a few unused bits at the bottom   */
/*  of each item's pointer, in which information may be stored.              */
/*                                                                           */
/*  Don't change this routine unless you understand it.                      */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::poolinit(struct memorypool *pool,
              int bytecount,
              int itemcount,
              int firstitemcount,
              int alignment)
{
    /* Find the proper alignment, which must be at least as large as:   */
    /*   - The parameter `alignment'.                                   */
    /*   - sizeof(void *), so the stack of dead items can be maintained */
    /*       without unaligned accesses.                                */
    if (alignment > sizeof(void *))
    {
        pool->alignbytes = alignment;
    }
    else
    {
        pool->alignbytes = sizeof(void *);
    }
    pool->itembytes =
        ((bytecount - 1) / pool->alignbytes + 1) * pool->alignbytes;
    pool->itemsperblock = itemcount;
    if (firstitemcount == 0)
    {
        pool->itemsfirstblock = itemcount;
    }
    else
    {
        pool->itemsfirstblock = firstitemcount;
    }

    /* Allocate a block of items.  Space for `itemsfirstblock' items and one  */
    /*   pointer (to point to the next block) are allocated, as well as space */
    /*   to ensure alignment of the items.                                    */
    pool->firstblock =
        (void **)trimalloc(pool->itemsfirstblock * pool->itembytes +
                           (int)sizeof(void *) + pool->alignbytes);
    /* Set the next block pointer to NULL. */
    *(pool->firstblock) = (void *)NULL;
    poolrestart(pool);
}

/*****************************************************************************/
/*                                                                           */
/*  pooldeinit()   Free to the operating system all memory taken by a pool.  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::pooldeinit(struct memorypool *pool)
{
    while (pool->firstblock != (void **)NULL)
    {
        pool->nowblock = (void **)*(pool->firstblock);
        trifree((void *)pool->firstblock);
        pool->firstblock = pool->nowblock;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  poolalloc()   Allocate space for an item.                                */
/*                                                                           */
/*****************************************************************************/

void *DelaunayTriangle::poolalloc(struct memorypool *pool)
{
    void *newitem;
    void **newblock;
    unsigned long alignptr;

    /* First check the linked list of dead items.  If the list is not   */
    /*   empty, allocate an item from the list rather than a fresh one. */
    if (pool->deaditemstack != (void *)NULL)
    {
        newitem = pool->deaditemstack; /* Take first item in list. */
        pool->deaditemstack = *(void **)pool->deaditemstack;
    }
    else
    {
        /* Check if there are any free items left in the current block. */
        if (pool->unallocateditems == 0)
        {
            /* Check if another block must be allocated. */
            if (*(pool->nowblock) == (void *)NULL)
            {
                /* Allocate a new block of items, pointed to by the previous
                 * block. */
                newblock =
                    (void **)trimalloc(pool->itemsperblock * pool->itembytes +
                                       (int)sizeof(void *) + pool->alignbytes);
                *(pool->nowblock) = (void *)newblock;
                /* The next block pointer is NULL. */
                *newblock = (void *)NULL;
            }

            /* Move to the new block. */
            pool->nowblock = (void **)*(pool->nowblock);
            /* Find the first item in the block.    */
            /*   Increment by the size of (void *). */
            alignptr = (unsigned long)(pool->nowblock + 1);
            /* Align the item on an `alignbytes'-byte boundary. */
            pool->nextitem =
                (void *)(alignptr + (unsigned long)pool->alignbytes -
                         (alignptr % (unsigned long)pool->alignbytes));
            /* There are lots of unallocated items left in this block. */
            pool->unallocateditems = pool->itemsperblock;
        }

        /* Allocate a new item. */
        newitem = pool->nextitem;
        /* Advance `nextitem' pointer to next free item in block. */
        pool->nextitem = (void *)((char *)pool->nextitem + pool->itembytes);
        pool->unallocateditems--;
        pool->maxitems++;
    }
    pool->items++;
    return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  pooldealloc()   Deallocate space for an item.                            */
/*                                                                           */
/*  The deallocated space is stored in a queue for later reuse.              */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::pooldealloc(struct memorypool *pool, void *dyingitem)
{
    /* Push freshly killed item onto stack. */
    *((void **)dyingitem) = pool->deaditemstack;
    pool->deaditemstack   = dyingitem;
    pool->items--;
}

/*****************************************************************************/
/*                                                                           */
/*  traversalinit()   Prepare to traverse the entire list of items.          */
/*                                                                           */
/*  This routine is used in conjunction with traverse().                     */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::traversalinit(struct memorypool *pool)
{
    unsigned long alignptr;

    /* Begin the traversal in the first block. */
    pool->pathblock = pool->firstblock;
    /* Find the first item in the block.  Increment by the size of (void *). */
    alignptr = (unsigned long)(pool->pathblock + 1);
    /* Align with item on an `alignbytes'-byte boundary. */
    pool->pathitem = (void *)(alignptr + (unsigned long)pool->alignbytes -
                              (alignptr % (unsigned long)pool->alignbytes));
    /* Set the number of items left in the current block. */
    pool->pathitemsleft = pool->itemsfirstblock;
}

/*****************************************************************************/
/*                                                                           */
/*  traverse()   Find the next item in the list.                             */
/*                                                                           */
/*  This routine is used in conjunction with traversalinit().  Be forewarned */
/*  that this routine successively returns all items in the list, including  */
/*  deallocated ones on the deaditemqueue.  It's up to you to figure out     */
/*  which ones are actually dead.  Why?  I don't want to allocate extra      */
/*  space just to demarcate dead items.  It can usually be done more         */
/*  space-efficiently by a routine that knows something about the structure  */
/*  of the item.                                                             */
/*                                                                           */
/*****************************************************************************/

void *DelaunayTriangle::traverse(struct memorypool *pool)
{
    void *newitem;
    unsigned long alignptr;

    /* Stop upon exhausting the list of items. */
    if (pool->pathitem == pool->nextitem)
    {
        return (void *)NULL;
    }

    /* Check whether any untraversed items remain in the current block. */
    if (pool->pathitemsleft == 0)
    {
        /* Find the next block. */
        pool->pathblock = (void **)*(pool->pathblock);
        /* Find the first item in the block.  Increment by the size of (void *).
         */
        alignptr = (unsigned long)(pool->pathblock + 1);
        /* Align with item on an `alignbytes'-byte boundary. */
        pool->pathitem = (void *)(alignptr + (unsigned long)pool->alignbytes -
                                  (alignptr % (unsigned long)pool->alignbytes));
        /* Set the number of items left in the current block. */
        pool->pathitemsleft = pool->itemsperblock;
    }

    newitem = pool->pathitem;
    /* Find the next item in the block. */
    pool->pathitem = (void *)((char *)pool->pathitem + pool->itembytes);
    pool->pathitemsleft--;
    return newitem;
}

/*****************************************************************************/
/*                                                                           */
/*  dummyinit()   Initialize the triangle that fills "outer space" and the   */
/*                omnipresent subsegment.                                    */
/*                                                                           */
/*  The triangle that fills "outer space," called `dummytri', is pointed to  */
/*  by every triangle and subsegment on a boundary (be it outer or inner) of */
/*  the triangulation.  Also, `dummytri' points to one of the triangles on   */
/*  the convex hull (until the holes and concavities are carved), making it  */
/*  possible to find a starting triangle for point location.                 */
/*                                                                           */
/*  The omnipresent subsegment, `dummysub', is pointed to by every triangle  */
/*  or subsegment that doesn't have a full complement of double subsegments */
/*  to point to.                                                             */
/*                                                                           */
/*  `dummytri' and `dummysub' are generally required to fulfill only a few   */
/*  invariants:  their vertices must remain NULL and `dummytri' must always  */
/*  be bonded (at offset zero) to some triangle on the convex hull of the    */
/*  mesh, via a boundary edge.  Otherwise, the connections of `dummytri' and */
/*  `dummysub' may change willy-nilly.  This makes it possible to avoid      */
/*  writing a good deal of special-case code (in the edge flip, for example) */
/*  for dealing with the boundary of the mesh, places where no subsegment is */
/*  present, and so forth.  Other entities are frequently bonded to          */
/*  `dummytri' and `dummysub' as if they were double mesh entities, with no */
/*  harm done.                                                               */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::dummyinit(struct mesh *m,
               struct behavior *b,
               int trianglebytes,
               int subsegbytes)
{
    unsigned long alignptr;

    /* Set up `dummytri', the `triangle' that occupies "outer space." */
    m->dummytribase =
        (triangle *)trimalloc(trianglebytes + m->triangles.alignbytes);
    /* Align `dummytri' on a `triangles.alignbytes'-byte boundary. */
    alignptr = (unsigned long)m->dummytribase;
    m->dummytri =
        (triangle *)(alignptr + (unsigned long)m->triangles.alignbytes -
                     (alignptr % (unsigned long)m->triangles.alignbytes));
    /* Initialize the three adjoining triangles to be "outer space."  These  */
    /*   will eventually be changed by various bonding operations, but their */
    /*   values don't really matter, as long as they can legally be          */
    /*   dereferenced.                                                       */
    m->dummytri[0] = (triangle)m->dummytri;
    m->dummytri[1] = (triangle)m->dummytri;
    m->dummytri[2] = (triangle)m->dummytri;
    /* Three NULL vertices. */
    m->dummytri[3] = (triangle)NULL;
    m->dummytri[4] = (triangle)NULL;
    m->dummytri[5] = (triangle)NULL;

    if (b->usesegments)
    {
        /* Set up `dummysub', the omnipresent subsegment pointed to by any */
        /*   triangle side or subsegment end that isn't attached to a double */
        /*   subsegment.                                                   */
        m->dummysubbase =
            (subseg *)trimalloc(subsegbytes + m->subsegs.alignbytes);
        /* Align `dummysub' on a `subsegs.alignbytes'-byte boundary. */
        alignptr = (unsigned long)m->dummysubbase;
        m->dummysub =
            (subseg *)(alignptr + (unsigned long)m->subsegs.alignbytes -
                       (alignptr % (unsigned long)m->subsegs.alignbytes));
        /* Initialize the two adjoining subsegments to be the omnipresent */
        /*   subsegment.  These will eventually be changed by various bonding */
        /*   operations, but their values don't really matter, as long as they
         */
        /*   can legally be dereferenced. */
        m->dummysub[0] = (subseg)m->dummysub;
        m->dummysub[1] = (subseg)m->dummysub;
        /* Four NULL vertices. */
        m->dummysub[2] = (subseg)NULL;
        m->dummysub[3] = (subseg)NULL;
        m->dummysub[4] = (subseg)NULL;
        m->dummysub[5] = (subseg)NULL;
        /* Initialize the two adjoining triangles to be "outer space." */
        m->dummysub[6] = (subseg)m->dummytri;
        m->dummysub[7] = (subseg)m->dummytri;
        /* Set the boundary marker to zero. */
        *(int *)(m->dummysub + 8) = 0;

        /* Initialize the three adjoining subsegments of `dummytri' to be */
        /*   the omnipresent subsegment.                                  */
        m->dummytri[6] = (triangle)m->dummysub;
        m->dummytri[7] = (triangle)m->dummysub;
        m->dummytri[8] = (triangle)m->dummysub;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  initializevertexpool()   Calculate the size of the vertex data structure */
/*                           and initialize its memory pool.                 */
/*                                                                           */
/*  This routine also computes the `vertexmarkindex' and `vertex2triindex'   */
/*  indices used to find values within each vertex.                          */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::initializevertexpool(struct mesh *m, struct behavior *b)
{
    int vertexsize;

    /* The index within each vertex at which the boundary marker is found,    */
    /*   followed by the vertex type.  Ensure the vertex marker is aligned to */
    /*   a sizeof(int)-byte address.                                          */
    m->vertexmarkindex =
        ((m->mesh_dim + m->nextras) * sizeof(double) + sizeof(int) - 1) /
        sizeof(int);
    vertexsize = (m->vertexmarkindex + 2) * sizeof(int);
    if (b->poly)
    {
        /* The index within each vertex at which a triangle pointer is found. */
        /*   Ensure the pointer is aligned to a sizeof(triangle)-byte address.
         */
        m->vertex2triindex =
            (vertexsize + sizeof(triangle) - 1) / sizeof(triangle);
        vertexsize = (m->vertex2triindex + 1) * sizeof(triangle);
    }

    /* Initialize the pool of vertices. */
    poolinit(&m->vertices,
             vertexsize,
             VERTEXPERBLOCK,
             m->invertices > VERTEXPERBLOCK ? m->invertices : VERTEXPERBLOCK,
             sizeof(double));
}

/*****************************************************************************/
/*                                                                           */
/*  initializetrisubpools()   Calculate the sizes of the triangle and        */
/*                            subsegment data structures and initialize      */
/*                            their memory pools.                            */
/*                                                                           */
/*  This routine also computes the `highorderindex', `elemattribindex', and  */
/*  `areaboundindex' indices used to find values within each triangle.       */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::initializetrisubpools(struct mesh *m, struct behavior *b)
{
    int trisize;

    /* The index within each triangle at which the extra nodes (above three)  */
    /*   associated with high order elements are found.  There are three      */
    /*   pointers to other triangles, three pointers to corners, and possibly */
    /*   three pointers to subsegments before the extra nodes.                */
    m->highorderindex = 6 + (b->usesegments * 3);
    /* The number of bytes occupied by a triangle. */
    trisize = (3 + (m->highorderindex - 3)) *
              sizeof(triangle);
    /* The index within each triangle at which its attributes are found, */
    /*   where the index is measured in doubles.                           */
    m->elemattribindex = (trisize + sizeof(double) - 1) / sizeof(double);
    /* The index within each triangle at which the maximum area constraint  */
    /*   is found, where the index is measured in doubles.  Note that if the  */
    /*   `regionattrib' flag is set, an additional attribute will be added. */
    m->areaboundindex = m->elemattribindex + m->eextras;
    /* If triangle attributes or an area bound are needed, increase the number
     */
    /*   of bytes occupied by a triangle. */
    trisize = m->areaboundindex * sizeof(double);

    /* Having determined the memory size of a triangle, initialize the pool. */
    poolinit(&m->triangles,
             trisize,
             TRIPERBLOCK,
             (2 * m->invertices - 2) > TRIPERBLOCK ? (2 * m->invertices - 2)
                                                   : TRIPERBLOCK,
             4);

    if (b->usesegments)
    {
        /* Initialize the pool of subsegments.  Take into account all eight */
        /*   pointers and one boundary marker.                              */
        poolinit(&m->subsegs,
                 8 * sizeof(triangle) + sizeof(int),
                 SUBSEGPERBLOCK,
                 SUBSEGPERBLOCK,
                 4);

        /* Initialize the "outer space" triangle and omnipresent subsegment. */
        dummyinit(m, b, m->triangles.itembytes, m->subsegs.itembytes);
    }
    else
    {
        /* Initialize the "outer space" triangle. */
        dummyinit(m, b, m->triangles.itembytes, 0);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  triangledealloc()   Deallocate space for a triangle, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::triangledealloc(struct mesh *m, triangle *dyingtriangle)
{
    /* Mark the triangle as dead.  This makes it possible to detect dead */
    /*   triangles when traversing the list of all triangles.            */
    killtri(dyingtriangle);
    pooldealloc(&m->triangles, (void *)dyingtriangle);
}

/*****************************************************************************/
/*                                                                           */
/*  triangletraverse()   Traverse the triangles, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

triangle *DelaunayTriangle::triangletraverse(struct mesh *m)
{
    triangle *newtriangle;

    do
    {
        newtriangle = (triangle *)traverse(&m->triangles);
        if (newtriangle == (triangle *)NULL)
        {
            return (triangle *)NULL;
        }
    } while (deadtri(newtriangle)); /* Skip dead ones. */
    return newtriangle;
}

/*****************************************************************************/
/*                                                                           */
/*  subsegdealloc()   Deallocate space for a subsegment, marking it dead.    */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::subsegdealloc(struct mesh *m, subseg *dyingsubseg)
{
    /* Mark the subsegment as dead.  This makes it possible to detect dead */
    /*   subsegments when traversing the list of all subsegments.          */
    killsubseg(dyingsubseg);
    pooldealloc(&m->subsegs, (void *)dyingsubseg);
}

/*****************************************************************************/
/*                                                                           */
/*  subsegtraverse()   Traverse the subsegments, skipping dead ones.         */
/*                                                                           */
/*****************************************************************************/

subseg *DelaunayTriangle::subsegtraverse(struct mesh *m)
{
    subseg *newsubseg;

    do
    {
        newsubseg = (subseg *)traverse(&m->subsegs);
        if (newsubseg == (subseg *)NULL)
        {
            return (subseg *)NULL;
        }
    } while (deadsubseg(newsubseg)); /* Skip dead ones. */
    return newsubseg;
}

/*****************************************************************************/
/*                                                                           */
/*  vertexdealloc()   Deallocate space for a vertex, marking it dead.        */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::vertexdealloc(struct mesh *m, vertex dyingvertex)
{
    /* Mark the vertex as dead.  This makes it possible to detect dead */
    /*   vertices when traversing the list of all vertices.            */
    setvertextype(dyingvertex, DEADVERTEX);
    pooldealloc(&m->vertices, (void *)dyingvertex);
}

/*****************************************************************************/
/*                                                                           */
/*  vertextraverse()   Traverse the vertices, skipping dead ones.            */
/*                                                                           */
/*****************************************************************************/

vertex DelaunayTriangle::vertextraverse(struct mesh *m)
{
    vertex newvertex;

    do
    {
        newvertex = (vertex)traverse(&m->vertices);
        if (newvertex == (vertex)NULL)
        {
            return (vertex)NULL;
        }
    } while (vertextype(newvertex) == DEADVERTEX); /* Skip dead ones. */
    return newvertex;
}

/*****************************************************************************/
/*                                                                           */
/*  badsubsegdealloc()   Deallocate space for a bad subsegment, marking it   */
/*                       dead.                                               */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::badsubsegdealloc(struct mesh *m, struct badsubseg *dyingseg)
{
    /* Set subsegment's origin to NULL.  This makes it possible to detect dead
     */
    /*   badsubsegs when traversing the list of all badsubsegs             . */
    dyingseg->subsegorg = (vertex)NULL;
    pooldealloc(&m->badsubsegs, (void *)dyingseg);
}

/*****************************************************************************/
/*                                                                           */
/*  badsubsegtraverse()   Traverse the bad subsegments, skipping dead ones.  */
/*                                                                           */
/*****************************************************************************/

struct badsubseg *DelaunayTriangle::badsubsegtraverse(struct mesh *m)
{
    struct badsubseg *newseg;

    do
    {
        newseg = (struct badsubseg *)traverse(&m->badsubsegs);
        if (newseg == (struct badsubseg *)NULL)
        {
            return (struct badsubseg *)NULL;
        }
    } while (newseg->subsegorg == (vertex)NULL); /* Skip dead ones. */
    return newseg;
}

/*****************************************************************************/
/*                                                                           */
/*  getvertex()   Get a specific vertex, by number, from the list.           */
/*                                                                           */
/*  The first vertex is number 'firstnumber'.                                */
/*                                                                           */
/*  Note that this takes O(n) time (with a small constant, if VERTEXPERBLOCK */
/*  is large).  I don't care to take the trouble to make it work in constant */
/*  time.                                                                    */
/*                                                                           */
/*****************************************************************************/

vertex DelaunayTriangle::getvertex(struct mesh *m, struct behavior *b, int number)
{
    void **getblock;
    char *foundvertex;
    unsigned long alignptr;
    int current;

    getblock = m->vertices.firstblock;
    current  = 0;

    /* Find the right block. */
    if (current + m->vertices.itemsfirstblock <= number)
    {
        getblock = (void **)*getblock;
        current += m->vertices.itemsfirstblock;
        while (current + m->vertices.itemsperblock <= number)
        {
            getblock = (void **)*getblock;
            current += m->vertices.itemsperblock;
        }
    }

    /* Now find the right vertex. */
    alignptr    = (unsigned long)(getblock + 1);
    foundvertex = (char *)(alignptr + (unsigned long)m->vertices.alignbytes -
                           (alignptr % (unsigned long)m->vertices.alignbytes));
    return (vertex)(foundvertex + m->vertices.itembytes * (number - current));
}

/*****************************************************************************/
/*                                                                           */
/*  triangledeinit()   Free all remaining allocated memory.                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::triangledeinit(struct mesh *m, struct behavior *b)
{
    pooldeinit(&m->triangles);
    trifree((void *)m->dummytribase);
    if (b->usesegments)
    {
        pooldeinit(&m->subsegs);
        trifree((void *)m->dummysubbase);
    }
    pooldeinit(&m->vertices);
    if (b->quality)
    {
        pooldeinit(&m->badsubsegs);
        if ((b->minangle > 0.0) || b->usertest)
        {
            pooldeinit(&m->badtriangles);
            pooldeinit(&m->flipstackers);
        }
    }
}

/**                                                                         **/
/**                                                                         **/
/********* Memory management routines end here                       *********/

/********* Constructors begin here                                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  maketriangle()   Create a new triangle with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::maketriangle(struct mesh *m, struct behavior *b, struct otri *newotri)
{
    int i;

    newotri->tri = (triangle *)poolalloc(&m->triangles);
    /* Initialize the three adjoining triangles to be "outer space". */
    newotri->tri[0] = (triangle)m->dummytri;
    newotri->tri[1] = (triangle)m->dummytri;
    newotri->tri[2] = (triangle)m->dummytri;
    /* Three NULL vertices. */
    newotri->tri[3] = (triangle)NULL;
    newotri->tri[4] = (triangle)NULL;
    newotri->tri[5] = (triangle)NULL;
    if (b->usesegments)
    {
        /* Initialize the three adjoining subsegments to be the omnipresent */
        /*   subsegment.                                                    */
        newotri->tri[6] = (triangle)m->dummysub;
        newotri->tri[7] = (triangle)m->dummysub;
        newotri->tri[8] = (triangle)m->dummysub;
    }
    for (i = 0; i < m->eextras; i++)
    {
        setelemattribute(*newotri, i, 0.0);
    }

    newotri->orient = 0;
}

/*****************************************************************************/
/*                                                                           */
/*  makesubseg()   Create a new subsegment with orientation zero.            */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::makesubseg(struct mesh *m, struct osub *newsubseg)
{
    newsubseg->ss = (subseg *)poolalloc(&m->subsegs);
    /* Initialize the two adjoining subsegments to be the omnipresent */
    /*   subsegment.                                                  */
    newsubseg->ss[0] = (subseg)m->dummysub;
    newsubseg->ss[1] = (subseg)m->dummysub;
    /* Four NULL vertices. */
    newsubseg->ss[2] = (subseg)NULL;
    newsubseg->ss[3] = (subseg)NULL;
    newsubseg->ss[4] = (subseg)NULL;
    newsubseg->ss[5] = (subseg)NULL;
    /* Initialize the two adjoining triangles to be "outer space." */
    newsubseg->ss[6] = (subseg)m->dummytri;
    newsubseg->ss[7] = (subseg)m->dummytri;
    /* Set the boundary marker to zero. */
    setmark(*newsubseg, 0);

    newsubseg->ssorient = 0;
}

/**                                                                         **/
/**                                                                         **/
/********* Constructors end here                                     *********/

/********* Geometric primitives begin here                           *********/
/**                                                                         **/
/**                                                                         **/

/* The adaptive exact arithmetic geometric predicates implemented herein are */
/*   described in detail in my paper, "Adaptive Precision Floating-Point     */
/*   Arithmetic and Fast Robust Geometric Predicates."  See the header for a */
/*   full citation.                                                          */

/* Which of the following two methods of finding the absolute values is      */
/*   fastest is compiler-dependent.  A few compilers can inline and optimize */
/*   the fabs() call; but most will incur the overhead of a function call,   */
/*   which is disastrously slow.  A faster way on IEEE machines might be to  */
/*   mask the appropriate bit, but that's difficult to do in C without       */
/*   forcing the value to be stored to memory (rather than be kept in the    */
/*   register to which the optimizer assigned it).                           */

#define Absolute(a) ((a) >= 0.0 ? (a) : -(a))
/* #define Absolute(a)  fabs(a) */

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */
/*                                                                           */
/* The operations Fast_Two_Sum(), Fast_Two_Diff(), Two_Sum(), Two_Diff(),    */
/*   Split(), and Two_Product() are all implemented as described in the      */
/*   reference.  Each of these macros requires certain variables to be       */
/*   defined in the calling routine.  The variables `bvirt', `c', `abig',    */
/*   `_i', `_j', `_k', `_l', `_m', and `_n' are declared `' because   */
/*   they store the result of an operation that may incur roundoff error.    */
/*   The input parameter `x' (or the highest numbered `x_' parameter) must   */
/*   also be declared `'.                                             */

#define Fast_Two_Sum_Tail(a, b, x, y)                                          \
    bvirt = x - a;                                                             \
    y     = b - bvirt

#define Fast_Two_Sum(a, b, x, y)                                               \
    x = (double)(a + b);                                                       \
    Fast_Two_Sum_Tail(a, b, x, y)

#define Two_Sum_Tail(a, b, x, y)                                               \
    bvirt  = (double)(x - a);                                                  \
    avirt  = x - bvirt;                                                        \
    bround = b - bvirt;                                                        \
    around = a - avirt;                                                        \
    y      = around + bround

#define Two_Sum(a, b, x, y)                                                    \
    x = (double)(a + b);                                                       \
    Two_Sum_Tail(a, b, x, y)

#define Two_Diff_Tail(a, b, x, y)                                              \
    bvirt  = (double)(a - x);                                                  \
    avirt  = x + bvirt;                                                        \
    bround = bvirt - b;                                                        \
    around = a - avirt;                                                        \
    y      = around + bround

#define Two_Diff(a, b, x, y)                                                   \
    x = (double)(a - b);                                                       \
    Two_Diff_Tail(a, b, x, y)

#define Split(a, ahi, alo)                                                     \
    c    = (double)(splitter * a);                                             \
    abig = (double)(c - a);                                                    \
    ahi  = c - abig;                                                           \
    alo  = a - ahi

#define Two_Product_Tail(a, b, x, y)                                           \
    Split(a, ahi, alo);                                                        \
    Split(b, bhi, blo);                                                        \
    err1 = x - (ahi * bhi);                                                    \
    err2 = err1 - (alo * bhi);                                                 \
    err3 = err2 - (ahi * blo);                                                 \
    y    = (alo * blo) - err3

#define Two_Product(a, b, x, y)                                                \
    x = (double)(a * b);                                                       \
    Two_Product_Tail(a, b, x, y)

/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#define Two_Product_Presplit(a, b, bhi, blo, x, y)                             \
    x = (double)(a * b);                                                       \
    Split(a, ahi, alo);                                                        \
    err1 = x - (ahi * bhi);                                                    \
    err2 = err1 - (alo * bhi);                                                 \
    err3 = err2 - (ahi * blo);                                                 \
    y    = (alo * blo) - err3

/* Square() can be done more quickly than Two_Product().                     */

#define Square_Tail(a, x, y)                                                   \
    Split(a, ahi, alo);                                                        \
    err1 = x - (ahi * ahi);                                                    \
    err3 = err1 - ((ahi + ahi) * alo);                                         \
    y    = (alo * alo) - err3

#define Square(a, x, y)                                                        \
    x = (double)(a * a);                                                       \
    Square_Tail(a, x, y)

/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#define Two_One_Sum(a1, a0, b, x2, x1, x0)                                     \
    Two_Sum(a0, b, _i, x0);                                                    \
    Two_Sum(a1, _i, x2, x1)

#define Two_One_Diff(a1, a0, b, x2, x1, x0)                                    \
    Two_Diff(a0, b, _i, x0);                                                   \
    Two_Sum(a1, _i, x2, x1)

#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0)                            \
    Two_One_Sum(a1, a0, b0, _j, _0, x0);                                       \
    Two_One_Sum(_j, _0, b1, x3, x2, x1)

#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0)                           \
    Two_One_Diff(a1, a0, b0, _j, _0, x0);                                      \
    Two_One_Diff(_j, _0, b1, x3, x2, x1)

/* Macro for multiplying a two-component expansion by a single component.    */

#define Two_One_Product(a1, a0, b, x3, x2, x1, x0)                             \
    Split(b, bhi, blo);                                                        \
    Two_Product_Presplit(a0, b, bhi, blo, _i, x0);                             \
    Two_Product_Presplit(a1, b, bhi, blo, _j, _0);                             \
    Two_Sum(_i, _0, _k, x1);                                                   \
    Fast_Two_Sum(_j, _k, x3, x2)

/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like double arithmetic. */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::exactinit()
{
    double half;
    double check, lastcheck;
    int every_other;

    every_other = 1;
    half        = 0.5;
    epsilon     = 1.0;
    splitter    = 1.0;
    check       = 1.0;
    /* Repeatedly divide `epsilon' by two until it is too small to add to */
    /*   one without causing roundoff.  (Also check if the sum is equal to */
    /*   the previous sum, for machines that round up instead of using exact */
    /*   rounding.  Not that these routines will work on such machines.) */
    do
    {
        lastcheck = check;
        epsilon *= half;
        if (every_other)
        {
            splitter *= 2.0;
        }
        every_other = !every_other;
        check       = 1.0 + epsilon;
    } while ((check != 1.0) && (check != lastcheck));
    splitter += 1.0;
    /* Error bounds for orientation and incircle tests. */
    resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
    ccwerrboundA   = (3.0 + 16.0 * epsilon) * epsilon;
    ccwerrboundB   = (2.0 + 12.0 * epsilon) * epsilon;
    ccwerrboundC   = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    iccerrboundA   = (10.0 + 96.0 * epsilon) * epsilon;
    iccerrboundB   = (4.0 + 48.0 * epsilon) * epsilon;
    iccerrboundC   = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
    o3derrboundA   = (7.0 + 56.0 * epsilon) * epsilon;
    o3derrboundB   = (3.0 + 28.0 * epsilon) * epsilon;
    o3derrboundC   = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelimTRI()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See my Robust Predicates paper for details.             */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

int DelaunayTriangle::fast_expansion_sum_zeroelimTRI(
    int elen, double *e, int flen, double *f, double *h)
{
    double Q;
     double Qnew;
     double hh;
     double bvirt;
    double avirt, bround, around;
    int eindex, findex, hindex;
    double enow, fnow;

    enow   = e[0];
    fnow   = f[0];
    eindex = findex = 0;
    if ((fnow > enow) == (fnow > -enow))
    {
        Q    = enow;
        enow = e[++eindex];
    }
    else
    {
        Q    = fnow;
        fnow = f[++findex];
    }
    hindex = 0;
    if ((eindex < elen) && (findex < flen))
    {
        if ((fnow > enow) == (fnow > -enow))
        {
            Fast_Two_Sum(enow, Q, Qnew, hh);
            enow = e[++eindex];
        }
        else
        {
            Fast_Two_Sum(fnow, Q, Qnew, hh);
            fnow = f[++findex];
        }
        Q = Qnew;
        if (hh != 0.0)
        {
            h[hindex++] = hh;
        }
        while ((eindex < elen) && (findex < flen))
        {
            if ((fnow > enow) == (fnow > -enow))
            {
                Two_Sum(Q, enow, Qnew, hh);
                enow = e[++eindex];
            }
            else
            {
                Two_Sum(Q, fnow, Qnew, hh);
                fnow = f[++findex];
            }
            Q = Qnew;
            if (hh != 0.0)
            {
                h[hindex++] = hh;
            }
        }
    }
    while (eindex < elen)
    {
        Two_Sum(Q, enow, Qnew, hh);
        enow = e[++eindex];
        Q    = Qnew;
        if (hh != 0.0)
        {
            h[hindex++] = hh;
        }
    }
    while (findex < flen)
    {
        Two_Sum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
        Q    = Qnew;
        if (hh != 0.0)
        {
            h[hindex++] = hh;
        }
    }
    if ((Q != 0.0) || (hindex == 0))
    {
        h[hindex++] = Q;
    }
    return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelimTRI()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See my Robust Predicates paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/

int DelaunayTriangle::scale_expansion_zeroelimTRI(int elen, double *e, double b, double *h)
{
     double Q, sum;
    double hh;
     double product1;
    double product0;
    int eindex, hindex;
    double enow;
     double bvirt;
    double avirt, bround, around;
     double c;
     double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;

    Split(b, bhi, blo);
    Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
    hindex = 0;
    if (hh != 0)
    {
        h[hindex++] = hh;
    }
    for (eindex = 1; eindex < elen; eindex++)
    {
        enow = e[eindex];
        Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
        Two_Sum(Q, product0, sum, hh);
        if (hh != 0)
        {
            h[hindex++] = hh;
        }
        Fast_Two_Sum(product1, sum, Q, hh);
        if (hh != 0)
        {
            h[hindex++] = hh;
        }
    }
    if ((Q != 0.0) || (hindex == 0))
    {
        h[hindex++] = Q;
    }
    return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  estimateTRI()   Produce a one-word estimateTRI of an expansion's value.        */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

double DelaunayTriangle::estimateTRI(int elen, double *e)
{
    double Q;
    int eindex;

    Q = e[0];
    for (eindex = 1; eindex < elen; eindex++)
    {
        Q += e[eindex];
    }
    return Q;
}

/*****************************************************************************/
/*                                                                           */
/*  counterclockwise()   Return a positive value if the points pa, pb, and   */
/*                       pc occur in counterclockwise order; a negative      */
/*                       value if they occur in clockwise order; and zero    */
/*                       if they are collinear.  The result is also a rough  */
/*                       approximation of twice the signed area of the       */
/*                       triangle defined by the three points.               */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are collinear or nearly so.            */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

double DelaunayTriangle::counterclockwiseadapt(vertex pa, vertex pb, vertex pc, double detsum)
{
     double acx, acy, bcx, bcy;
    double acxtail, acytail, bcxtail, bcytail;
     double detleft, detright;
    double detlefttail, detrighttail;
    double det, errbound;
    double B[4], C1[8], C2[12], D[16];
     double B3;
    int C1length, C2length, Dlength;
    double u[4];
     double u3;
     double s1, t1;
    double s0, t0;

     double bvirt;
    double avirt, bround, around;
     double c;
     double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
     double _i, _j;
    double _0;

    acx = (double)(pa[0] - pc[0]);
    bcx = (double)(pb[0] - pc[0]);
    acy = (double)(pa[1] - pc[1]);
    bcy = (double)(pb[1] - pc[1]);

    Two_Product(acx, bcy, detleft, detlefttail);
    Two_Product(acy, bcx, detright, detrighttail);

    Two_Two_Diff(
        detleft, detlefttail, detright, detrighttail, B3, B[2], B[1], B[0]);
    B[3] = B3;

    det      = estimateTRI(4, B);
    errbound = ccwerrboundB * detsum;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    Two_Diff_Tail(pa[0], pc[0], acx, acxtail);
    Two_Diff_Tail(pb[0], pc[0], bcx, bcxtail);
    Two_Diff_Tail(pa[1], pc[1], acy, acytail);
    Two_Diff_Tail(pb[1], pc[1], bcy, bcytail);

    if ((acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) &&
        (bcytail == 0.0))
    {
        return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * Absolute(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    Two_Product(acxtail, bcy, s1, s0);
    Two_Product(acytail, bcx, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3]     = u3;
    C1length = fast_expansion_sum_zeroelimTRI(4, B, 4, u, C1);

    Two_Product(acx, bcytail, s1, s0);
    Two_Product(acy, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3]     = u3;
    C2length = fast_expansion_sum_zeroelimTRI(C1length, C1, 4, u, C2);

    Two_Product(acxtail, bcytail, s1, s0);
    Two_Product(acytail, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3]    = u3;
    Dlength = fast_expansion_sum_zeroelimTRI(C2length, C2, 4, u, D);

    return (D[Dlength - 1]);
}

double DelaunayTriangle::counterclockwise(
    struct mesh *m, struct behavior *b, vertex pa, vertex pb, vertex pc)
{
    double detleft, detright, det;
    double detsum, errbound;

    m->counterclockcount++;

    detleft  = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    det      = detleft - detright;

    if (detleft > 0.0)
    {
        if (detright <= 0.0)
        {
            return det;
        }
        else
        {
            detsum = detleft + detright;
        }
    }
    else if (detleft < 0.0)
    {
        if (detright >= 0.0)
        {
            return det;
        }
        else
        {
            detsum = -detleft - detright;
        }
    }
    else
    {
        return det;
    }

    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    return counterclockwiseadapt(pa, pb, pc, detsum);
}

/*****************************************************************************/
/*                                                                           */
/*  incircle()   Return a positive value if the point pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are cocircular or nearly so.           */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

double DelaunayTriangle::incircleadaptTRI(
    vertex pa, vertex pb, vertex pc, vertex pd, double permanent)
{
     double adx, bdx, cdx, ady, bdy, cdy;
    double det, errbound;

     double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
    double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
    double bc[4], ca[4], ab[4];
     double bc3, ca3, ab3;
    double axbc[8], axxbc[16], aybc[8], ayybc[16], adet[32];
    int axbclen, axxbclen, aybclen, ayybclen, alen;
    double bxca[8], bxxca[16], byca[8], byyca[16], bdet[32];
    int bxcalen, bxxcalen, bycalen, byycalen, blen;
    double cxab[8], cxxab[16], cyab[8], cyyab[16], cdet[32];
    int cxablen, cxxablen, cyablen, cyyablen, clen;
    double abdet[64];
    int ablen;
    double fin1[1152], fin2[1152];
    double *finnow, *finother, *finswap;
    int finlength;

    double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
     double adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
    double adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
    double aa[4], bb[4], cc[4];
     double aa3, bb3, cc3;
     double ti1, tj1;
    double ti0, tj0;
    double u[4], v[4];
     double u3, v3;
    double temp8[8], temp16a[16], temp16b[16], temp16c[16];
    double temp32a[32], temp32b[32], temp48[48], temp64[64];
    int temp8len, temp16alen, temp16blen, temp16clen;
    int temp32alen, temp32blen, temp48len, temp64len;
    double axtbb[8], axtcc[8], aytbb[8], aytcc[8];
    int axtbblen, axtcclen, aytbblen, aytcclen;
    double bxtaa[8], bxtcc[8], bytaa[8], bytcc[8];
    int bxtaalen, bxtcclen, bytaalen, bytcclen;
    double cxtaa[8], cxtbb[8], cytaa[8], cytbb[8];
    int cxtaalen, cxtbblen, cytaalen, cytbblen;
    double axtbc[8], aytbc[8], bxtca[8], bytca[8], cxtab[8], cytab[8];
    int axtbclen, aytbclen, bxtcalen, bytcalen, cxtablen, cytablen;
    double axtbct[16], aytbct[16], bxtcat[16], bytcat[16], cxtabt[16],
        cytabt[16];
    int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
    double axtbctt[8], aytbctt[8], bxtcatt[8];
    double bytcatt[8], cxtabtt[8], cytabtt[8];
    int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
    double abt[8], bct[8], cat[8];
    int abtlen, bctlen, catlen;
    double abtt[4], bctt[4], catt[4];
    int abttlen, bcttlen, cattlen;
     double abtt3, bctt3, catt3;
    double negate;

     double bvirt;
    double avirt, bround, around;
     double c;
     double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
     double _i, _j;
    double _0;

    adx = (double)(pa[0] - pd[0]);
    bdx = (double)(pb[0] - pd[0]);
    cdx = (double)(pc[0] - pd[0]);
    ady = (double)(pa[1] - pd[1]);
    bdy = (double)(pb[1] - pd[1]);
    cdy = (double)(pc[1] - pd[1]);

    Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
    Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
    Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
    bc[3]    = bc3;
    axbclen  = scale_expansion_zeroelimTRI(4, bc, adx, axbc);
    axxbclen = scale_expansion_zeroelimTRI(axbclen, axbc, adx, axxbc);
    aybclen  = scale_expansion_zeroelimTRI(4, bc, ady, aybc);
    ayybclen = scale_expansion_zeroelimTRI(aybclen, aybc, ady, ayybc);
    alen = fast_expansion_sum_zeroelimTRI(axxbclen, axxbc, ayybclen, ayybc, adet);

    Two_Product(cdx, ady, cdxady1, cdxady0);
    Two_Product(adx, cdy, adxcdy1, adxcdy0);
    Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
    ca[3]    = ca3;
    bxcalen  = scale_expansion_zeroelimTRI(4, ca, bdx, bxca);
    bxxcalen = scale_expansion_zeroelimTRI(bxcalen, bxca, bdx, bxxca);
    bycalen  = scale_expansion_zeroelimTRI(4, ca, bdy, byca);
    byycalen = scale_expansion_zeroelimTRI(bycalen, byca, bdy, byyca);
    blen = fast_expansion_sum_zeroelimTRI(bxxcalen, bxxca, byycalen, byyca, bdet);

    Two_Product(adx, bdy, adxbdy1, adxbdy0);
    Two_Product(bdx, ady, bdxady1, bdxady0);
    Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
    ab[3]    = ab3;
    cxablen  = scale_expansion_zeroelimTRI(4, ab, cdx, cxab);
    cxxablen = scale_expansion_zeroelimTRI(cxablen, cxab, cdx, cxxab);
    cyablen  = scale_expansion_zeroelimTRI(4, ab, cdy, cyab);
    cyyablen = scale_expansion_zeroelimTRI(cyablen, cyab, cdy, cyyab);
    clen = fast_expansion_sum_zeroelimTRI(cxxablen, cxxab, cyyablen, cyyab, cdet);

    ablen     = fast_expansion_sum_zeroelimTRI(alen, adet, blen, bdet, abdet);
    finlength = fast_expansion_sum_zeroelimTRI(ablen, abdet, clen, cdet, fin1);

    det      = estimateTRI(finlength, fin1);
    errbound = iccerrboundB * permanent;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
    Two_Diff_Tail(pa[1], pd[1], ady, adytail);
    Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
    Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
    Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
    Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
    if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) &&
        (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0))
    {
        return det;
    }

    errbound = iccerrboundC * permanent + resulterrbound * Absolute(det);
    det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) -
                                       (bdy * cdxtail + cdx * bdytail)) +
            2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)) +
           ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) -
                                       (cdy * adxtail + adx * cdytail)) +
            2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) +
           ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) -
                                       (ady * bdxtail + bdx * adytail)) +
            2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    finnow   = fin1;
    finother = fin2;

    if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) ||
        (cdytail != 0.0))
    {
        Square(adx, adxadx1, adxadx0);
        Square(ady, adyady1, adyady0);
        Two_Two_Sum(
            adxadx1, adxadx0, adyady1, adyady0, aa3, aa[2], aa[1], aa[0]);
        aa[3] = aa3;
    }
    if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) ||
        (adytail != 0.0))
    {
        Square(bdx, bdxbdx1, bdxbdx0);
        Square(bdy, bdybdy1, bdybdy0);
        Two_Two_Sum(
            bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb3, bb[2], bb[1], bb[0]);
        bb[3] = bb3;
    }
    if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) ||
        (bdytail != 0.0))
    {
        Square(cdx, cdxcdx1, cdxcdx0);
        Square(cdy, cdycdy1, cdycdy0);
        Two_Two_Sum(
            cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc3, cc[2], cc[1], cc[0]);
        cc[3] = cc3;
    }

    if (adxtail != 0.0)
    {
        axtbclen = scale_expansion_zeroelimTRI(4, bc, adxtail, axtbc);
        temp16alen =
            scale_expansion_zeroelimTRI(axtbclen, axtbc, 2.0 * adx, temp16a);

        axtcclen   = scale_expansion_zeroelimTRI(4, cc, adxtail, axtcc);
        temp16blen = scale_expansion_zeroelimTRI(axtcclen, axtcc, bdy, temp16b);

        axtbblen   = scale_expansion_zeroelimTRI(4, bb, adxtail, axtbb);
        temp16clen = scale_expansion_zeroelimTRI(axtbblen, axtbb, -cdy, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (adytail != 0.0)
    {
        aytbclen = scale_expansion_zeroelimTRI(4, bc, adytail, aytbc);
        temp16alen =
            scale_expansion_zeroelimTRI(aytbclen, aytbc, 2.0 * ady, temp16a);

        aytbblen   = scale_expansion_zeroelimTRI(4, bb, adytail, aytbb);
        temp16blen = scale_expansion_zeroelimTRI(aytbblen, aytbb, cdx, temp16b);

        aytcclen   = scale_expansion_zeroelimTRI(4, cc, adytail, aytcc);
        temp16clen = scale_expansion_zeroelimTRI(aytcclen, aytcc, -bdx, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (bdxtail != 0.0)
    {
        bxtcalen = scale_expansion_zeroelimTRI(4, ca, bdxtail, bxtca);
        temp16alen =
            scale_expansion_zeroelimTRI(bxtcalen, bxtca, 2.0 * bdx, temp16a);

        bxtaalen   = scale_expansion_zeroelimTRI(4, aa, bdxtail, bxtaa);
        temp16blen = scale_expansion_zeroelimTRI(bxtaalen, bxtaa, cdy, temp16b);

        bxtcclen   = scale_expansion_zeroelimTRI(4, cc, bdxtail, bxtcc);
        temp16clen = scale_expansion_zeroelimTRI(bxtcclen, bxtcc, -ady, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (bdytail != 0.0)
    {
        bytcalen = scale_expansion_zeroelimTRI(4, ca, bdytail, bytca);
        temp16alen =
            scale_expansion_zeroelimTRI(bytcalen, bytca, 2.0 * bdy, temp16a);

        bytcclen   = scale_expansion_zeroelimTRI(4, cc, bdytail, bytcc);
        temp16blen = scale_expansion_zeroelimTRI(bytcclen, bytcc, adx, temp16b);

        bytaalen   = scale_expansion_zeroelimTRI(4, aa, bdytail, bytaa);
        temp16clen = scale_expansion_zeroelimTRI(bytaalen, bytaa, -cdx, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (cdxtail != 0.0)
    {
        cxtablen = scale_expansion_zeroelimTRI(4, ab, cdxtail, cxtab);
        temp16alen =
            scale_expansion_zeroelimTRI(cxtablen, cxtab, 2.0 * cdx, temp16a);

        cxtbblen   = scale_expansion_zeroelimTRI(4, bb, cdxtail, cxtbb);
        temp16blen = scale_expansion_zeroelimTRI(cxtbblen, cxtbb, ady, temp16b);

        cxtaalen   = scale_expansion_zeroelimTRI(4, aa, cdxtail, cxtaa);
        temp16clen = scale_expansion_zeroelimTRI(cxtaalen, cxtaa, -bdy, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (cdytail != 0.0)
    {
        cytablen = scale_expansion_zeroelimTRI(4, ab, cdytail, cytab);
        temp16alen =
            scale_expansion_zeroelimTRI(cytablen, cytab, 2.0 * cdy, temp16a);

        cytaalen   = scale_expansion_zeroelimTRI(4, aa, cdytail, cytaa);
        temp16blen = scale_expansion_zeroelimTRI(cytaalen, cytaa, bdx, temp16b);

        cytbblen   = scale_expansion_zeroelimTRI(4, bb, cdytail, cytbb);
        temp16clen = scale_expansion_zeroelimTRI(cytbblen, cytbb, -adx, temp16c);

        temp32alen = fast_expansion_sum_zeroelimTRI(
            temp16alen, temp16a, temp16blen, temp16b, temp32a);
        temp48len = fast_expansion_sum_zeroelimTRI(
            temp16clen, temp16c, temp32alen, temp32a, temp48);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, temp48len, temp48, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }

    if ((adxtail != 0.0) || (adytail != 0.0))
    {
        if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) ||
            (cdytail != 0.0))
        {
            Two_Product(bdxtail, cdy, ti1, ti0);
            Two_Product(bdx, cdytail, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3]   = u3;
            negate = -bdy;
            Two_Product(cdxtail, negate, ti1, ti0);
            negate = -bdytail;
            Two_Product(cdx, negate, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3]   = v3;
            bctlen = fast_expansion_sum_zeroelimTRI(4, u, 4, v, bct);

            Two_Product(bdxtail, cdytail, ti1, ti0);
            Two_Product(cdxtail, bdytail, tj1, tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, bctt3, bctt[2], bctt[1], bctt[0]);
            bctt[3] = bctt3;
            bcttlen = 4;
        }
        else
        {
            bct[0]  = 0.0;
            bctlen  = 1;
            bctt[0] = 0.0;
            bcttlen = 1;
        }

        if (adxtail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(axtbclen, axtbc, adxtail, temp16a);
            axtbctlen = scale_expansion_zeroelimTRI(bctlen, bct, adxtail, axtbct);
            temp32alen =
                scale_expansion_zeroelimTRI(axtbctlen, axtbct, 2.0 * adx, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (bdytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, cc, adxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, bdytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
            if (cdytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, bb, -adxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, cdytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }

            temp32alen =
                scale_expansion_zeroelimTRI(axtbctlen, axtbct, adxtail, temp32a);
            axtbcttlen =
                scale_expansion_zeroelimTRI(bcttlen, bctt, adxtail, axtbctt);
            temp16alen = scale_expansion_zeroelimTRI(
                axtbcttlen, axtbctt, 2.0 * adx, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(axtbcttlen, axtbctt, adxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
        if (adytail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(aytbclen, aytbc, adytail, temp16a);
            aytbctlen = scale_expansion_zeroelimTRI(bctlen, bct, adytail, aytbct);
            temp32alen =
                scale_expansion_zeroelimTRI(aytbctlen, aytbct, 2.0 * ady, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;

            temp32alen =
                scale_expansion_zeroelimTRI(aytbctlen, aytbct, adytail, temp32a);
            aytbcttlen =
                scale_expansion_zeroelimTRI(bcttlen, bctt, adytail, aytbctt);
            temp16alen = scale_expansion_zeroelimTRI(
                aytbcttlen, aytbctt, 2.0 * ady, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(aytbcttlen, aytbctt, adytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
    }
    if ((bdxtail != 0.0) || (bdytail != 0.0))
    {
        if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) ||
            (adytail != 0.0))
        {
            Two_Product(cdxtail, ady, ti1, ti0);
            Two_Product(cdx, adytail, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3]   = u3;
            negate = -cdy;
            Two_Product(adxtail, negate, ti1, ti0);
            negate = -cdytail;
            Two_Product(adx, negate, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3]   = v3;
            catlen = fast_expansion_sum_zeroelimTRI(4, u, 4, v, cat);

            Two_Product(cdxtail, adytail, ti1, ti0);
            Two_Product(adxtail, cdytail, tj1, tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, catt3, catt[2], catt[1], catt[0]);
            catt[3] = catt3;
            cattlen = 4;
        }
        else
        {
            cat[0]  = 0.0;
            catlen  = 1;
            catt[0] = 0.0;
            cattlen = 1;
        }

        if (bdxtail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(bxtcalen, bxtca, bdxtail, temp16a);
            bxtcatlen = scale_expansion_zeroelimTRI(catlen, cat, bdxtail, bxtcat);
            temp32alen =
                scale_expansion_zeroelimTRI(bxtcatlen, bxtcat, 2.0 * bdx, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (cdytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, aa, bdxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, cdytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
            if (adytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, cc, -bdxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, adytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }

            temp32alen =
                scale_expansion_zeroelimTRI(bxtcatlen, bxtcat, bdxtail, temp32a);
            bxtcattlen =
                scale_expansion_zeroelimTRI(cattlen, catt, bdxtail, bxtcatt);
            temp16alen = scale_expansion_zeroelimTRI(
                bxtcattlen, bxtcatt, 2.0 * bdx, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(bxtcattlen, bxtcatt, bdxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
        if (bdytail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(bytcalen, bytca, bdytail, temp16a);
            bytcatlen = scale_expansion_zeroelimTRI(catlen, cat, bdytail, bytcat);
            temp32alen =
                scale_expansion_zeroelimTRI(bytcatlen, bytcat, 2.0 * bdy, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;

            temp32alen =
                scale_expansion_zeroelimTRI(bytcatlen, bytcat, bdytail, temp32a);
            bytcattlen =
                scale_expansion_zeroelimTRI(cattlen, catt, bdytail, bytcatt);
            temp16alen = scale_expansion_zeroelimTRI(
                bytcattlen, bytcatt, 2.0 * bdy, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(bytcattlen, bytcatt, bdytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
    }
    if ((cdxtail != 0.0) || (cdytail != 0.0))
    {
        if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) ||
            (bdytail != 0.0))
        {
            Two_Product(adxtail, bdy, ti1, ti0);
            Two_Product(adx, bdytail, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u3, u[2], u[1], u[0]);
            u[3]   = u3;
            negate = -ady;
            Two_Product(bdxtail, negate, ti1, ti0);
            negate = -adytail;
            Two_Product(bdx, negate, tj1, tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v3, v[2], v[1], v[0]);
            v[3]   = v3;
            abtlen = fast_expansion_sum_zeroelimTRI(4, u, 4, v, abt);

            Two_Product(adxtail, bdytail, ti1, ti0);
            Two_Product(bdxtail, adytail, tj1, tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, abtt3, abtt[2], abtt[1], abtt[0]);
            abtt[3] = abtt3;
            abttlen = 4;
        }
        else
        {
            abt[0]  = 0.0;
            abtlen  = 1;
            abtt[0] = 0.0;
            abttlen = 1;
        }

        if (cdxtail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(cxtablen, cxtab, cdxtail, temp16a);
            cxtabtlen = scale_expansion_zeroelimTRI(abtlen, abt, cdxtail, cxtabt);
            temp32alen =
                scale_expansion_zeroelimTRI(cxtabtlen, cxtabt, 2.0 * cdx, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (adytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, bb, cdxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, adytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
            if (bdytail != 0.0)
            {
                temp8len = scale_expansion_zeroelimTRI(4, aa, -cdxtail, temp8);
                temp16alen =
                    scale_expansion_zeroelimTRI(temp8len, temp8, bdytail, temp16a);
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, temp16alen, temp16a, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }

            temp32alen =
                scale_expansion_zeroelimTRI(cxtabtlen, cxtabt, cdxtail, temp32a);
            cxtabttlen =
                scale_expansion_zeroelimTRI(abttlen, abtt, cdxtail, cxtabtt);
            temp16alen = scale_expansion_zeroelimTRI(
                cxtabttlen, cxtabtt, 2.0 * cdx, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(cxtabttlen, cxtabtt, cdxtail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
        if (cdytail != 0.0)
        {
            temp16alen =
                scale_expansion_zeroelimTRI(cytablen, cytab, cdytail, temp16a);
            cytabtlen = scale_expansion_zeroelimTRI(abtlen, abt, cdytail, cytabt);
            temp32alen =
                scale_expansion_zeroelimTRI(cytabtlen, cytabt, 2.0 * cdy, temp32a);
            temp48len = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp32alen, temp32a, temp48);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp48len, temp48, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;

            temp32alen =
                scale_expansion_zeroelimTRI(cytabtlen, cytabt, cdytail, temp32a);
            cytabttlen =
                scale_expansion_zeroelimTRI(abttlen, abtt, cdytail, cytabtt);
            temp16alen = scale_expansion_zeroelimTRI(
                cytabttlen, cytabtt, 2.0 * cdy, temp16a);
            temp16blen =
                scale_expansion_zeroelimTRI(cytabttlen, cytabtt, cdytail, temp16b);
            temp32blen = fast_expansion_sum_zeroelimTRI(
                temp16alen, temp16a, temp16blen, temp16b, temp32b);
            temp64len = fast_expansion_sum_zeroelimTRI(
                temp32alen, temp32a, temp32blen, temp32b, temp64);
            finlength = fast_expansion_sum_zeroelimTRI(
                finlength, finnow, temp64len, temp64, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
        }
    }

    return finnow[finlength - 1];
}

double DelaunayTriangle::incircle(struct mesh *m,
                struct behavior *b,
                vertex pa,
                vertex pb,
                vertex pc,
                vertex pd)
{
    double adx, bdx, cdx, ady, bdy, cdy;
    double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
    double alift, blift, clift;
    double det;
    double permanent, errbound;

    m->incirclecount++;

    adx = pa[0] - pd[0];
    bdx = pb[0] - pd[0];
    cdx = pc[0] - pd[0];
    ady = pa[1] - pd[1];
    bdy = pb[1] - pd[1];
    cdy = pc[1] - pd[1];

    bdxcdy = bdx * cdy;
    cdxbdy = cdx * bdy;
    alift  = adx * adx + ady * ady;

    cdxady = cdx * ady;
    adxcdy = adx * cdy;
    blift  = bdx * bdx + bdy * bdy;

    adxbdy = adx * bdy;
    bdxady = bdx * ady;
    clift  = cdx * cdx + cdy * cdy;

    det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) +
          clift * (adxbdy - bdxady);

    permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift +
                (Absolute(cdxady) + Absolute(adxcdy)) * blift +
                (Absolute(adxbdy) + Absolute(bdxady)) * clift;
    errbound = iccerrboundA * permanent;
    if ((det > errbound) || (-det > errbound))
    {
        return det;
    }

    return incircleadaptTRI(pa, pb, pc, pd, permanent);
}

/*****************************************************************************/
/*                                                                           */
/*  orient3d()   Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
/*                                                                           */
/*  Uses exact arithmetic if necessary to ensure a correct answer.  The      */
/*  result returned is the determinant of a matrix.  This determinant is     */
/*  computed adaptively, in the sense that exact arithmetic is used only to  */
/*  the degree it is needed to ensure that the returned value has the        */
/*  correct sign.  Hence, this function is usually quite fast, but will run  */
/*  more slowly when the input points are coplanar or nearly so.             */
/*                                                                           */
/*  See my Robust Predicates paper for details.                              */
/*                                                                           */
/*****************************************************************************/

double DelaunayTriangle::orient3dadapt(vertex pa,
                     vertex pb,
                     vertex pc,
                     vertex pd,
                     double aheight,
                     double bheight,
                     double cheight,
                     double dheight,
                     double permanent)
{
     double adx, bdx, cdx, ady, bdy, cdy, adheight, bdheight, cdheight;
    double det, errbound;

     double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
    double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
    double bc[4], ca[4], ab[4];
     double bc3, ca3, ab3;
    double adet[8], bdet[8], cdet[8];
    int alen, blen, clen;
    double abdet[16];
    int ablen;
    double *finnow, *finother, *finswap;
    double fin1[192], fin2[192];
    int finlength;

    double adxtail, bdxtail, cdxtail;
    double adytail, bdytail, cdytail;
    double adheighttail, bdheighttail, cdheighttail;
     double at_blarge, at_clarge;
     double bt_clarge, bt_alarge;
     double ct_alarge, ct_blarge;
    double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
    int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
     double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
     double adxt_cdy1, adxt_bdy1, bdxt_ady1;
    double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
    double adxt_cdy0, adxt_bdy0, bdxt_ady0;
     double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
     double adyt_cdx1, adyt_bdx1, bdyt_adx1;
    double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
    double adyt_cdx0, adyt_bdx0, bdyt_adx0;
    double bct[8], cat[8], abt[8];
    int bctlen, catlen, abtlen;
     double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
     double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
    double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
    double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
    double u[4], v[12], w[16];
     double u3;
    int vlength, wlength;
    double negate;

     double bvirt;
    double avirt, bround, around;
     double c;
     double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
     double _i, _j, _k;
    double _0;

    adx      = (double)(pa[0] - pd[0]);
    bdx      = (double)(pb[0] - pd[0]);
    cdx      = (double)(pc[0] - pd[0]);
    ady      = (double)(pa[1] - pd[1]);
    bdy      = (double)(pb[1] - pd[1]);
    cdy      = (double)(pc[1] - pd[1]);
    adheight = (double)(aheight - dheight);
    bdheight = (double)(bheight - dheight);
    cdheight = (double)(cheight - dheight);

    Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
    Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
    Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
    bc[3] = bc3;
    alen  = scale_expansion_zeroelimTRI(4, bc, adheight, adet);

    Two_Product(cdx, ady, cdxady1, cdxady0);
    Two_Product(adx, cdy, adxcdy1, adxcdy0);
    Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
    ca[3] = ca3;
    blen  = scale_expansion_zeroelimTRI(4, ca, bdheight, bdet);

    Two_Product(adx, bdy, adxbdy1, adxbdy0);
    Two_Product(bdx, ady, bdxady1, bdxady0);
    Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
    ab[3] = ab3;
    clen  = scale_expansion_zeroelimTRI(4, ab, cdheight, cdet);

    ablen     = fast_expansion_sum_zeroelimTRI(alen, adet, blen, bdet, abdet);
    finlength = fast_expansion_sum_zeroelimTRI(ablen, abdet, clen, cdet, fin1);

    det      = estimateTRI(finlength, fin1);
    errbound = o3derrboundB * permanent;
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
    Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
    Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
    Two_Diff_Tail(pa[1], pd[1], ady, adytail);
    Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
    Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
    Two_Diff_Tail(aheight, dheight, adheight, adheighttail);
    Two_Diff_Tail(bheight, dheight, bdheight, bdheighttail);
    Two_Diff_Tail(cheight, dheight, cdheight, cdheighttail);

    if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) &&
        (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0) &&
        (adheighttail == 0.0) && (bdheighttail == 0.0) && (cdheighttail == 0.0))
    {
        return det;
    }

    errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
    det += (adheight * ((bdx * cdytail + cdy * bdxtail) -
                        (bdy * cdxtail + cdx * bdytail)) +
            adheighttail * (bdx * cdy - bdy * cdx)) +
           (bdheight * ((cdx * adytail + ady * cdxtail) -
                        (cdy * adxtail + adx * cdytail)) +
            bdheighttail * (cdx * ady - cdy * adx)) +
           (cdheight * ((adx * bdytail + bdy * adxtail) -
                        (ady * bdxtail + bdx * adytail)) +
            cdheighttail * (adx * bdy - ady * bdx));
    if ((det >= errbound) || (-det >= errbound))
    {
        return det;
    }

    finnow   = fin1;
    finother = fin2;

    if (adxtail == 0.0)
    {
        if (adytail == 0.0)
        {
            at_b[0] = 0.0;
            at_blen = 1;
            at_c[0] = 0.0;
            at_clen = 1;
        }
        else
        {
            negate = -adytail;
            Two_Product(negate, bdx, at_blarge, at_b[0]);
            at_b[1] = at_blarge;
            at_blen = 2;
            Two_Product(adytail, cdx, at_clarge, at_c[0]);
            at_c[1] = at_clarge;
            at_clen = 2;
        }
    }
    else
    {
        if (adytail == 0.0)
        {
            Two_Product(adxtail, bdy, at_blarge, at_b[0]);
            at_b[1] = at_blarge;
            at_blen = 2;
            negate  = -adxtail;
            Two_Product(negate, cdy, at_clarge, at_c[0]);
            at_c[1] = at_clarge;
            at_clen = 2;
        }
        else
        {
            Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
            Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
            Two_Two_Diff(adxt_bdy1,
                         adxt_bdy0,
                         adyt_bdx1,
                         adyt_bdx0,
                         at_blarge,
                         at_b[2],
                         at_b[1],
                         at_b[0]);
            at_b[3] = at_blarge;
            at_blen = 4;
            Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
            Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
            Two_Two_Diff(adyt_cdx1,
                         adyt_cdx0,
                         adxt_cdy1,
                         adxt_cdy0,
                         at_clarge,
                         at_c[2],
                         at_c[1],
                         at_c[0]);
            at_c[3] = at_clarge;
            at_clen = 4;
        }
    }
    if (bdxtail == 0.0)
    {
        if (bdytail == 0.0)
        {
            bt_c[0] = 0.0;
            bt_clen = 1;
            bt_a[0] = 0.0;
            bt_alen = 1;
        }
        else
        {
            negate = -bdytail;
            Two_Product(negate, cdx, bt_clarge, bt_c[0]);
            bt_c[1] = bt_clarge;
            bt_clen = 2;
            Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
            bt_a[1] = bt_alarge;
            bt_alen = 2;
        }
    }
    else
    {
        if (bdytail == 0.0)
        {
            Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
            bt_c[1] = bt_clarge;
            bt_clen = 2;
            negate  = -bdxtail;
            Two_Product(negate, ady, bt_alarge, bt_a[0]);
            bt_a[1] = bt_alarge;
            bt_alen = 2;
        }
        else
        {
            Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
            Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
            Two_Two_Diff(bdxt_cdy1,
                         bdxt_cdy0,
                         bdyt_cdx1,
                         bdyt_cdx0,
                         bt_clarge,
                         bt_c[2],
                         bt_c[1],
                         bt_c[0]);
            bt_c[3] = bt_clarge;
            bt_clen = 4;
            Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
            Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
            Two_Two_Diff(bdyt_adx1,
                         bdyt_adx0,
                         bdxt_ady1,
                         bdxt_ady0,
                         bt_alarge,
                         bt_a[2],
                         bt_a[1],
                         bt_a[0]);
            bt_a[3] = bt_alarge;
            bt_alen = 4;
        }
    }
    if (cdxtail == 0.0)
    {
        if (cdytail == 0.0)
        {
            ct_a[0] = 0.0;
            ct_alen = 1;
            ct_b[0] = 0.0;
            ct_blen = 1;
        }
        else
        {
            negate = -cdytail;
            Two_Product(negate, adx, ct_alarge, ct_a[0]);
            ct_a[1] = ct_alarge;
            ct_alen = 2;
            Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
            ct_b[1] = ct_blarge;
            ct_blen = 2;
        }
    }
    else
    {
        if (cdytail == 0.0)
        {
            Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
            ct_a[1] = ct_alarge;
            ct_alen = 2;
            negate  = -cdxtail;
            Two_Product(negate, bdy, ct_blarge, ct_b[0]);
            ct_b[1] = ct_blarge;
            ct_blen = 2;
        }
        else
        {
            Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
            Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
            Two_Two_Diff(cdxt_ady1,
                         cdxt_ady0,
                         cdyt_adx1,
                         cdyt_adx0,
                         ct_alarge,
                         ct_a[2],
                         ct_a[1],
                         ct_a[0]);
            ct_a[3] = ct_alarge;
            ct_alen = 4;
            Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
            Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
            Two_Two_Diff(cdyt_bdx1,
                         cdyt_bdx0,
                         cdxt_bdy1,
                         cdxt_bdy0,
                         ct_blarge,
                         ct_b[2],
                         ct_b[1],
                         ct_b[0]);
            ct_b[3] = ct_blarge;
            ct_blen = 4;
        }
    }

    bctlen  = fast_expansion_sum_zeroelimTRI(bt_clen, bt_c, ct_blen, ct_b, bct);
    wlength = scale_expansion_zeroelimTRI(bctlen, bct, adheight, w);
    finlength =
        fast_expansion_sum_zeroelimTRI(finlength, finnow, wlength, w, finother);
    finswap  = finnow;
    finnow   = finother;
    finother = finswap;

    catlen  = fast_expansion_sum_zeroelimTRI(ct_alen, ct_a, at_clen, at_c, cat);
    wlength = scale_expansion_zeroelimTRI(catlen, cat, bdheight, w);
    finlength =
        fast_expansion_sum_zeroelimTRI(finlength, finnow, wlength, w, finother);
    finswap  = finnow;
    finnow   = finother;
    finother = finswap;

    abtlen  = fast_expansion_sum_zeroelimTRI(at_blen, at_b, bt_alen, bt_a, abt);
    wlength = scale_expansion_zeroelimTRI(abtlen, abt, cdheight, w);
    finlength =
        fast_expansion_sum_zeroelimTRI(finlength, finnow, wlength, w, finother);
    finswap  = finnow;
    finnow   = finother;
    finother = finswap;

    if (adheighttail != 0.0)
    {
        vlength   = scale_expansion_zeroelimTRI(4, bc, adheighttail, v);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, vlength, v, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (bdheighttail != 0.0)
    {
        vlength   = scale_expansion_zeroelimTRI(4, ca, bdheighttail, v);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, vlength, v, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (cdheighttail != 0.0)
    {
        vlength   = scale_expansion_zeroelimTRI(4, ab, cdheighttail, v);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, vlength, v, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }

    if (adxtail != 0.0)
    {
        if (bdytail != 0.0)
        {
            Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
            Two_One_Product(
                adxt_bdyt1, adxt_bdyt0, cdheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (cdheighttail != 0.0)
            {
                Two_One_Product(
                    adxt_bdyt1, adxt_bdyt0, cdheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
        if (cdytail != 0.0)
        {
            negate = -adxtail;
            Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
            Two_One_Product(
                adxt_cdyt1, adxt_cdyt0, bdheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (bdheighttail != 0.0)
            {
                Two_One_Product(
                    adxt_cdyt1, adxt_cdyt0, bdheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
    }
    if (bdxtail != 0.0)
    {
        if (cdytail != 0.0)
        {
            Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
            Two_One_Product(
                bdxt_cdyt1, bdxt_cdyt0, adheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (adheighttail != 0.0)
            {
                Two_One_Product(
                    bdxt_cdyt1, bdxt_cdyt0, adheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
        if (adytail != 0.0)
        {
            negate = -bdxtail;
            Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
            Two_One_Product(
                bdxt_adyt1, bdxt_adyt0, cdheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (cdheighttail != 0.0)
            {
                Two_One_Product(
                    bdxt_adyt1, bdxt_adyt0, cdheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
    }
    if (cdxtail != 0.0)
    {
        if (adytail != 0.0)
        {
            Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
            Two_One_Product(
                cdxt_adyt1, cdxt_adyt0, bdheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (bdheighttail != 0.0)
            {
                Two_One_Product(
                    cdxt_adyt1, cdxt_adyt0, bdheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
        if (bdytail != 0.0)
        {
            negate = -cdxtail;
            Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
            Two_One_Product(
                cdxt_bdyt1, cdxt_bdyt0, adheight, u3, u[2], u[1], u[0]);
            u[3] = u3;
            finlength =
                fast_expansion_sum_zeroelimTRI(finlength, finnow, 4, u, finother);
            finswap  = finnow;
            finnow   = finother;
            finother = finswap;
            if (adheighttail != 0.0)
            {
                Two_One_Product(
                    cdxt_bdyt1, cdxt_bdyt0, adheighttail, u3, u[2], u[1], u[0]);
                u[3]      = u3;
                finlength = fast_expansion_sum_zeroelimTRI(
                    finlength, finnow, 4, u, finother);
                finswap  = finnow;
                finnow   = finother;
                finother = finswap;
            }
        }
    }

    if (adheighttail != 0.0)
    {
        wlength   = scale_expansion_zeroelimTRI(bctlen, bct, adheighttail, w);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, wlength, w, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (bdheighttail != 0.0)
    {
        wlength   = scale_expansion_zeroelimTRI(catlen, cat, bdheighttail, w);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, wlength, w, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }
    if (cdheighttail != 0.0)
    {
        wlength   = scale_expansion_zeroelimTRI(abtlen, abt, cdheighttail, w);
        finlength = fast_expansion_sum_zeroelimTRI(
            finlength, finnow, wlength, w, finother);
        finswap  = finnow;
        finnow   = finother;
        finother = finswap;
    }

    return finnow[finlength - 1];
}

double DelaunayTriangle::orient3d(struct mesh *m,
                struct behavior *b,
                vertex pa,
                vertex pb,
                vertex pc,
                vertex pd,
                double aheight,
                double bheight,
                double cheight,
                double dheight)
{
    double adx, bdx, cdx, ady, bdy, cdy, adheight, bdheight, cdheight;
    double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
    double det;
    double permanent, errbound;

    m->orient3dcount++;

    adx      = pa[0] - pd[0];
    bdx      = pb[0] - pd[0];
    cdx      = pc[0] - pd[0];
    ady      = pa[1] - pd[1];
    bdy      = pb[1] - pd[1];
    cdy      = pc[1] - pd[1];
    adheight = aheight - dheight;
    bdheight = bheight - dheight;
    cdheight = cheight - dheight;

    bdxcdy = bdx * cdy;
    cdxbdy = cdx * bdy;

    cdxady = cdx * ady;
    adxcdy = adx * cdy;

    adxbdy = adx * bdy;
    bdxady = bdx * ady;

    det = adheight * (bdxcdy - cdxbdy) + bdheight * (cdxady - adxcdy) +
          cdheight * (adxbdy - bdxady);

    permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adheight) +
                (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdheight) +
                (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdheight);
    errbound = o3derrboundA * permanent;
    if ((det > errbound) || (-det > errbound))
    {
        return det;
    }

    return orient3dadapt(
        pa, pb, pc, pd, aheight, bheight, cheight, dheight, permanent);
}

/*****************************************************************************/
/*                                                                           */
/*  nonregular()   Return a positive value if the point pd is incompatible   */
/*                 with the circle or plane passing through pa, pb, and pc   */
/*                 (meaning that pd is inside the circle or below the        */
/*                 plane); a negative value if it is compatible; and zero if */
/*                 the four points are cocircular/coplanar.  The points pa,  */
/*                 pb, and pc must be in counterclockwise order, or the sign */
/*                 of the result will be reversed.                           */
/*                                                                           */
/*  If the -w switch is used, the points are lifted onto the parabolic       */
/*  lifting map, then they are dropped according to their weights, then the  */
/*  3D orientation test is applied.  If the -W switch is used, the points'   */
/*  heights are already provided, so the 3D orientation test is applied      */
/*  directly.  If neither switch is used, the incircle test is applied.      */
/*                                                                           */
/*****************************************************************************/

double DelaunayTriangle::nonregular(struct mesh *m,
                  struct behavior *b,
                  vertex pa,
                  vertex pb,
                  vertex pc,
                  vertex pd)
{
    if (b->weighted == 0)
    {
        return incircle(m, b, pa, pb, pc, pd);
    }
    else if (b->weighted == 1)
    {
        return orient3d(m,
                        b,
                        pa,
                        pb,
                        pc,
                        pd,
                        pa[0] * pa[0] + pa[1] * pa[1] - pa[2],
                        pb[0] * pb[0] + pb[1] * pb[1] - pb[2],
                        pc[0] * pc[0] + pc[1] * pc[1] - pc[2],
                        pd[0] * pd[0] + pd[1] * pd[1] - pd[2]);
    }
    else
    {
        return orient3d(m, b, pa, pb, pc, pd, pa[2], pb[2], pc[2], pd[2]);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  findcircumcenter()   Find the circumcenter of a triangle.                */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  (barycentric) coordinates.  The xi-eta coordinate system is defined in   */
/*  terms of the triangle:  the origin of the triangle is the origin of the  */
/*  coordinate system; the destination of the triangle is one unit along the */
/*  xi axis; and the apex of the triangle is one unit along the eta axis.    */
/*  This procedure also returns the square of the length of the triangle's   */
/*  shortest edge.                                                           */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::findcircumcenter(struct mesh *m,
                      struct behavior *b,
                      vertex torg,
                      vertex tdest,
                      vertex tapex,
                      vertex circumcenter,
                      double *xi,
                      double *eta,
                      int offcenter)
{
    double xdo, ydo, xao, yao;
    double dodist, aodist, dadist;
    double denominator;
    double dx, dy, dxoff, dyoff;

    m->circumcentercount++;

    /* Compute the circumcenter of the triangle. */
    xdo    = tdest[0] - torg[0];
    ydo    = tdest[1] - torg[1];
    xao    = tapex[0] - torg[0];
    yao    = tapex[1] - torg[1];
    dodist = xdo * xdo + ydo * ydo;
    aodist = xao * xao + yao * yao;
    dadist = (tdest[0] - tapex[0]) * (tdest[0] - tapex[0]) +
             (tdest[1] - tapex[1]) * (tdest[1] - tapex[1]);

    /* Use the counterclockwise() routine to ensure a positive (and */
    /*   reasonably accurate) result, avoiding any possibility of   */
    /*   division by zero.                                          */
    denominator = 0.5 / counterclockwise(m, b, tdest, tapex, torg);
    /* Don't count the above as an orientation test. */
    m->counterclockcount--;

    dx = (yao * dodist - ydo * aodist) * denominator;
    dy = (xdo * aodist - xao * dodist) * denominator;

    /* Find the (squared) length of the triangle's shortest edge.  This   */
    /*   serves as a conservative estimateTRI of the insertion radius of the */
    /*   circumcenter's parent.  The estimateTRI is used to ensure that      */
    /*   the algorithm terminates even if very small angles appear in     */
    /*   the input PSLG.                                                  */
    if ((dodist < aodist) && (dodist < dadist))
    {
        if (offcenter && (b->offconstant > 0.0))
        {
            /* Find the position of the off-center, as described by Alper Ungor.
             */
            dxoff = 0.5 * xdo - b->offconstant * ydo;
            dyoff = 0.5 * ydo + b->offconstant * xdo;
            /* If the off-center is closer to the origin than the */
            /*   circumcenter, use the off-center instead.        */
            if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
            {
                dx = dxoff;
                dy = dyoff;
            }
        }
    }
    else if (aodist < dadist)
    {
        if (offcenter && (b->offconstant > 0.0))
        {
            dxoff = 0.5 * xao + b->offconstant * yao;
            dyoff = 0.5 * yao - b->offconstant * xao;
            /* If the off-center is closer to the origin than the */
            /*   circumcenter, use the off-center instead.        */
            if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
            {
                dx = dxoff;
                dy = dyoff;
            }
        }
    }
    else
    {
        if (offcenter && (b->offconstant > 0.0))
        {
            dxoff = 0.5 * (tapex[0] - tdest[0]) -
                    b->offconstant * (tapex[1] - tdest[1]);
            dyoff = 0.5 * (tapex[1] - tdest[1]) +
                    b->offconstant * (tapex[0] - tdest[0]);
            /* If the off-center is closer to the destination than the */
            /*   circumcenter, use the off-center instead.             */
            if (dxoff * dxoff + dyoff * dyoff <
                (dx - xdo) * (dx - xdo) + (dy - ydo) * (dy - ydo))
            {
                dx = xdo + dxoff;
                dy = ydo + dyoff;
            }
        }
    }

    circumcenter[0] = torg[0] + dx;
    circumcenter[1] = torg[1] + dy;

    /* To interpolate vertex attributes for the new vertex inserted at */
    /*   the circumcenter, define a coordinate system with a xi-axis,  */
    /*   directed from the triangle's origin to its destination, and   */
    /*   an eta-axis, directed from its origin to its apex.            */
    /*   Calculate the xi and eta coordinates of the circumcenter.     */
    *xi  = (yao * dx - xao * dy) * (2.0 * denominator);
    *eta = (xdo * dy - ydo * dx) * (2.0 * denominator);
}

/**                                                                         **/
/**                                                                         **/
/********* Geometric primitives end here                             *********/

/*****************************************************************************/
/*                                                                           */
/*  triangleinit()   Initialize some variables.                              */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::triangleinit(struct mesh *m)
{
    poolzero(&m->vertices);
    poolzero(&m->triangles);
    poolzero(&m->subsegs);
    poolzero(&m->viri);
    poolzero(&m->badsubsegs);
    poolzero(&m->badtriangles);
    poolzero(&m->flipstackers);
    poolzero(&m->splaynodes);

    m->recenttri.tri = (triangle *)NULL; /* No triangle has been visited yet. */
    m->undeads       = 0;                /* No eliminated input vertices yet. */
    m->samples       = 1; /* Point location should take at least one sample. */
    m->checksegments = 0; /* There are no segments in the triangulation yet. */
    m->checkquality  = 0; /* The quality triangulation stage has not begun. */
    m->incirclecount = m->counterclockcount = m->orient3dcount = 0;
    m->hyperbolacount = m->circletopcount = m->circumcentercount = 0;
    randomseed                                                   = 1;

    exactinit(); /* Initialize exact arithmetic constants. */
}

/*****************************************************************************/
/*                                                                           */
/*  randomnation()   Generate a random number between 0 and `choices' - 1.   */
/*                                                                           */
/*  This is a simple linear congruential random number generator.  Hence, it */
/*  is a bad random number generator, but good enough for most randomized    */
/*  geometric algorithms.                                                    */
/*                                                                           */
/*****************************************************************************/

unsigned long DelaunayTriangle::randomnation(unsigned int choices)
{
    randomseed = (randomseed * 1366l + 150889l) % 714025l;
    return randomseed / (714025l / choices + 1);
}

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtriang()   Add a bad triangle data structure to the end of a   */
/*                       queue.                                              */
/*                                                                           */
/*  The queue is actually a set of 4096 queues.  I use multiple queues to    */
/*  give priority to smaller angles.  I originally implemented a heap, but   */
/*  the queues are faster by a larger margin than I'd suspected.             */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::enqueuebadtriang(struct mesh *m,
                      struct behavior *b,
                      struct badtriang *badtri)
{
    double length, multiplier;
    int exponent, expincrement;
    int queuenumber;
    int posexponent;
    int i;

    /* Determine the appropriate queue to put the bad triangle into.    */
    /*   Recall that the key is the square of its shortest edge length. */
    if (badtri->key >= 1.0)
    {
        length      = badtri->key;
        posexponent = 1;
    }
    else
    {
        /* `badtri->key' is 2.0 to a negative exponent, so we'll record that */
        /*   fact and use the reciprocal of `badtri->key', which is > 1.0.   */
        length      = 1.0 / badtri->key;
        posexponent = 0;
    }
    /* `length' is approximately 2.0 to what exponent?  The following code */
    /*   determines the answer in time logarithmic in the exponent.        */
    exponent = 0;
    while (length > 2.0)
    {
        /* Find an approximation by repeated squaring of two. */
        expincrement = 1;
        multiplier   = 0.5;
        while (length * multiplier * multiplier > 1.0)
        {
            expincrement *= 2;
            multiplier *= multiplier;
        }
        /* Reduce the value of `length', then iterate if necessary. */
        exponent += expincrement;
        length *= multiplier;
    }
    /* `length' is approximately squareroot(2.0) to what exponent? */
    exponent = 2.0 * exponent + (length > SQUAREROOTTWO);
    /* `exponent' is now in the range 0...2047 for IEEE double precision.   */
    /*   Choose a queue in the range 0...4095.  The shortest edges have the */
    /*   highest priority (queue 4095).                                     */
    if (posexponent)
    {
        queuenumber = 2047 - exponent;
    }
    else
    {
        queuenumber = 2048 + exponent;
    }

    /* Are we inserting into an empty queue? */
    if (m->queuefront[queuenumber] == (struct badtriang *)NULL)
    {
        /* Yes, we are inserting into an empty queue.     */
        /*   Will this become the highest-priority queue? */
        if (queuenumber > m->firstnonemptyq)
        {
            /* Yes, this is the highest-priority queue. */
            m->nextnonemptyq[queuenumber] = m->firstnonemptyq;
            m->firstnonemptyq             = queuenumber;
        }
        else
        {
            /* No, this is not the highest-priority queue. */
            /*   Find the queue with next higher priority. */
            i = queuenumber + 1;
            while (m->queuefront[i] == (struct badtriang *)NULL)
            {
                i++;
            }
            /* Mark the newly nonempty queue as following a higher-priority
             * queue. */
            m->nextnonemptyq[queuenumber] = m->nextnonemptyq[i];
            m->nextnonemptyq[i]           = queuenumber;
        }
        /* Put the bad triangle at the beginning of the (empty) queue. */
        m->queuefront[queuenumber] = badtri;
    }
    else
    {
        /* Add the bad triangle to the end of an already nonempty queue. */
        m->queuetail[queuenumber]->nexttriang = badtri;
    }
    /* Maintain a pointer to the last triangle of the queue. */
    m->queuetail[queuenumber] = badtri;
    /* Newly enqueued bad triangle has no successor in the queue. */
    badtri->nexttriang = (struct badtriang *)NULL;
}

/*****************************************************************************/
/*                                                                           */
/*  enqueuebadtri()   Add a bad triangle to the end of a queue.              */
/*                                                                           */
/*  Allocates a badtriang data structure for the triangle, then passes it to */
/*  enqueuebadtriang().                                                      */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::enqueuebadtri(struct mesh *m,
                   struct behavior *b,
                   struct otri *enqtri,
                   double minedge,
                   vertex enqapex,
                   vertex enqorg,
                   vertex enqdest)
{
    struct badtriang *newbad;

    /* Allocate space for the bad triangle. */
    newbad             = (struct badtriang *)poolalloc(&m->badtriangles);
    newbad->poortri    = encode(*enqtri);
    newbad->key        = minedge;
    newbad->triangapex = enqapex;
    newbad->triangorg  = enqorg;
    newbad->triangdest = enqdest;
    enqueuebadtriang(m, b, newbad);
}

/*****************************************************************************/
/*                                                                           */
/*  dequeuebadtriang()   Remove a triangle from the front of the queue.      */
/*                                                                           */
/*****************************************************************************/

struct badtriang *DelaunayTriangle::dequeuebadtriang(struct mesh *m)
{
    struct badtriang *result;

    /* If no queues are nonempty, return NULL. */
    if (m->firstnonemptyq < 0)
    {
        return (struct badtriang *)NULL;
    }
    /* Find the first triangle of the highest-priority queue. */
    result = m->queuefront[m->firstnonemptyq];
    /* Remove the triangle from the queue. */
    m->queuefront[m->firstnonemptyq] = result->nexttriang;
    /* If this queue is now empty, note the new highest-priority */
    /*   nonempty queue.                                         */
    if (result == m->queuetail[m->firstnonemptyq])
    {
        m->firstnonemptyq = m->nextnonemptyq[m->firstnonemptyq];
    }
    return result;
}

/*****************************************************************************/
/*                                                                           */
/*  checkseg4encroach()   Check a subsegment to see if it is encroached; add */
/*                        it to the list if it is.                           */
/*                                                                           */
/*  A subsegment is encroached if there is a vertex in its diametral lens.   */
/*  For Ruppert's algorithm (-D switch), the "diametral lens" is the         */
/*  diametral circle.  For Chew's algorithm (default), the diametral lens is */
/*  just big enough to enclose two isosceles triangles whose bases are the   */
/*  subsegment.  Each of the two isosceles triangles has two angles equal    */
/*  to `b->minangle'.                                                        */
/*                                                                           */
/*  Chew's algorithm does not require diametral lenses at all--but they save */
/*  time.  Any vertex inside a subsegment's diametral lens implies that the  */
/*  triangle adjoining the subsegment will be too skinny, so it's only a     */
/*  matter of time before the encroaching vertex is deleted by Chew's        */
/*  algorithm.  It's faster to simply not insert the doomed vertex in the    */
/*  first place, which is why I use diametral lenses with Chew's algorithm.  */
/*                                                                           */
/*  Returns a nonzero value if the subsegment is encroached.                 */
/*                                                                           */
/*****************************************************************************/

int DelaunayTriangle::checkseg4encroach(struct mesh *m,
                      struct behavior *b,
                      struct osub *testsubseg)
{
    struct otri neighbortri;
    struct osub testsym;
    struct badsubseg *encroachedseg;
    double dotproduct;
    int encroached;
    int sides;
    vertex eorg, edest, eapex;
    triangle ptr; /* Temporary variable used by stpivot(). */

    encroached = 0;
    sides      = 0;

    sorg(*testsubseg, eorg);
    sdest(*testsubseg, edest);
    /* Check one neighbor of the subsegment. */
    stpivot(*testsubseg, neighbortri);
    /* Does the neighbor exist, or is this a boundary edge? */
    if (neighbortri.tri != m->dummytri)
    {
        sides++;
        /* Find a vertex opposite this subsegment. */
        apex(neighbortri, eapex);
        /* Check whether the apex is in the diametral lens of the subsegment */
        /*   (the diametral circle if `conformdel' is set).  A dot product   */
        /*   of two sides of the triangle is used to check whether the angle */
        /*   at the apex is greater than (180 - 2 `minangle') degrees (for   */
        /*   lenses; 90 degrees for diametral circles).                      */
        dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                     (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
        if (dotproduct < 0.0)
        {
            if (dotproduct * dotproduct >=
                 (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
                     ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
                      (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
                     ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
                      (edest[1] - eapex[1]) * (edest[1] - eapex[1])))
            {
                encroached = 1;
            }
        }
    }
    /* Check the other neighbor of the subsegment. */
    ssym(*testsubseg, testsym);
    stpivot(testsym, neighbortri);
    /* Does the neighbor exist, or is this a boundary edge? */
    if (neighbortri.tri != m->dummytri)
    {
        sides++;
        /* Find the other vertex opposite this subsegment. */
        apex(neighbortri, eapex);
        /* Check whether the apex is in the diametral lens of the subsegment */
        /*   (or the diametral circle, if `conformdel' is set).              */
        dotproduct = (eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                     (eorg[1] - eapex[1]) * (edest[1] - eapex[1]);
        if (dotproduct < 0.0)
        {
            if (dotproduct * dotproduct >=
                 (2.0 * b->goodangle - 1.0) * (2.0 * b->goodangle - 1.0) *
                     ((eorg[0] - eapex[0]) * (eorg[0] - eapex[0]) +
                      (eorg[1] - eapex[1]) * (eorg[1] - eapex[1])) *
                     ((edest[0] - eapex[0]) * (edest[0] - eapex[0]) +
                      (edest[1] - eapex[1]) * (edest[1] - eapex[1])))
            {
                encroached += 2;
            }
        }
    }

    if (encroached && (!b->nobisect || ((b->nobisect == 1) && (sides == 2))))
    {
        /* Add the subsegment to the list of encroached subsegments. */
        /*   Be sure to get the orientation right.                   */
        encroachedseg = (struct badsubseg *)poolalloc(&m->badsubsegs);
        if (encroached == 1)
        {
            encroachedseg->encsubseg  = sencode(*testsubseg);
            encroachedseg->subsegorg  = eorg;
            encroachedseg->subsegdest = edest;
        }
        else
        {
            encroachedseg->encsubseg  = sencode(testsym);
            encroachedseg->subsegorg  = edest;
            encroachedseg->subsegdest = eorg;
        }
    }

    return encroached;
}

/*****************************************************************************/
/*                                                                           */
/*  testtriangle()   Test a triangle for quality and size.                   */
/*                                                                           */
/*  Tests a triangle to see if it satisfies the minimum angle condition and  */
/*  the maximum area condition.  Triangles that aren't up to spec are added  */
/*  to the bad triangle queue.                                               */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::testtriangle(struct mesh *m, struct behavior *b, struct otri *testtri)
{
    struct otri tri1, tri2;
    struct osub testsub;
    vertex torg, tdest, tapex;
    vertex base1, base2;
    vertex org1, dest1, org2, dest2;
    vertex joinvertex;
    double dxod, dyod, dxda, dyda, dxao, dyao;
    double dxod2, dyod2, dxda2, dyda2, dxao2, dyao2;
    double apexlen, orglen, destlen, minedge;
    double angle;
    double area=0.0;
    double dist1, dist2;
    subseg sptr;  /* Temporary variable used by tspivot(). */
    triangle ptr; /* Temporary variable used by oprev() and dnext(). */

    org(*testtri, torg);
    dest(*testtri, tdest);
    apex(*testtri, tapex);
    dxod  = torg[0] - tdest[0];
    dyod  = torg[1] - tdest[1];
    dxda  = tdest[0] - tapex[0];
    dyda  = tdest[1] - tapex[1];
    dxao  = tapex[0] - torg[0];
    dyao  = tapex[1] - torg[1];
    dxod2 = dxod * dxod;
    dyod2 = dyod * dyod;
    dxda2 = dxda * dxda;
    dyda2 = dyda * dyda;
    dxao2 = dxao * dxao;
    dyao2 = dyao * dyao;
    /* Find the lengths of the triangle's three edges. */
    apexlen = dxod2 + dyod2;
    orglen  = dxda2 + dyda2;
    destlen = dxao2 + dyao2;

    if ((apexlen < orglen) && (apexlen < destlen))
    {
        /* The edge opposite the apex is shortest. */
        minedge = apexlen;
        /* Find the square of the cosine of the angle at the apex. */
        angle = dxda * dxao + dyda * dyao;
        angle = angle * angle / (orglen * destlen);
        base1 = torg;
        base2 = tdest;
        otricopy(*testtri, tri1);
    }
    else if (orglen < destlen)
    {
        /* The edge opposite the origin is shortest. */
        minedge = orglen;
        /* Find the square of the cosine of the angle at the origin. */
        angle = dxod * dxao + dyod * dyao;
        angle = angle * angle / (apexlen * destlen);
        base1 = tdest;
        base2 = tapex;
        lnext(*testtri, tri1);
    }
    else
    {
        /* The edge opposite the destination is shortest. */
        minedge = destlen;
        /* Find the square of the cosine of the angle at the destination. */
        angle = dxod * dxda + dyod * dyda;
        angle = angle * angle / (apexlen * orglen);
        base1 = tapex;
        base2 = torg;
        lprev(*testtri, tri1);
    }

    if (b->usertest)
    {
        if (b->usertest)
        {
            /* Check whether the user thinks this triangle is too large. */
            if (triunsuitable(torg, tdest, tapex, area))
            {
                enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
                return;
            }
        }
    }

    /* Check whether the angle is smaller than permitted. */
    if (angle > b->goodangle)
    {
        /* Use the rules of Miller, Pav, and Walkington to decide that certain
         */
        /*   triangles should not be split, even if they have bad angles. */
        /*   A skinny triangle is not split if its shortest edge subtends a */
        /*   small input angle, and both endpoints of the edge lie on a */
        /*   concentric circular shell.  For convenience, I make a small */
        /*   adjustment to that rule:  I check if the endpoints of the edge */
        /*   both lie in segment interiors, equidistant from the apex where */
        /*   the two segments meet. */
        /* First, check if both points lie in segment interiors. */
        if ((vertextype(base1) == SEGMENTVERTEX) &&
            (vertextype(base2) == SEGMENTVERTEX))
        {
            /* Check if both points lie in a common segment.  If they do, the */
            /*   skinny triangle is enqueued to be split as usual.            */
            tspivot(tri1, testsub);
            if (testsub.ss == m->dummysub)
            {
                /* No common segment.  Find a subsegment that contains `torg'.
                 */
                otricopy(tri1, tri2);
                do
                {
                    oprevself(tri1);
                    tspivot(tri1, testsub);
                } while (testsub.ss == m->dummysub);
                /* Find the endpoints of the containing segment. */
                segorg(testsub, org1);
                segdest(testsub, dest1);
                /* Find a subsegment that contains `tdest'. */
                do
                {
                    dnextself(tri2);
                    tspivot(tri2, testsub);
                } while (testsub.ss == m->dummysub);
                /* Find the endpoints of the containing segment. */
                segorg(testsub, org2);
                segdest(testsub, dest2);
                /* Check if the two containing segments have an endpoint in
                 * common. */
                joinvertex = (vertex)NULL;
                if ((dest1[0] == org2[0]) && (dest1[1] == org2[1]))
                {
                    joinvertex = dest1;
                }
                else if ((org1[0] == dest2[0]) && (org1[1] == dest2[1]))
                {
                    joinvertex = org1;
                }
                if (joinvertex != (vertex)NULL)
                {
                    /* Compute the distance from the common endpoint (of the two
                     */
                    /*   segments) to each of the endpoints of the shortest
                     * edge. */
                    dist1 = ((base1[0] - joinvertex[0]) *
                                 (base1[0] - joinvertex[0]) +
                             (base1[1] - joinvertex[1]) *
                                 (base1[1] - joinvertex[1]));
                    dist2 = ((base2[0] - joinvertex[0]) *
                                 (base2[0] - joinvertex[0]) +
                             (base2[1] - joinvertex[1]) *
                                 (base2[1] - joinvertex[1]));
                    /* If the two distances are equal, don't split the triangle.
                     */
                    if ((dist1 < 1.001 * dist2) && (dist1 > 0.999 * dist2))
                    {
                        /* Return now to avoid enqueueing the bad triangle. */
                        return;
                    }
                }
            }
        }

        /* Add this triangle to the list of bad triangles. */
        enqueuebadtri(m, b, testtri, minedge, tapex, torg, tdest);
    }
}

/**                                                                         **/
/**                                                                         **/
/********* Mesh quality testing routines end here                    *********/

/********* Point location routines begin here                        *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  makevertexmap()   Construct a mapping from vertices to triangles to      */
/*                    improve the speed of point location for segment        */
/*                    insertion.                                             */
/*                                                                           */
/*  Traverses all the triangles, and provides each corner of each triangle   */
/*  with a pointer to that triangle.  Of course, pointers will be            */
/*  overwritten by other pointers because (almost) each vertex is a corner   */
/*  of several triangles, but in the end every vertex will point to some     */
/*  triangle that contains it.                                               */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::makevertexmap(struct mesh *m, struct behavior *b)
{
    struct otri triangleloop;
    vertex triorg;

    traversalinit(&m->triangles);
    triangleloop.tri = triangletraverse(m);
    while (triangleloop.tri != (triangle *)NULL)
    {
        /* Check all three vertices of the triangle. */
        for (triangleloop.orient = 0; triangleloop.orient < 3;
             triangleloop.orient++)
        {
            org(triangleloop, triorg);
            setvertex2tri(triorg, encode(triangleloop));
        }
        triangleloop.tri = triangletraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  preciselocate()   Find a triangle or edge containing a given point.      */
/*                                                                           */
/*  Begins its search from `searchtri'.  It is important that `searchtri'    */
/*  be a handle with the property that `searchpoint' is strictly to the left */
/*  of the edge denoted by `searchtri', or is collinear with that edge and   */
/*  does not intersect that edge.  (In particular, `searchpoint' should not  */
/*  be the origin or destination of that edge.)                              */
/*                                                                           */
/*  These conditions are imposed because preciselocate() is normally used in */
/*  one of two situations:                                                   */
/*                                                                           */
/*  (1)  To try to find the location to insert a new point.  Normally, we    */
/*       know an edge that the point is strictly to the left of.  In the     */
/*       incremental Delaunay algorithm, that edge is a bounding box edge.   */
/*       In Ruppert's Delaunay refinement algorithm for quality meshing,     */
/*       that edge is the shortest edge of the triangle whose circumcenter   */
/*       is being inserted.                                                  */
/*                                                                           */
/*  (2)  To try to find an existing point.  In this case, any edge on the    */
/*       convex hull is a good starting edge.  You must screen out the       */
/*       possibility that the vertex sought is an endpoint of the starting   */
/*       edge before you call preciselocate().                               */
/*                                                                           */
/*  On completion, `searchtri' is a triangle that contains `searchpoint'.    */
/*                                                                           */
/*  This implementation differs from that given by Guibas and Stolfi.  It    */
/*  walks from triangle to triangle, crossing an edge only if `searchpoint'  */
/*  is on the other side of the line containing that edge.  After entering   */
/*  a triangle, there are two edges by which one can leave that triangle.    */
/*  If both edges are valid (`searchpoint' is on the other side of both      */
/*  edges), one of the two is chosen by drawing a line perpendicular to      */
/*  the entry edge (whose endpoints are `forg' and `fdest') passing through  */
/*  `fapex'.  Depending on which side of this perpendicular `searchpoint'    */
/*  falls on, an exit edge is chosen.                                        */
/*                                                                           */
/*  This implementation is empirically faster than the Guibas and Stolfi     */
/*  point location routine (which I originally used), which tends to spiral  */
/*  in toward its target.                                                    */
/*                                                                           */
/*  Returns ONVERTEX if the point lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the point lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the point lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the point lies strictly within a triangle.         */
/*  `searchtri' is a handle on the triangle that contains the point.         */
/*                                                                           */
/*  Returns OUTSIDE if the point lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the point is to the right of.  This might      */
/*  occur when the circumcenter of a triangle falls just slightly outside    */
/*  the mesh due to floating-point roundoff error.  It also occurs when      */
/*  seeking a hole or region point that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  If `stopatsubsegment' is nonzero, the search will stop if it tries to    */
/*  walk through a subsegment, and will return OUTSIDE.                      */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*  However, it can still be used to find the circumcenter of a triangle, as */
/*  long as the search is begun from the triangle in question.               */
/*                                                                           */
/*****************************************************************************/

enum locateresult DelaunayTriangle::preciselocate(struct mesh *m,
                                struct behavior *b,
                                vertex searchpoint,
                                struct otri *searchtri,
                                int stopatsubsegment)
{
    struct otri backtracktri;
    struct osub checkedge;
    vertex forg, fdest, fapex;
    double orgorient, destorient;
    int moveleft;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Where are we? */
    org(*searchtri, forg);
    dest(*searchtri, fdest);
    apex(*searchtri, fapex);
    while (1)
    {
        /* Check whether the apex is the point we seek. */
        if ((fapex[0] == searchpoint[0]) && (fapex[1] == searchpoint[1]))
        {
            lprevself(*searchtri);
            return ONVERTEX;
        }
        /* Does the point lie on the other side of the line defined by the */
        /*   triangle edge opposite the triangle's destination?            */
        destorient = counterclockwise(m, b, forg, fapex, searchpoint);
        /* Does the point lie on the other side of the line defined by the */
        /*   triangle edge opposite the triangle's origin?                 */
        orgorient = counterclockwise(m, b, fapex, fdest, searchpoint);
        if (destorient > 0.0)
        {
            if (orgorient > 0.0)
            {
                /* Move left if the inner product of (fapex - searchpoint) and
                 */
                /*   (fdest - forg) is positive.  This is equivalent to drawing
                 */
                /*   a line perpendicular to the line (forg, fdest) and passing
                 */
                /*   through `fapex', and determining which side of this line */
                /*   `searchpoint' falls on. */
                moveleft =
                    (fapex[0] - searchpoint[0]) * (fdest[0] - forg[0]) +
                        (fapex[1] - searchpoint[1]) * (fdest[1] - forg[1]) >
                    0.0;
            }
            else
            {
                moveleft = 1;
            }
        }
        else
        {
            if (orgorient > 0.0)
            {
                moveleft = 0;
            }
            else
            {
                /* The point we seek must be on the boundary of or inside this
                 */
                /*   triangle. */
                if (destorient == 0.0)
                {
                    lprevself(*searchtri);
                    return ONEDGE;
                }
                if (orgorient == 0.0)
                {
                    lnextself(*searchtri);
                    return ONEDGE;
                }
                return INTRIANGLE;
            }
        }

        /* Move to another triangle.  Leave a trace `backtracktri' in case */
        /*   floating-point roundoff or some such bogey causes us to walk  */
        /*   off a boundary of the triangulation.                          */
        if (moveleft)
        {
            lprev(*searchtri, backtracktri);
            fdest = fapex;
        }
        else
        {
            lnext(*searchtri, backtracktri);
            forg = fapex;
        }
        sym(backtracktri, *searchtri);

        if (m->checksegments && stopatsubsegment)
        {
            /* Check for walking through a subsegment. */
            tspivot(backtracktri, checkedge);
            if (checkedge.ss != m->dummysub)
            {
                /* Go back to the last triangle. */
                otricopy(backtracktri, *searchtri);
                return OUTSIDE;
            }
        }
        /* Check for walking right out of the triangulation. */
        if (searchtri->tri == m->dummytri)
        {
            /* Go back to the last triangle. */
            otricopy(backtracktri, *searchtri);
            return OUTSIDE;
        }

        apex(*searchtri, fapex);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  locate()   Find a triangle or edge containing a given point.             */
/*                                                                           */
/*  Searching begins from one of:  the input `searchtri', a recently         */
/*  encountered triangle `recenttri', or from a triangle chosen from a       */
/*  random sample.  The choice is made by determining which triangle's       */
/*  origin is closest to the point we are searching for.  Normally,          */
/*  `searchtri' should be a handle on the convex hull of the triangulation.  */
/*                                                                           */
/*  Details on the random sampling method can be found in the Mucke, Saias,  */
/*  and Zhu paper cited in the header of this code.                          */
/*                                                                           */
/*  On completion, `searchtri' is a triangle that contains `searchpoint'.    */
/*                                                                           */
/*  Returns ONVERTEX if the point lies on an existing vertex.  `searchtri'   */
/*  is a handle whose origin is the existing vertex.                         */
/*                                                                           */
/*  Returns ONEDGE if the point lies on a mesh edge.  `searchtri' is a       */
/*  handle whose primary edge is the edge on which the point lies.           */
/*                                                                           */
/*  Returns INTRIANGLE if the point lies strictly within a triangle.         */
/*  `searchtri' is a handle on the triangle that contains the point.         */
/*                                                                           */
/*  Returns OUTSIDE if the point lies outside the mesh.  `searchtri' is a    */
/*  handle whose primary edge the point is to the right of.  This might      */
/*  occur when the circumcenter of a triangle falls just slightly outside    */
/*  the mesh due to floating-point roundoff error.  It also occurs when      */
/*  seeking a hole or region point that a foolish user has placed outside    */
/*  the mesh.                                                                */
/*                                                                           */
/*  WARNING:  This routine is designed for convex triangulations, and will   */
/*  not generally work after the holes and concavities have been carved.     */
/*                                                                           */
/*****************************************************************************/

enum locateresult DelaunayTriangle::locate(struct mesh *m,
                         struct behavior *b,
                         vertex searchpoint,
                         struct otri *searchtri)
{
    void **sampleblock;
    char *firsttri;
    struct otri sampletri;
    vertex torg, tdest;
    unsigned long alignptr;
    double searchdist, dist;
    double ahead;
    long samplesperblock, totalsamplesleft, samplesleft;
    long population, totalpopulation;
    triangle ptr; /* Temporary variable used by sym(). */

    /* Record the distance from the suggested starting triangle to the */
    /*   point we seek.                                                */
    org(*searchtri, torg);
    searchdist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
                 (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);

    /* If a recently encountered triangle has been recorded and has not been */
    /*   deallocated, test it as a good starting point.                      */
    if (m->recenttri.tri != (triangle *)NULL)
    {
        if (!deadtri(m->recenttri.tri))
        {
            org(m->recenttri, torg);
            if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1]))
            {
                otricopy(m->recenttri, *searchtri);
                return ONVERTEX;
            }
            dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
                   (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
            if (dist < searchdist)
            {
                otricopy(m->recenttri, *searchtri);
                searchdist = dist;
            }
        }
    }

    /* The number of random samples taken is proportional to the cube root of */
    /*   the number of triangles in the mesh.  The next bit of code assumes   */
    /*   that the number of triangles increases monotonically (or at least    */
    /*   doesn't decrease enough to matter).                                  */
    while (SAMPLEFACTOR * m->samples * m->samples * m->samples <
           m->triangles.items)
    {
        m->samples++;
    }

    /* We'll draw ceiling(samples * TRIPERBLOCK / maxitems) random samples  */
    /*   from each block of triangles (except the first)--until we meet the */
    /*   sample quota.  The ceiling means that blocks at the end might be   */
    /*   neglected, but I don't care.                                       */
    samplesperblock =
        (m->samples * TRIPERBLOCK - 1) / m->triangles.maxitems + 1;
    /* We'll draw ceiling(samples * itemsfirstblock / maxitems) random samples
     */
    /*   from the first block of triangles. */
    samplesleft = (m->samples * m->triangles.itemsfirstblock - 1) /
                      m->triangles.maxitems +
                  1;
    totalsamplesleft = m->samples;
    population       = m->triangles.itemsfirstblock;
    totalpopulation  = m->triangles.maxitems;
    sampleblock      = m->triangles.firstblock;
    sampletri.orient = 0;
    while (totalsamplesleft > 0)
    {
        /* If we're in the last block, `population' needs to be corrected. */
        if (population > totalpopulation)
        {
            population = totalpopulation;
        }
        /* Find a pointer to the first triangle in the block. */
        alignptr = (unsigned long)(sampleblock + 1);
        firsttri =
            (char *)(alignptr + (unsigned long)m->triangles.alignbytes -
                     (alignptr % (unsigned long)m->triangles.alignbytes));

        /* Choose `samplesleft' randomly sampled triangles in this block. */
        do
        {
            sampletri.tri =
                (triangle *)(firsttri +
                             (randomnation((unsigned int)population) *
                              m->triangles.itembytes));
            if (!deadtri(sampletri.tri))
            {
                org(sampletri, torg);
                dist = (searchpoint[0] - torg[0]) * (searchpoint[0] - torg[0]) +
                       (searchpoint[1] - torg[1]) * (searchpoint[1] - torg[1]);
                if (dist < searchdist)
                {
                    otricopy(sampletri, *searchtri);
                    searchdist = dist;
                }
            }

            samplesleft--;
            totalsamplesleft--;
        } while ((samplesleft > 0) && (totalsamplesleft > 0));

        if (totalsamplesleft > 0)
        {
            sampleblock = (void **)*sampleblock;
            samplesleft = samplesperblock;
            totalpopulation -= population;
            population = TRIPERBLOCK;
        }
    }

    /* Where are we? */
    org(*searchtri, torg);
    dest(*searchtri, tdest);
    /* Check the starting triangle's vertices. */
    if ((torg[0] == searchpoint[0]) && (torg[1] == searchpoint[1]))
    {
        return ONVERTEX;
    }
    if ((tdest[0] == searchpoint[0]) && (tdest[1] == searchpoint[1]))
    {
        lnextself(*searchtri);
        return ONVERTEX;
    }
    /* Orient `searchtri' to fit the preconditions of calling preciselocate().
     */
    ahead = counterclockwise(m, b, torg, tdest, searchpoint);
    if (ahead < 0.0)
    {
        /* Turn around so that `searchpoint' is to the left of the */
        /*   edge specified by `searchtri'.                        */
        symself(*searchtri);
    }
    else if (ahead == 0.0)
    {
        /* Check if `searchpoint' is between `torg' and `tdest'. */
        if (((torg[0] < searchpoint[0]) == (searchpoint[0] < tdest[0])) &&
            ((torg[1] < searchpoint[1]) == (searchpoint[1] < tdest[1])))
        {
            return ONEDGE;
        }
    }
    return preciselocate(m, b, searchpoint, searchtri, 0);
}

/**                                                                         **/
/**                                                                         **/
/********* Point location routines end here                          *********/

/********* Mesh transformation routines begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  insertsubseg()   Create a new subsegment and insert it between two       */
/*                   triangles.                                              */
/*                                                                           */
/*  The new subsegment is inserted at the edge described by the handle       */
/*  `tri'.  Its vertices are properly initialized.  The marker `subsegmark'  */
/*  is applied to the subsegment and, if appropriate, its vertices.          */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::insertsubseg(struct mesh *m,
                  struct behavior *b,
                  struct otri *tri,
                  int subsegmark)
{
    struct otri oppotri;
    struct osub newsubseg;
    vertex triorg, tridest;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    org(*tri, triorg);
    dest(*tri, tridest);
    /* Mark vertices if possible. */
    if (vertexmark(triorg) == 0)
    {
        setvertexmark(triorg, subsegmark);
    }
    if (vertexmark(tridest) == 0)
    {
        setvertexmark(tridest, subsegmark);
    }
    /* Check if there's already a subsegment here. */
    tspivot(*tri, newsubseg);
    if (newsubseg.ss == m->dummysub)
    {
        /* Make new subsegment and initialize its vertices. */
        makesubseg(m, &newsubseg);
        setsorg(newsubseg, tridest);
        setsdest(newsubseg, triorg);
        setsegorg(newsubseg, tridest);
        setsegdest(newsubseg, triorg);
        /* Bond new subsegment to the two triangles it is sandwiched between. */
        /*   Note that the facing triangle `oppotri' might be equal to        */
        /*   `dummytri' (outer space), but the new subsegment is bonded to it */
        /*   all the same.                                                    */
        tsbond(*tri, newsubseg);
        sym(*tri, oppotri);
        ssymself(newsubseg);
        tsbond(oppotri, newsubseg);
        setmark(newsubseg, subsegmark);
    }
    else
    {
        if (mark(newsubseg) == 0)
        {
            setmark(newsubseg, subsegmark);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  Terminology                                                              */
/*                                                                           */
/*  A "local transformation" replaces a small set of triangles with another  */
/*  set of triangles.  This may or may not involve inserting or deleting a   */
/*  vertex.                                                                  */
/*                                                                           */
/*  The term "casing" is used to describe the set of triangles that are      */
/*  attached to the triangles being transformed, but are not transformed     */
/*  themselves.  Think of the casing as a fixed hollow structure inside      */
/*  which all the action happens.  A "casing" is only defined relative to    */
/*  a single transformation; each occurrence of a transformation will        */
/*  involve a different casing.                                              */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  flip()   Transform two triangles to two different triangles by flipping  */
/*           an edge counterclockwise within a quadrilateral.                */
/*                                                                           */
/*  Imagine the original triangles, abc and bad, oriented so that the        */
/*  shared edge ab lies in a horizontal plane, with the vertex b on the left */
/*  and the vertex a on the right.  The vertex c lies below the edge, and    */
/*  the vertex d lies above the edge.  The `flipedge' handle holds the edge  */
/*  ab of triangle abc, and is directed left, from vertex a to vertex b.     */
/*                                                                           */
/*  The triangles abc and bad are deleted and replaced by the triangles cdb  */
/*  and dca.  The triangles that represent abc and bad are NOT deallocated;  */
/*  they are reused for dca and cdb, respectively.  Hence, any handles that  */
/*  may have held the original triangles are still valid, although not       */
/*  directed as they were before.                                            */
/*                                                                           */
/*  Upon completion of this routine, the `flipedge' handle holds the edge    */
/*  dc of triangle dca, and is directed down, from vertex d to vertex c.     */
/*  (Hence, the two triangles have rotated counterclockwise.)                */
/*                                                                           */
/*  WARNING:  This transformation is geometrically valid only if the         */
/*  quadrilateral adbc is convex.  Furthermore, this transformation is       */
/*  valid only if there is not a subsegment between the triangles abc and    */
/*  bad.  This routine does not check either of these preconditions, and     */
/*  it is the responsibility of the calling routine to ensure that they are  */
/*  met.  If they are not, the streets shall be filled with wailing and      */
/*  gnashing of teeth.                                                       */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::flip(struct mesh *m, struct behavior *b, struct otri *flipedge)
{
    struct otri botleft, botright;
    struct otri topleft, topright;
    struct otri top;
    struct otri botlcasing, botrcasing;
    struct otri toplcasing, toprcasing;
    struct osub botlsubseg, botrsubseg;
    struct osub toplsubseg, toprsubseg;
    vertex leftvertex, rightvertex, botvertex;
    vertex farvertex;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Identify the vertices of the quadrilateral. */
    org(*flipedge, rightvertex);
    dest(*flipedge, leftvertex);
    apex(*flipedge, botvertex);
    sym(*flipedge, top);
    apex(top, farvertex);

    /* Identify the casing of the quadrilateral. */
    lprev(top, topleft);
    sym(topleft, toplcasing);
    lnext(top, topright);
    sym(topright, toprcasing);
    lnext(*flipedge, botleft);
    sym(botleft, botlcasing);
    lprev(*flipedge, botright);
    sym(botright, botrcasing);
    /* Rotate the quadrilateral one-quarter turn counterclockwise. */
    bond(topleft, botlcasing);
    bond(botleft, botrcasing);
    bond(botright, toprcasing);
    bond(topright, toplcasing);

    if (m->checksegments)
    {
        /* Check for subsegments and rebond them to the quadrilateral. */
        tspivot(topleft, toplsubseg);
        tspivot(botleft, botlsubseg);
        tspivot(botright, botrsubseg);
        tspivot(topright, toprsubseg);
        if (toplsubseg.ss == m->dummysub)
        {
            tsdissolve(topright);
        }
        else
        {
            tsbond(topright, toplsubseg);
        }
        if (botlsubseg.ss == m->dummysub)
        {
            tsdissolve(topleft);
        }
        else
        {
            tsbond(topleft, botlsubseg);
        }
        if (botrsubseg.ss == m->dummysub)
        {
            tsdissolve(botleft);
        }
        else
        {
            tsbond(botleft, botrsubseg);
        }
        if (toprsubseg.ss == m->dummysub)
        {
            tsdissolve(botright);
        }
        else
        {
            tsbond(botright, toprsubseg);
        }
    }

    /* New vertex assignments for the rotated quadrilateral. */
    setorg(*flipedge, farvertex);
    setdest(*flipedge, botvertex);
    setapex(*flipedge, rightvertex);
    setorg(top, botvertex);
    setdest(top, farvertex);
    setapex(top, leftvertex);
}

/*****************************************************************************/
/*                                                                           */
/*  unflip()   Transform two triangles to two different triangles by         */
/*             flipping an edge clockwise within a quadrilateral.  Reverses  */
/*             the flip() operation so that the data structures representing */
/*             the triangles are back where they were before the flip().     */
/*                                                                           */
/*  Imagine the original triangles, abc and bad, oriented so that the        */
/*  shared edge ab lies in a horizontal plane, with the vertex b on the left */
/*  and the vertex a on the right.  The vertex c lies below the edge, and    */
/*  the vertex d lies above the edge.  The `flipedge' handle holds the edge  */
/*  ab of triangle abc, and is directed left, from vertex a to vertex b.     */
/*                                                                           */
/*  The triangles abc and bad are deleted and replaced by the triangles cdb  */
/*  and dca.  The triangles that represent abc and bad are NOT deallocated;  */
/*  they are reused for cdb and dca, respectively.  Hence, any handles that  */
/*  may have held the original triangles are still valid, although not       */
/*  directed as they were before.                                            */
/*                                                                           */
/*  Upon completion of this routine, the `flipedge' handle holds the edge    */
/*  cd of triangle cdb, and is directed up, from vertex c to vertex d.       */
/*  (Hence, the two triangles have rotated clockwise.)                       */
/*                                                                           */
/*  WARNING:  This transformation is geometrically valid only if the         */
/*  quadrilateral adbc is convex.  Furthermore, this transformation is       */
/*  valid only if there is not a subsegment between the triangles abc and    */
/*  bad.  This routine does not check either of these preconditions, and     */
/*  it is the responsibility of the calling routine to ensure that they are  */
/*  met.  If they are not, the streets shall be filled with wailing and      */
/*  gnashing of teeth.                                                       */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::unflip(struct mesh *m, struct behavior *b, struct otri *flipedge)
{
    struct otri botleft, botright;
    struct otri topleft, topright;
    struct otri top;
    struct otri botlcasing, botrcasing;
    struct otri toplcasing, toprcasing;
    struct osub botlsubseg, botrsubseg;
    struct osub toplsubseg, toprsubseg;
    vertex leftvertex, rightvertex, botvertex;
    vertex farvertex;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Identify the vertices of the quadrilateral. */
    org(*flipedge, rightvertex);
    dest(*flipedge, leftvertex);
    apex(*flipedge, botvertex);
    sym(*flipedge, top);
    apex(top, farvertex);

    /* Identify the casing of the quadrilateral. */
    lprev(top, topleft);
    sym(topleft, toplcasing);
    lnext(top, topright);
    sym(topright, toprcasing);
    lnext(*flipedge, botleft);
    sym(botleft, botlcasing);
    lprev(*flipedge, botright);
    sym(botright, botrcasing);
    /* Rotate the quadrilateral one-quarter turn clockwise. */
    bond(topleft, toprcasing);
    bond(botleft, toplcasing);
    bond(botright, botlcasing);
    bond(topright, botrcasing);

    if (m->checksegments)
    {
        /* Check for subsegments and rebond them to the quadrilateral. */
        tspivot(topleft, toplsubseg);
        tspivot(botleft, botlsubseg);
        tspivot(botright, botrsubseg);
        tspivot(topright, toprsubseg);
        if (toplsubseg.ss == m->dummysub)
        {
            tsdissolve(botleft);
        }
        else
        {
            tsbond(botleft, toplsubseg);
        }
        if (botlsubseg.ss == m->dummysub)
        {
            tsdissolve(botright);
        }
        else
        {
            tsbond(botright, botlsubseg);
        }
        if (botrsubseg.ss == m->dummysub)
        {
            tsdissolve(topright);
        }
        else
        {
            tsbond(topright, botrsubseg);
        }
        if (toprsubseg.ss == m->dummysub)
        {
            tsdissolve(topleft);
        }
        else
        {
            tsbond(topleft, toprsubseg);
        }
    }

    /* New vertex assignments for the rotated quadrilateral. */
    setorg(*flipedge, botvertex);
    setdest(*flipedge, farvertex);
    setapex(*flipedge, leftvertex);
    setorg(top, farvertex);
    setdest(top, botvertex);
    setapex(top, rightvertex);
}

/*****************************************************************************/
/*                                                                           */
/*  insertvertex()   Insert a vertex into a Delaunay triangulation,          */
/*                   performing flips as necessary to maintain the Delaunay  */
/*                   property.                                               */
/*                                                                           */
/*  The point `insertvertex' is located.  If `searchtri.tri' is not NULL,    */
/*  the search for the containing triangle begins from `searchtri'.  If      */
/*  `searchtri.tri' is NULL, a full point location procedure is called.      */
/*  If `insertvertex' is found inside a triangle, the triangle is split into */
/*  three; if `insertvertex' lies on an edge, the edge is split in two,      */
/*  thereby splitting the two adjacent triangles into four.  Edge flips are  */
/*  used to restore the Delaunay property.  If `insertvertex' lies on an     */
/*  existing vertex, no action is taken, and the value DUPLICATEVERTEX is    */
/*  returned.  On return, `searchtri' is set to a handle whose origin is the */
/*  existing vertex.                                                         */
/*                                                                           */
/*  Normally, the parameter `splitseg' is set to NULL, implying that no      */
/*  subsegment should be split.  In this case, if `insertvertex' is found to */
/*  lie on a segment, no action is taken, and the value VIOLATINGVERTEX is   */
/*  returned.  On return, `searchtri' is set to a handle whose primary edge  */
/*  is the violated subsegment.                                              */
/*                                                                           */
/*  If the calling routine wishes to split a subsegment by inserting a       */
/*  vertex in it, the parameter `splitseg' should be that subsegment.  In    */
/*  this case, `searchtri' MUST be the triangle handle reached by pivoting   */
/*  from that subsegment; no point location is done.                         */
/*                                                                           */
/*  `segmentflaws' and `triflaws' are flags that indicate whether or not     */
/*  there should be checks for the creation of encroached subsegments or bad */
/*  quality triangles.  If a newly inserted vertex encroaches upon           */
/*  subsegments, these subsegments are added to the list of subsegments to   */
/*  be split if `segmentflaws' is set.  If bad triangles are created, these  */
/*  are added to the queue if `triflaws' is set.                             */
/*                                                                           */
/*  If a duplicate vertex or violated segment does not prevent the vertex    */
/*  from being inserted, the return value will be ENCROACHINGVERTEX if the   */
/*  vertex encroaches upon a subsegment (and checking is enabled), or        */
/*  SUCCESSFULVERTEX otherwise.  In either case, `searchtri' is set to a     */
/*  handle whose origin is the newly inserted vertex.                        */
/*                                                                           */
/*  insertvertex() does not use flip() for reasons of speed; some            */
/*  information can be reused from edge flip to edge flip, like the          */
/*  locations of subsegments.                                                */
/*                                                                           */
/*****************************************************************************/

enum insertvertexresult DelaunayTriangle::insertvertex(struct mesh *m,
                                     struct behavior *b,
                                     vertex newvertex,
                                     struct otri *searchtri,
                                     struct osub *splitseg,
                                     int segmentflaws,
                                     int triflaws)
{
    struct otri horiz;
    struct otri top;
    struct otri botleft, botright;
    struct otri topleft, topright;
    struct otri newbotleft, newbotright;
    struct otri newtopright;
    struct otri botlcasing, botrcasing;
    struct otri toplcasing, toprcasing;
    struct otri testtri;
    struct osub botlsubseg, botrsubseg;
    struct osub toplsubseg, toprsubseg;
    struct osub brokensubseg;
    struct osub checksubseg;
    struct osub rightsubseg;
    struct osub newsubseg;
    struct badsubseg *encroached;
    struct flipstacker *newflip;
    vertex first;
    vertex leftvertex, rightvertex, botvertex, topvertex, farvertex;
    vertex segmentorg, segmentdest;
    double attrib;
    enum insertvertexresult success;
    enum locateresult intersect;
    int doflip;
    int mirrorflag;
    int enq;
    int i;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by spivot() and tspivot(). */

    if (splitseg == (struct osub *)NULL)
    {
        /* Find the location of the vertex to be inserted.  Check if a good */
        /*   starting triangle has already been provided by the caller.     */
        if (searchtri->tri == m->dummytri)
        {
            /* Find a boundary triangle. */
            horiz.tri    = m->dummytri;
            horiz.orient = 0;
            symself(horiz);
            /* Search for a triangle containing `newvertex'. */
            intersect = locate(m, b, newvertex, &horiz);
        }
        else
        {
            /* Start searching from the triangle provided by the caller. */
            otricopy(*searchtri, horiz);
            intersect = preciselocate(m, b, newvertex, &horiz, 1);
        }
    }
    else
    {
        /* The calling routine provides the subsegment in which */
        /*   the vertex is inserted.                             */
        otricopy(*searchtri, horiz);
        intersect = ONEDGE;
    }

    if (intersect == ONVERTEX)
    {
        /* There's already a vertex there.  Return in `searchtri' a triangle */
        /*   whose origin is the existing vertex.                            */
        otricopy(horiz, *searchtri);
        otricopy(horiz, m->recenttri);
        return DUPLICATEVERTEX;
    }
    if ((intersect == ONEDGE) || (intersect == OUTSIDE))
    {
        /* The vertex falls on an edge or boundary. */
        if (m->checksegments && (splitseg == (struct osub *)NULL))
        {
            /* Check whether the vertex falls on a subsegment. */
            tspivot(horiz, brokensubseg);
            if (brokensubseg.ss != m->dummysub)
            {
                /* The vertex falls on a subsegment, and hence will not be
                 * inserted. */
                if (segmentflaws)
                {
                    enq = b->nobisect != 2;
                    if (enq && (b->nobisect == 1))
                    {
                        /* This subsegment may be split only if it is an */
                        /*   internal boundary.                          */
                        sym(horiz, testtri);
                        enq = testtri.tri != m->dummytri;
                    }
                    if (enq)
                    {
                        /* Add the subsegment to the list of encroached
                         * subsegments. */
                        encroached =
                            (struct badsubseg *)poolalloc(&m->badsubsegs);
                        encroached->encsubseg = sencode(brokensubseg);
                        sorg(brokensubseg, encroached->subsegorg);
                        sdest(brokensubseg, encroached->subsegdest);

                    }
                }
                /* Return a handle whose primary edge contains the vertex, */
                /*   which has not been inserted.                          */
                otricopy(horiz, *searchtri);
                otricopy(horiz, m->recenttri);
                return VIOLATINGVERTEX;
            }
        }

        /* Insert the vertex on an edge, dividing one triangle into two (if */
        /*   the edge lies on a boundary) or two triangles into four.       */
        lprev(horiz, botright);
        sym(botright, botrcasing);
        sym(horiz, topright);
        /* Is there a second triangle?  (Or does this edge lie on a boundary?)
         */
        mirrorflag = topright.tri != m->dummytri;
        if (mirrorflag)
        {
            lnextself(topright);
            sym(topright, toprcasing);
            maketriangle(m, b, &newtopright);
        }
        else
        {
            /* Splitting a boundary edge increases the number of boundary edges.
             */
            m->hullsize++;
        }
        maketriangle(m, b, &newbotright);

        /* Set the vertices of changed and new triangles. */
        org(horiz, rightvertex);
        dest(horiz, leftvertex);
        apex(horiz, botvertex);
        setorg(newbotright, botvertex);
        setdest(newbotright, rightvertex);
        setapex(newbotright, newvertex);
        setorg(horiz, newvertex);
        for (i = 0; i < m->eextras; i++)
        {
            /* Set the element attributes of a new triangle. */
            setelemattribute(newbotright, i, elemattribute(botright, i));
        }

        if (mirrorflag)
        {
            dest(topright, topvertex);
            setorg(newtopright, rightvertex);
            setdest(newtopright, topvertex);
            setapex(newtopright, newvertex);
            setorg(topright, newvertex);
            for (i = 0; i < m->eextras; i++)
            {
                /* Set the element attributes of another new triangle. */
                setelemattribute(newtopright, i, elemattribute(topright, i));
            }
        }

        /* There may be subsegments that need to be bonded */
        /*   to the new triangle(s).                       */
        if (m->checksegments)
        {
            tspivot(botright, botrsubseg);
            if (botrsubseg.ss != m->dummysub)
            {
                tsdissolve(botright);
                tsbond(newbotright, botrsubseg);
            }
            if (mirrorflag)
            {
                tspivot(topright, toprsubseg);
                if (toprsubseg.ss != m->dummysub)
                {
                    tsdissolve(topright);
                    tsbond(newtopright, toprsubseg);
                }
            }
        }

        /* Bond the new triangle(s) to the surrounding triangles. */
        bond(newbotright, botrcasing);
        lprevself(newbotright);
        bond(newbotright, botright);
        lprevself(newbotright);
        if (mirrorflag)
        {
            bond(newtopright, toprcasing);
            lnextself(newtopright);
            bond(newtopright, topright);
            lnextself(newtopright);
            bond(newtopright, newbotright);
        }

        if (splitseg != (struct osub *)NULL)
        {
            /* Split the subsegment into two. */
            setsdest(*splitseg, newvertex);
            segorg(*splitseg, segmentorg);
            segdest(*splitseg, segmentdest);
            ssymself(*splitseg);
            spivot(*splitseg, rightsubseg);
            insertsubseg(m, b, &newbotright, mark(*splitseg));
            tspivot(newbotright, newsubseg);
            setsegorg(newsubseg, segmentorg);
            setsegdest(newsubseg, segmentdest);
            sbond(*splitseg, newsubseg);
            ssymself(newsubseg);
            sbond(newsubseg, rightsubseg);
            ssymself(*splitseg);
            /* Transfer the subsegment's boundary marker to the vertex */
            /*   if required.                                          */
            if (vertexmark(newvertex) == 0)
            {
                setvertexmark(newvertex, mark(*splitseg));
            }
        }

        if (m->checkquality)
        {
            poolrestart(&m->flipstackers);
            m->lastflip = (struct flipstacker *)poolalloc(&m->flipstackers);
            m->lastflip->flippedtri = encode(horiz);
            m->lastflip->prevflip   = (struct flipstacker *)1;
        }

        /* Position `horiz' on the first edge to check for */
        /*   the Delaunay property.                        */
        lnextself(horiz);
    }
    else
    {
        /* Insert the vertex in a triangle, splitting it into three. */
        lnext(horiz, botleft);
        lprev(horiz, botright);
        sym(botleft, botlcasing);
        sym(botright, botrcasing);
        maketriangle(m, b, &newbotleft);
        maketriangle(m, b, &newbotright);

        /* Set the vertices of changed and new triangles. */
        org(horiz, rightvertex);
        dest(horiz, leftvertex);
        apex(horiz, botvertex);
        setorg(newbotleft, leftvertex);
        setdest(newbotleft, botvertex);
        setapex(newbotleft, newvertex);
        setorg(newbotright, botvertex);
        setdest(newbotright, rightvertex);
        setapex(newbotright, newvertex);
        setapex(horiz, newvertex);
        for (i = 0; i < m->eextras; i++)
        {
            /* Set the element attributes of the new triangles. */
            attrib = elemattribute(horiz, i);
            setelemattribute(newbotleft, i, attrib);
            setelemattribute(newbotright, i, attrib);
        }

        /* There may be subsegments that need to be bonded */
        /*   to the new triangles.                         */
        if (m->checksegments)
        {
            tspivot(botleft, botlsubseg);
            if (botlsubseg.ss != m->dummysub)
            {
                tsdissolve(botleft);
                tsbond(newbotleft, botlsubseg);
            }
            tspivot(botright, botrsubseg);
            if (botrsubseg.ss != m->dummysub)
            {
                tsdissolve(botright);
                tsbond(newbotright, botrsubseg);
            }
        }

        /* Bond the new triangles to the surrounding triangles. */
        bond(newbotleft, botlcasing);
        bond(newbotright, botrcasing);
        lnextself(newbotleft);
        lprevself(newbotright);
        bond(newbotleft, newbotright);
        lnextself(newbotleft);
        bond(botleft, newbotleft);
        lprevself(newbotright);
        bond(botright, newbotright);

        if (m->checkquality)
        {
            poolrestart(&m->flipstackers);
            m->lastflip = (struct flipstacker *)poolalloc(&m->flipstackers);
            m->lastflip->flippedtri = encode(horiz);
            m->lastflip->prevflip   = (struct flipstacker *)NULL;
        }

    }

    /* The insertion is successful by default, unless an encroached */
    /*   subsegment is found.                                       */
    success = SUCCESSFULVERTEX;
    /* Circle around the newly inserted vertex, checking each edge opposite */
    /*   it for the Delaunay property.  Non-Delaunay edges are flipped.     */
    /*   `horiz' is always the edge being checked.  `first' marks where to  */
    /*   stop circling.                                                     */
    org(horiz, first);
    rightvertex = first;
    dest(horiz, leftvertex);
    /* Circle until finished. */
    while (1)
    {
        /* By default, the edge will be flipped. */
        doflip = 1;

        if (m->checksegments)
        {
            /* Check for a subsegment, which cannot be flipped. */
            tspivot(horiz, checksubseg);
            if (checksubseg.ss != m->dummysub)
            {
                /* The edge is a subsegment and cannot be flipped. */
                doflip = 0;
                if (segmentflaws)
                {
                    /* Does the new vertex encroach upon this subsegment? */
                    if (checkseg4encroach(m, b, &checksubseg))
                    {
                        success = ENCROACHINGVERTEX;
                    }
                }
            }
        }

        if (doflip)
        {
            /* Check if the edge is a boundary edge. */
            sym(horiz, top);
            if (top.tri == m->dummytri)
            {
                /* The edge is a boundary edge and cannot be flipped. */
                doflip = 0;
            }
            else
            {
                /* Find the vertex on the other side of the edge. */
                apex(top, farvertex);
                /* In the incremental Delaunay triangulation algorithm, any of
                 */
                /*   `leftvertex', `rightvertex', and `farvertex' could be
                 * vertices */
                /*   of the triangular bounding box.  These vertices must be */
                /*   treated as if they are infinitely distant, even though
                 * their   */
                /*   "coordinates" are not. */
                if ((leftvertex == m->infvertex1) ||
                    (leftvertex == m->infvertex2) ||
                    (leftvertex == m->infvertex3))
                {
                    /* `leftvertex' is infinitely distant.  Check the convexity
                     * of  */
                    /*   the boundary of the triangulation.  'farvertex' might
                     * be   */
                    /*   infinite as well, but trust me, this same condition
                     * should */
                    /*   be applied. */
                    doflip = counterclockwise(
                                 m, b, newvertex, rightvertex, farvertex) > 0.0;
                }
                else if ((rightvertex == m->infvertex1) ||
                         (rightvertex == m->infvertex2) ||
                         (rightvertex == m->infvertex3))
                {
                    /* `rightvertex' is infinitely distant.  Check the convexity
                     * of */
                    /*   the boundary of the triangulation.  'farvertex' might
                     * be   */
                    /*   infinite as well, but trust me, this same condition
                     * should */
                    /*   be applied. */
                    doflip = counterclockwise(
                                 m, b, farvertex, leftvertex, newvertex) > 0.0;
                }
                else if ((farvertex == m->infvertex1) ||
                         (farvertex == m->infvertex2) ||
                         (farvertex == m->infvertex3))
                {
                    /* `farvertex' is infinitely distant and cannot be inside */
                    /*   the circumcircle of the triangle `horiz'.            */
                    doflip = 0;
                }
                else
                {
                    /* Test whether the edge is locally Delaunay. */
                    doflip = incircle(m,
                                      b,
                                      leftvertex,
                                      newvertex,
                                      rightvertex,
                                      farvertex) > 0.0;
                }
                if (doflip)
                {
                    /* We made it!  Flip the edge `horiz' by rotating its
                     * containing */
                    /*   quadrilateral (the two triangles adjacent to `horiz').
                     */
                    /* Identify the casing of the quadrilateral. */
                    lprev(top, topleft);
                    sym(topleft, toplcasing);
                    lnext(top, topright);
                    sym(topright, toprcasing);
                    lnext(horiz, botleft);
                    sym(botleft, botlcasing);
                    lprev(horiz, botright);
                    sym(botright, botrcasing);
                    /* Rotate the quadrilateral one-quarter turn
                     * counterclockwise. */
                    bond(topleft, botlcasing);
                    bond(botleft, botrcasing);
                    bond(botright, toprcasing);
                    bond(topright, toplcasing);
                    if (m->checksegments)
                    {
                        /* Check for subsegments and rebond them to the
                         * quadrilateral. */
                        tspivot(topleft, toplsubseg);
                        tspivot(botleft, botlsubseg);
                        tspivot(botright, botrsubseg);
                        tspivot(topright, toprsubseg);
                        if (toplsubseg.ss == m->dummysub)
                        {
                            tsdissolve(topright);
                        }
                        else
                        {
                            tsbond(topright, toplsubseg);
                        }
                        if (botlsubseg.ss == m->dummysub)
                        {
                            tsdissolve(topleft);
                        }
                        else
                        {
                            tsbond(topleft, botlsubseg);
                        }
                        if (botrsubseg.ss == m->dummysub)
                        {
                            tsdissolve(botleft);
                        }
                        else
                        {
                            tsbond(botleft, botrsubseg);
                        }
                        if (toprsubseg.ss == m->dummysub)
                        {
                            tsdissolve(botright);
                        }
                        else
                        {
                            tsbond(botright, toprsubseg);
                        }
                    }
                    /* New vertex assignments for the rotated quadrilateral. */
                    setorg(horiz, farvertex);
                    setdest(horiz, newvertex);
                    setapex(horiz, rightvertex);
                    setorg(top, newvertex);
                    setdest(top, farvertex);
                    setapex(top, leftvertex);
                    for (i = 0; i < m->eextras; i++)
                    {
                        /* Take the average of the two triangles' attributes. */
                        attrib = 0.5 * (elemattribute(top, i) +
                                        elemattribute(horiz, i));
                        setelemattribute(top, i, attrib);
                        setelemattribute(horiz, i, attrib);
                    }

                    if (m->checkquality)
                    {
                        newflip =
                            (struct flipstacker *)poolalloc(&m->flipstackers);
                        newflip->flippedtri = encode(horiz);
                        newflip->prevflip   = m->lastflip;
                        m->lastflip         = newflip;
                    }

                    /* On the next iterations, consider the two edges that were
                     */
                    /*   exposed (this is, are now visible to the newly inserted
                     */
                    /*   vertex) by the edge flip. */
                    lprevself(horiz);
                    leftvertex = farvertex;
                }
            }
        }
        if (!doflip)
        {
/* The handle `horiz' is accepted as locally Delaunay. */
            if (triflaws)
            {
                /* Check the triangle `horiz' for quality. */
                testtriangle(m, b, &horiz);
            }
            /* Look for the next edge around the newly inserted vertex. */
            lnextself(horiz);
            sym(horiz, testtri);
            /* Check for finishing a complete revolution about the new vertex,
             * or */
            /*   falling outside  of the triangulation.  The latter will happen
             */
            /*   when a vertex is inserted at a boundary. */
            if ((leftvertex == first) || (testtri.tri == m->dummytri))
            {
                /* We're done.  Return a triangle whose origin is the new
                 * vertex. */
                lnext(horiz, *searchtri);
                lnext(horiz, m->recenttri);
                return success;
            }
            /* Finish finding the next edge around the newly inserted vertex. */
            lnext(testtri, horiz);
            rightvertex = leftvertex;
            dest(horiz, leftvertex);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  triangulatepolygon()   Find the Delaunay triangulation of a polygon that */
/*                         has a certain "nice" shape.  This includes the    */
/*                         polygons that result from deletion of a vertex or */
/*                         insertion of a segment.                           */
/*                                                                           */
/*  This is a conceptually difficult routine.  The starting assumption is    */
/*  that we have a polygon with n sides.  n - 1 of these sides are currently */
/*  represented as edges in the mesh.  One side, called the "base", need not */
/*  be.                                                                      */
/*                                                                           */
/*  Inside the polygon is a structure I call a "fan", consisting of n - 1    */
/*  triangles that share a common origin.  For each of these triangles, the  */
/*  edge opposite the origin is one of the sides of the polygon.  The        */
/*  primary edge of each triangle is the edge directed from the origin to    */
/*  the destination; note that this is not the same edge that is a side of   */
/*  the polygon.  `firstedge' is the primary edge of the first triangle.     */
/*  From there, the triangles follow in counterclockwise order about the     */
/*  polygon, until `lastedge', the primary edge of the last triangle.        */
/*  `firstedge' and `lastedge' are probably connected to other triangles     */
/*  beyond the extremes of the fan, but their identity is not important, as  */
/*  long as the fan remains connected to them.                               */
/*                                                                           */
/*  Imagine the polygon oriented so that its base is at the bottom.  This    */
/*  puts `firstedge' on the far right, and `lastedge' on the far left.       */
/*  The right vertex of the base is the destination of `firstedge', and the  */
/*  left vertex of the base is the apex of `lastedge'.                       */
/*                                                                           */
/*  The challenge now is to find the right sequence of edge flips to         */
/*  transform the fan into a Delaunay triangulation of the polygon.  Each    */
/*  edge flip effectively removes one triangle from the fan, committing it   */
/*  to the polygon.  The resulting polygon has one fewer edge.  If `doflip'  */
/*  is set, the final flip will be performed, resulting in a fan of one      */
/*  (useless?) triangle.  If `doflip' is not set, the final flip is not      */
/*  performed, resulting in a fan of two triangles, and an unfinished        */
/*  triangular polygon that is not yet filled out with a single triangle.    */
/*  On completion of the routine, `lastedge' is the last remaining triangle, */
/*  or the leftmost of the last two.                                         */
/*                                                                           */
/*  Although the flips are performed in the order described above, the       */
/*  decisions about what flips to perform are made in precisely the reverse  */
/*  order.  The recursive triangulatepolygon() procedure makes a decision,   */
/*  uses up to two recursive calls to triangulate the "subproblems"          */
/*  (polygons with fewer edges), and then performs an edge flip.             */
/*                                                                           */
/*  The "decision" it makes is which vertex of the polygon should be         */
/*  connected to the base.  This decision is made by testing every possible  */
/*  vertex.  Once the best vertex is found, the two edges that connect this  */
/*  vertex to the base become the bases for two smaller polygons.  These     */
/*  are triangulated recursively.  Unfortunately, this approach can take     */
/*  O(n^2) time not only in the worst case, but in many common cases.  It's  */
/*  rarely a big deal for vertex deletion, where n is rarely larger than     */
/*  ten, but it could be a big deal for segment insertion, especially if     */
/*  there's a lot of long segments that each cut many triangles.  I ought to */
/*  code a faster algorithm some day.                                        */
/*                                                                           */
/*  The `edgecount' parameter is the number of sides of the polygon,         */
/*  including its base.  `triflaws' is a flag that determines whether the    */
/*  new triangles should be tested for quality, and enqueued if they are     */
/*  bad.                                                                     */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::triangulatepolygon(struct mesh *m,
                        struct behavior *b,
                        struct otri *firstedge,
                        struct otri *lastedge,
                        int edgecount,
                        int doflip,
                        int triflaws)
{
    struct otri testtri;
    struct otri besttri;
    struct otri tempedge;
    vertex leftbasevertex, rightbasevertex;
    vertex testvertex;
    vertex bestvertex;
    int bestnumber;
    int i;
    triangle ptr; /* Temporary variable used by sym(), onext(), and oprev(). */

    /* Identify the base vertices. */
    apex(*lastedge, leftbasevertex);
    dest(*firstedge, rightbasevertex);

    /* Find the best vertex to connect the base to. */
    onext(*firstedge, besttri);
    dest(besttri, bestvertex);
    otricopy(besttri, testtri);
    bestnumber = 1;
    for (i = 2; i <= edgecount - 2; i++)
    {
        onextself(testtri);
        dest(testtri, testvertex);
        /* Is this a better vertex? */
        if (incircle(
                m, b, leftbasevertex, rightbasevertex, bestvertex, testvertex) >
            0.0)
        {
            otricopy(testtri, besttri);
            bestvertex = testvertex;
            bestnumber = i;
        }
    }

    if (bestnumber > 1)
    {
        /* Recursively triangulate the smaller polygon on the right. */
        oprev(besttri, tempedge);
        triangulatepolygon(
            m, b, firstedge, &tempedge, bestnumber + 1, 1, triflaws);
    }
    if (bestnumber < edgecount - 2)
    {
        /* Recursively triangulate the smaller polygon on the left. */
        sym(besttri, tempedge);
        triangulatepolygon(
            m, b, &besttri, lastedge, edgecount - bestnumber, 1, triflaws);
        /* Find `besttri' again; it may have been lost to edge flips. */
        sym(tempedge, besttri);
    }
    if (doflip)
    {
        /* Do one final edge flip. */
        flip(m, b, &besttri);
        if (triflaws)
        {
            /* Check the quality of the newly committed triangle. */
            sym(besttri, testtri);
            testtriangle(m, b, &testtri);
        }
    }
    /* Return the base triangle. */
    otricopy(besttri, *lastedge);
}

/*****************************************************************************/
/*                                                                           */
/*  deletevertex()   Delete a vertex from a Delaunay triangulation, ensuring */
/*                   that the triangulation remains Delaunay.                */
/*                                                                           */
/*  The origin of `deltri' is deleted.  The union of the triangles adjacent  */
/*  to this vertex is a polygon, for which the Delaunay triangulation is     */
/*  found.  Two triangles are removed from the mesh.                         */
/*                                                                           */
/*  Only interior vertices that do not lie on segments or boundaries may be  */
/*  deleted.                                                                 */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::deletevertex(struct mesh *m, struct behavior *b, struct otri *deltri)
{
    struct otri countingtri;
    struct otri firstedge, lastedge;
    struct otri deltriright;
    struct otri lefttri, righttri;
    struct otri leftcasing, rightcasing;
    struct osub leftsubseg, rightsubseg;
    vertex delvertex;
    vertex neworg;
    int edgecount;
    triangle ptr; /* Temporary variable used by sym(), onext(), and oprev(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    org(*deltri, delvertex);

    vertexdealloc(m, delvertex);

    /* Count the degree of the vertex being deleted. */
    onext(*deltri, countingtri);
    edgecount = 1;
    while (!otriequal(*deltri, countingtri))
    {
        edgecount++;
        onextself(countingtri);
    }

    if (edgecount > 3)
    {
        /* Triangulate the polygon defined by the union of all triangles */
        /*   adjacent to the vertex being deleted.  Check the quality of */
        /*   the resulting triangles.                                    */
        onext(*deltri, firstedge);
        oprev(*deltri, lastedge);
        triangulatepolygon(
            m, b, &firstedge, &lastedge, edgecount, 0, !b->nobisect);
    }
    /* Splice out two triangles. */
    lprev(*deltri, deltriright);
    dnext(*deltri, lefttri);
    sym(lefttri, leftcasing);
    oprev(deltriright, righttri);
    sym(righttri, rightcasing);
    bond(*deltri, leftcasing);
    bond(deltriright, rightcasing);
    tspivot(lefttri, leftsubseg);
    if (leftsubseg.ss != m->dummysub)
    {
        tsbond(*deltri, leftsubseg);
    }
    tspivot(righttri, rightsubseg);
    if (rightsubseg.ss != m->dummysub)
    {
        tsbond(deltriright, rightsubseg);
    }

    /* Set the new origin of `deltri' and check its quality. */
    org(lefttri, neworg);
    setorg(*deltri, neworg);
    if (!b->nobisect)
    {
        testtriangle(m, b, deltri);
    }

    /* Delete the two spliced-out triangles. */
    triangledealloc(m, lefttri.tri);
    triangledealloc(m, righttri.tri);
}

/*****************************************************************************/
/*                                                                           */
/*  undovertex()   Undo the most recent vertex insertion.                    */
/*                                                                           */
/*  Walks through the list of transformations (flips and a vertex insertion) */
/*  in the reverse of the order in which they were done, and undoes them.    */
/*  The inserted vertex is removed from the triangulation and deallocated.   */
/*  Two triangles (possibly just one) are also deallocated.                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::undovertex(struct mesh *m, struct behavior *b)
{
    struct otri fliptri;
    struct otri botleft, botright, topright;
    struct otri botlcasing, botrcasing, toprcasing;
    struct otri gluetri;
    struct osub botlsubseg, botrsubseg, toprsubseg;
    vertex botvertex, rightvertex;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Walk through the list of transformations (flips and a vertex insertion)
     */
    /*   in the reverse of the order in which they were done, and undo them. */
    while (m->lastflip != (struct flipstacker *)NULL)
    {
        /* Find a triangle involved in the last unreversed transformation. */
        decode(m->lastflip->flippedtri, fliptri);

        /* We are reversing one of three transformations:  a trisection of one
         */
        /*   triangle into three (by inserting a vertex in the triangle), a */
        /*   bisection of two triangles into four (by inserting a vertex in an
         */
        /*   edge), or an edge flip. */
        if (m->lastflip->prevflip == (struct flipstacker *)NULL)
        {
            /* Restore a triangle that was split into three triangles, */
            /*   so it is again one triangle.                          */
            dprev(fliptri, botleft);
            lnextself(botleft);
            onext(fliptri, botright);
            lprevself(botright);
            sym(botleft, botlcasing);
            sym(botright, botrcasing);
            dest(botleft, botvertex);

            setapex(fliptri, botvertex);
            lnextself(fliptri);
            bond(fliptri, botlcasing);
            tspivot(botleft, botlsubseg);
            tsbond(fliptri, botlsubseg);
            lnextself(fliptri);
            bond(fliptri, botrcasing);
            tspivot(botright, botrsubseg);
            tsbond(fliptri, botrsubseg);

            /* Delete the two spliced-out triangles. */
            triangledealloc(m, botleft.tri);
            triangledealloc(m, botright.tri);
        }
        else if (m->lastflip->prevflip == (struct flipstacker *)1)
        {
            /* Restore two triangles that were split into four triangles, */
            /*   so they are again two triangles.                         */
            lprev(fliptri, gluetri);
            sym(gluetri, botright);
            lnextself(botright);
            sym(botright, botrcasing);
            dest(botright, rightvertex);

            setorg(fliptri, rightvertex);
            bond(gluetri, botrcasing);
            tspivot(botright, botrsubseg);
            tsbond(gluetri, botrsubseg);

            /* Delete the spliced-out triangle. */
            triangledealloc(m, botright.tri);

            sym(fliptri, gluetri);
            if (gluetri.tri != m->dummytri)
            {
                lnextself(gluetri);
                dnext(gluetri, topright);
                sym(topright, toprcasing);

                setorg(gluetri, rightvertex);
                bond(gluetri, toprcasing);
                tspivot(topright, toprsubseg);
                tsbond(gluetri, toprsubseg);

                /* Delete the spliced-out triangle. */
                triangledealloc(m, topright.tri);
            }

            /* This is the end of the list, sneakily encoded. */
            m->lastflip->prevflip = (struct flipstacker *)NULL;
        }
        else
        {
            /* Undo an edge flip. */
            unflip(m, b, &fliptri);
        }

        /* Go on and process the next transformation. */
        m->lastflip = m->lastflip->prevflip;
    }
}

/**                                                                         **/
/**                                                                         **/
/********* Mesh transformation routines end here                     *********/

/********* Divide-and-conquer Delaunay triangulation begins here     *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  The divide-and-conquer bounding box                                      */
/*                                                                           */
/*  I originally implemented the divide-and-conquer and incremental Delaunay */
/*  triangulations using the edge-based data structure presented by Guibas   */
/*  and Stolfi.  Switching to a triangle-based data structure doubled the    */
/*  speed.  However, I had to think of a few extra tricks to maintain the    */
/*  elegance of the original algorithms.                                     */
/*                                                                           */
/*  The "bounding box" used by my variant of the divide-and-conquer          */
/*  algorithm uses one triangle for each edge of the convex hull of the      */
/*  triangulation.  These bounding triangles all share a common apical       */
/*  vertex, which is represented by NULL and which represents nothing.       */
/*  The bounding triangles are linked in a circular fan about this NULL      */
/*  vertex, and the edges on the convex hull of the triangulation appear     */
/*  opposite the NULL vertex.  You might find it easiest to imagine that     */
/*  the NULL vertex is a point in 3D space behind the center of the          */
/*  triangulation, and that the bounding triangles form a sort of cone.      */
/*                                                                           */
/*  This bounding box makes it easy to represent degenerate cases.  For      */
/*  instance, the triangulation of two vertices is a single edge.  This edge */
/*  is represented by two bounding box triangles, one on each "side" of the  */
/*  edge.  These triangles are also linked together in a fan about the NULL  */
/*  vertex.                                                                  */
/*                                                                           */
/*  The bounding box also makes it easy to traverse the convex hull, as the  */
/*  divide-and-conquer algorithm needs to do.                                */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  vertexsort()   Sort an array of vertices by x-coordinate, using the      */
/*                 y-coordinate as a secondary key.                          */
/*                                                                           */
/*  Uses quicksort.  Randomized O(n log n) time.  No, I did not make any of  */
/*  the usual quicksort mistakes.                                            */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::vertexsort(vertex *sortarray, unsigned int arraysize)
{
    int left, right;
    int pivot;
    double pivotx, pivoty;
    vertex temp;

    if (arraysize == 2)
    {
        /* Recursive base case. */
        if ((sortarray[0][0] > sortarray[1][0]) ||
            ((sortarray[0][0] == sortarray[1][0]) &&
             (sortarray[0][1] > sortarray[1][1])))
        {
            temp         = sortarray[1];
            sortarray[1] = sortarray[0];
            sortarray[0] = temp;
        }
        return;
    }
    /* Choose a random pivot to split the array. */
    pivot  = (int)randomnation((unsigned int)arraysize);
    pivotx = sortarray[pivot][0];
    pivoty = sortarray[pivot][1];
    /* Split the array. */
    left  = -1;
    right = arraysize;
    while (left < right)
    {
        /* Search for a vertex whose x-coordinate is too large for the left. */
        do
        {
            left++;
        } while ((left <= right) && ((sortarray[left][0] < pivotx) ||
                                     ((sortarray[left][0] == pivotx) &&
                                      (sortarray[left][1] < pivoty))));
        /* Search for a vertex whose x-coordinate is too small for the right. */
        do
        {
            right--;
        } while ((left <= right) && ((sortarray[right][0] > pivotx) ||
                                     ((sortarray[right][0] == pivotx) &&
                                      (sortarray[right][1] > pivoty))));
        if (left < right)
        {
            /* Swap the left and right vertices. */
            temp             = sortarray[left];
            sortarray[left]  = sortarray[right];
            sortarray[right] = temp;
        }
    }
    if (left > 1)
    {
        /* Recursively sort the left subset. */
        vertexsort(sortarray, left);
    }
    if (right < arraysize - 2)
    {
        /* Recursively sort the right subset. */
        vertexsort(&sortarray[right + 1], arraysize - right - 1);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  vertexmedian()   An order statistic algorithm, almost.  Shuffles an      */
/*                   array of vertices so that the first `median' vertices   */
/*                   occur lexicographically before the remaining vertices.  */
/*                                                                           */
/*  Uses the x-coordinate as the primary key if axis == 0; the y-coordinate  */
/*  if axis == 1.  Very similar to the vertexsort() procedure, but runs in   */
/*  randomized linear time.                                                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::vertexmedian(vertex *sortarray, int arraysize, int median, int axis)
{
    int left, right;
    int pivot;
    double pivot1, pivot2;
    vertex temp;

    if (arraysize == 2)
    {
        /* Recursive base case. */
        if ((sortarray[0][axis] > sortarray[1][axis]) ||
            ((sortarray[0][axis] == sortarray[1][axis]) &&
             (sortarray[0][1 - axis] > sortarray[1][1 - axis])))
        {
            temp         = sortarray[1];
            sortarray[1] = sortarray[0];
            sortarray[0] = temp;
        }
        return;
    }
    /* Choose a random pivot to split the array. */
    pivot  = (int)randomnation((unsigned int)arraysize);
    pivot1 = sortarray[pivot][axis];
    pivot2 = sortarray[pivot][1 - axis];
    /* Split the array. */
    left  = -1;
    right = arraysize;
    while (left < right)
    {
        /* Search for a vertex whose x-coordinate is too large for the left. */
        do
        {
            left++;
        } while ((left <= right) && ((sortarray[left][axis] < pivot1) ||
                                     ((sortarray[left][axis] == pivot1) &&
                                      (sortarray[left][1 - axis] < pivot2))));
        /* Search for a vertex whose x-coordinate is too small for the right. */
        do
        {
            right--;
        } while ((left <= right) && ((sortarray[right][axis] > pivot1) ||
                                     ((sortarray[right][axis] == pivot1) &&
                                      (sortarray[right][1 - axis] > pivot2))));
        if (left < right)
        {
            /* Swap the left and right vertices. */
            temp             = sortarray[left];
            sortarray[left]  = sortarray[right];
            sortarray[right] = temp;
        }
    }
    /* Unlike in vertexsort(), at most one of the following */
    /*   conditionals is true.                             */
    if (left > median)
    {
        /* Recursively shuffle the left subset. */
        vertexmedian(sortarray, left, median, axis);
    }
    if (right < median - 1)
    {
        /* Recursively shuffle the right subset. */
        vertexmedian(&sortarray[right + 1],
                     arraysize - right - 1,
                     median - right - 1,
                     axis);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  alternateaxes()   Sorts the vertices as appropriate for the divide-and-  */
/*                    conquer algorithm with alternating cuts.               */
/*                                                                           */
/*  Partitions by x-coordinate if axis == 0; by y-coordinate if axis == 1.   */
/*  For the base case, subsets containing only two or three vertices are     */
/*  always sorted by x-coordinate.                                           */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::alternateaxes(vertex *sortarray, int arraysize, int axis)
{
    int divider;

    divider = arraysize >> 1;
    if (arraysize <= 3)
    {
        /* Recursive base case:  subsets of two or three vertices will be    */
        /*   handled specially, and should always be sorted by x-coordinate. */
        axis = 0;
    }
    /* Partition with a horizontal or vertical cut. */
    vertexmedian(sortarray, arraysize, divider, axis);
    /* Recursively partition the subsets with a cross cut. */
    if (arraysize - divider >= 2)
    {
        if (divider >= 2)
        {
            alternateaxes(sortarray, divider, 1 - axis);
        }
        alternateaxes(&sortarray[divider], arraysize - divider, 1 - axis);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  mergehulls()   Merge two adjacent Delaunay triangulations into a         */
/*                 single Delaunay triangulation.                            */
/*                                                                           */
/*  This is similar to the algorithm given by Guibas and Stolfi, but uses    */
/*  a triangle-based, rather than edge-based, data structure.                */
/*                                                                           */
/*  The algorithm walks up the gap between the two triangulations, knitting  */
/*  them together.  As they are merged, some of their bounding triangles     */
/*  are converted into double triangles of the triangulation.  The procedure */
/*  pulls each hull's bounding triangles apart, then knits them together     */
/*  like the teeth of two gears.  The Delaunay property determines, at each  */
/*  step, whether the next "tooth" is a bounding triangle of the left hull   */
/*  or the right.  When a bounding triangle becomes real, its apex is        */
/*  changed from NULL to a double vertex. */
/*                                                                           */
/*  Only two new triangles need to be allocated.  These become new bounding  */
/*  triangles at the top and bottom of the seam.  They are used to connect   */
/*  the remaining bounding triangles (those that have not been converted     */
/*  into double triangles) into a single fan. */
/*                                                                           */
/*  On entry, `farleft' and `innerleft' are bounding triangles of the left   */
/*  triangulation.  The origin of `farleft' is the leftmost vertex, and      */
/*  the destination of `innerleft' is the rightmost vertex of the            */
/*  triangulation.  Similarly, `innerright' and `farright' are bounding      */
/*  triangles of the right triangulation.  The origin of `innerright' and    */
/*  destination of `farright' are the leftmost and rightmost vertices.       */
/*                                                                           */
/*  On completion, the origin of `farleft' is the leftmost vertex of the     */
/*  merged triangulation, and the destination of `farright' is the rightmost */
/*  vertex.                                                                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::mergehulls(struct mesh *m,
                struct behavior *b,
                struct otri *farleft,
                struct otri *innerleft,
                struct otri *innerright,
                struct otri *farright,
                int axis)
{
    struct otri leftcand, rightcand;
    struct otri baseedge;
    struct otri nextedge;
    struct otri sidecasing, topcasing, outercasing;
    struct otri checkedge;
    vertex innerleftdest;
    vertex innerrightorg;
    vertex innerleftapex, innerrightapex;
    vertex farleftpt, farrightpt;
    vertex farleftapex, farrightapex;
    vertex lowerleft, lowerright;
    vertex upperleft, upperright;
    vertex nextapex;
    vertex checkvertex;
    int changemade;
    int badedge;
    int leftfinished, rightfinished;
    triangle ptr; /* Temporary variable used by sym(). */

    dest(*innerleft, innerleftdest);
    apex(*innerleft, innerleftapex);
    org(*innerright, innerrightorg);
    apex(*innerright, innerrightapex);
    /* Special treatment for horizontal cuts. */
    if (axis == 1)
    {
        org(*farleft, farleftpt);
        apex(*farleft, farleftapex);
        dest(*farright, farrightpt);
        apex(*farright, farrightapex);
        /* The pointers to the extremal vertices are shifted to point to the */
        /*   topmost and bottommost vertex of each hull, rather than the     */
        /*   leftmost and rightmost vertices.                                */
        while (farleftapex[1] < farleftpt[1])
        {
            lnextself(*farleft);
            symself(*farleft);
            farleftpt = farleftapex;
            apex(*farleft, farleftapex);
        }
        sym(*innerleft, checkedge);
        apex(checkedge, checkvertex);
        while (checkvertex[1] > innerleftdest[1])
        {
            lnext(checkedge, *innerleft);
            innerleftapex = innerleftdest;
            innerleftdest = checkvertex;
            sym(*innerleft, checkedge);
            apex(checkedge, checkvertex);
        }
        while (innerrightapex[1] < innerrightorg[1])
        {
            lnextself(*innerright);
            symself(*innerright);
            innerrightorg = innerrightapex;
            apex(*innerright, innerrightapex);
        }
        sym(*farright, checkedge);
        apex(checkedge, checkvertex);
        while (checkvertex[1] > farrightpt[1])
        {
            lnext(checkedge, *farright);
            farrightapex = farrightpt;
            farrightpt   = checkvertex;
            sym(*farright, checkedge);
            apex(checkedge, checkvertex);
        }
    }
    /* Find a line tangent to and below both hulls. */
    do
    {
        changemade = 0;
        /* Make innerleftdest the "bottommost" vertex of the left hull. */
        if (counterclockwise(
                m, b, innerleftdest, innerleftapex, innerrightorg) > 0.0)
        {
            lprevself(*innerleft);
            symself(*innerleft);
            innerleftdest = innerleftapex;
            apex(*innerleft, innerleftapex);
            changemade = 1;
        }
        /* Make innerrightorg the "bottommost" vertex of the right hull. */
        if (counterclockwise(
                m, b, innerrightapex, innerrightorg, innerleftdest) > 0.0)
        {
            lnextself(*innerright);
            symself(*innerright);
            innerrightorg = innerrightapex;
            apex(*innerright, innerrightapex);
            changemade = 1;
        }
    } while (changemade);
    /* Find the two candidates to be the next "gear tooth." */
    sym(*innerleft, leftcand);
    sym(*innerright, rightcand);
    /* Create the bottom new bounding triangle. */
    maketriangle(m, b, &baseedge);
    /* Connect it to the bounding boxes of the left and right triangulations. */
    bond(baseedge, *innerleft);
    lnextself(baseedge);
    bond(baseedge, *innerright);
    lnextself(baseedge);
    setorg(baseedge, innerrightorg);
    setdest(baseedge, innerleftdest);
    /* Apex is intentionally left NULL. */

    /* Fix the extreme triangles if necessary. */
    org(*farleft, farleftpt);
    if (innerleftdest == farleftpt)
    {
        lnext(baseedge, *farleft);
    }
    dest(*farright, farrightpt);
    if (innerrightorg == farrightpt)
    {
        lprev(baseedge, *farright);
    }
    /* The vertices of the current knitting edge. */
    lowerleft  = innerleftdest;
    lowerright = innerrightorg;
    /* The candidate vertices for knitting. */
    apex(leftcand, upperleft);
    apex(rightcand, upperright);
    /* Walk up the gap between the two triangulations, knitting them together.
     */
    while (1)
    {
        /* Have we reached the top?  (This isn't quite the right question, */
        /*   because even though the left triangulation might seem finished now,
         */
        /*   moving up on the right triangulation might reveal a new vertex of
         */
        /*   the left triangulation.  And vice-versa.) */
        leftfinished =
            counterclockwise(m, b, upperleft, lowerleft, lowerright) <= 0.0;
        rightfinished =
            counterclockwise(m, b, upperright, lowerleft, lowerright) <= 0.0;
        if (leftfinished && rightfinished)
        {
            /* Create the top new bounding triangle. */
            maketriangle(m, b, &nextedge);
            setorg(nextedge, lowerleft);
            setdest(nextedge, lowerright);
            /* Apex is intentionally left NULL. */
            /* Connect it to the bounding boxes of the two triangulations. */
            bond(nextedge, baseedge);
            lnextself(nextedge);
            bond(nextedge, rightcand);
            lnextself(nextedge);
            bond(nextedge, leftcand);

            /* Special treatment for horizontal cuts. */
            if (axis == 1)
            {
                org(*farleft, farleftpt);
                apex(*farleft, farleftapex);
                dest(*farright, farrightpt);
                apex(*farright, farrightapex);
                sym(*farleft, checkedge);
                apex(checkedge, checkvertex);
                /* The pointers to the extremal vertices are restored to the  */
                /*   leftmost and rightmost vertices (rather than topmost and */
                /*   bottommost).                                             */
                while (checkvertex[0] < farleftpt[0])
                {
                    lprev(checkedge, *farleft);
                    farleftapex = farleftpt;
                    farleftpt   = checkvertex;
                    sym(*farleft, checkedge);
                    apex(checkedge, checkvertex);
                }
                while (farrightapex[0] > farrightpt[0])
                {
                    lprevself(*farright);
                    symself(*farright);
                    farrightpt = farrightapex;
                    apex(*farright, farrightapex);
                }
            }
            return;
        }
        /* Consider eliminating edges from the left triangulation. */
        if (!leftfinished)
        {
            /* What vertex would be exposed if an edge were deleted? */
            lprev(leftcand, nextedge);
            symself(nextedge);
            apex(nextedge, nextapex);
            /* If nextapex is NULL, then no vertex would be exposed; the */
            /*   triangulation would have been eaten right through.      */
            if (nextapex != (vertex)NULL)
            {
                /* Check whether the edge is Delaunay. */
                badedge =
                    incircle(m, b, lowerleft, lowerright, upperleft, nextapex) >
                    0.0;
                while (badedge)
                {
                    /* Eliminate the edge with an edge flip.  As a result, the
                     */
                    /*   left triangulation will have one more boundary
                     * triangle. */
                    lnextself(nextedge);
                    sym(nextedge, topcasing);
                    lnextself(nextedge);
                    sym(nextedge, sidecasing);
                    bond(nextedge, topcasing);
                    bond(leftcand, sidecasing);
                    lnextself(leftcand);
                    sym(leftcand, outercasing);
                    lprevself(nextedge);
                    bond(nextedge, outercasing);
                    /* Correct the vertices to reflect the edge flip. */
                    setorg(leftcand, lowerleft);
                    setdest(leftcand, NULL);
                    setapex(leftcand, nextapex);
                    setorg(nextedge, NULL);
                    setdest(nextedge, upperleft);
                    setapex(nextedge, nextapex);
                    /* Consider the newly exposed vertex. */
                    upperleft = nextapex;
                    /* What vertex would be exposed if another edge were
                     * deleted? */
                    otricopy(sidecasing, nextedge);
                    apex(nextedge, nextapex);
                    if (nextapex != (vertex)NULL)
                    {
                        /* Check whether the edge is Delaunay. */
                        badedge = incircle(m,
                                           b,
                                           lowerleft,
                                           lowerright,
                                           upperleft,
                                           nextapex) > 0.0;
                    }
                    else
                    {
                        /* Avoid eating right through the triangulation. */
                        badedge = 0;
                    }
                }
            }
        }
        /* Consider eliminating edges from the right triangulation. */
        if (!rightfinished)
        {
            /* What vertex would be exposed if an edge were deleted? */
            lnext(rightcand, nextedge);
            symself(nextedge);
            apex(nextedge, nextapex);
            /* If nextapex is NULL, then no vertex would be exposed; the */
            /*   triangulation would have been eaten right through.      */
            if (nextapex != (vertex)NULL)
            {
                /* Check whether the edge is Delaunay. */
                badedge =
                    incircle(
                        m, b, lowerleft, lowerright, upperright, nextapex) >
                    0.0;
                while (badedge)
                {
                    /* Eliminate the edge with an edge flip.  As a result, the
                     */
                    /*   right triangulation will have one more boundary
                     * triangle. */
                    lprevself(nextedge);
                    sym(nextedge, topcasing);
                    lprevself(nextedge);
                    sym(nextedge, sidecasing);
                    bond(nextedge, topcasing);
                    bond(rightcand, sidecasing);
                    lprevself(rightcand);
                    sym(rightcand, outercasing);
                    lnextself(nextedge);
                    bond(nextedge, outercasing);
                    /* Correct the vertices to reflect the edge flip. */
                    setorg(rightcand, NULL);
                    setdest(rightcand, lowerright);
                    setapex(rightcand, nextapex);
                    setorg(nextedge, upperright);
                    setdest(nextedge, NULL);
                    setapex(nextedge, nextapex);
                    /* Consider the newly exposed vertex. */
                    upperright = nextapex;
                    /* What vertex would be exposed if another edge were
                     * deleted? */
                    otricopy(sidecasing, nextedge);
                    apex(nextedge, nextapex);
                    if (nextapex != (vertex)NULL)
                    {
                        /* Check whether the edge is Delaunay. */
                        badedge = incircle(m,
                                           b,
                                           lowerleft,
                                           lowerright,
                                           upperright,
                                           nextapex) > 0.0;
                    }
                    else
                    {
                        /* Avoid eating right through the triangulation. */
                        badedge = 0;
                    }
                }
            }
        }
        if (leftfinished ||
            (!rightfinished &&
             (incircle(m, b, upperleft, lowerleft, lowerright, upperright) >
              0.0)))
        {
            /* Knit the triangulations, adding an edge from `lowerleft' */
            /*   to `upperright'.                                       */
            bond(baseedge, rightcand);
            lprev(rightcand, baseedge);
            setdest(baseedge, lowerleft);
            lowerright = upperright;
            sym(baseedge, rightcand);
            apex(rightcand, upperright);
        }
        else
        {
            /* Knit the triangulations, adding an edge from `upperleft' */
            /*   to `lowerright'.                                       */
            bond(baseedge, leftcand);
            lnext(leftcand, baseedge);
            setorg(baseedge, lowerright);
            lowerleft = upperleft;
            sym(baseedge, leftcand);
            apex(leftcand, upperleft);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  divconqrecurse()   Recursively form a Delaunay triangulation by the      */
/*                     divide-and-conquer method.                            */
/*                                                                           */
/*  Recursively breaks down the problem into smaller pieces, which are       */
/*  knitted together by mergehulls().  The base cases (problems of two or    */
/*  three vertices) are handled specially here.                              */
/*                                                                           */
/*  On completion, `farleft' and `farright' are bounding triangles such that */
/*  the origin of `farleft' is the leftmost vertex (breaking ties by         */
/*  choosing the highest leftmost vertex), and the destination of            */
/*  `farright' is the rightmost vertex (breaking ties by choosing the        */
/*  lowest rightmost vertex).                                                */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::divconqrecurse(struct mesh *m,
                    struct behavior *b,
                    vertex *sortarray,
                    int vertices,
                    int axis,
                    struct otri *farleft,
                    struct otri *farright)
{
    struct otri midtri, tri1, tri2, tri3;
    struct otri innerleft, innerright;
    double area;
    int divider;

    if (vertices == 2)
    {
        /* The triangulation of two vertices is an edge.  An edge is */
        /*   represented by two bounding triangles.                  */
        maketriangle(m, b, farleft);
        setorg(*farleft, sortarray[0]);
        setdest(*farleft, sortarray[1]);
        /* The apex is intentionally left NULL. */
        maketriangle(m, b, farright);
        setorg(*farright, sortarray[1]);
        setdest(*farright, sortarray[0]);
        /* The apex is intentionally left NULL. */
        bond(*farleft, *farright);
        lprevself(*farleft);
        lnextself(*farright);
        bond(*farleft, *farright);
        lprevself(*farleft);
        lnextself(*farright);
        bond(*farleft, *farright);

        /* Ensure that the origin of `farleft' is sortarray[0]. */
        lprev(*farright, *farleft);
        return;
    }
    else if (vertices == 3)
    {
        /* The triangulation of three vertices is either a triangle (with */
        /*   three bounding triangles) or two edges (with four bounding   */
        /*   triangles).  In either case, four triangles are created.     */
        maketriangle(m, b, &midtri);
        maketriangle(m, b, &tri1);
        maketriangle(m, b, &tri2);
        maketriangle(m, b, &tri3);
        area = counterclockwise(m, b, sortarray[0], sortarray[1], sortarray[2]);
        if (area == 0.0)
        {
            /* Three collinear vertices; the triangulation is two edges. */
            setorg(midtri, sortarray[0]);
            setdest(midtri, sortarray[1]);
            setorg(tri1, sortarray[1]);
            setdest(tri1, sortarray[0]);
            setorg(tri2, sortarray[2]);
            setdest(tri2, sortarray[1]);
            setorg(tri3, sortarray[1]);
            setdest(tri3, sortarray[2]);
            /* All apices are intentionally left NULL. */
            bond(midtri, tri1);
            bond(tri2, tri3);
            lnextself(midtri);
            lprevself(tri1);
            lnextself(tri2);
            lprevself(tri3);
            bond(midtri, tri3);
            bond(tri1, tri2);
            lnextself(midtri);
            lprevself(tri1);
            lnextself(tri2);
            lprevself(tri3);
            bond(midtri, tri1);
            bond(tri2, tri3);
            /* Ensure that the origin of `farleft' is sortarray[0]. */
            otricopy(tri1, *farleft);
            /* Ensure that the destination of `farright' is sortarray[2]. */
            otricopy(tri2, *farright);
        }
        else
        {
            /* The three vertices are not collinear; the triangulation is one */
            /*   triangle, namely `midtri'.                                   */
            setorg(midtri, sortarray[0]);
            setdest(tri1, sortarray[0]);
            setorg(tri3, sortarray[0]);
            /* Apices of tri1, tri2, and tri3 are left NULL. */
            if (area > 0.0)
            {
                /* The vertices are in counterclockwise order. */
                setdest(midtri, sortarray[1]);
                setorg(tri1, sortarray[1]);
                setdest(tri2, sortarray[1]);
                setapex(midtri, sortarray[2]);
                setorg(tri2, sortarray[2]);
                setdest(tri3, sortarray[2]);
            }
            else
            {
                /* The vertices are in clockwise order. */
                setdest(midtri, sortarray[2]);
                setorg(tri1, sortarray[2]);
                setdest(tri2, sortarray[2]);
                setapex(midtri, sortarray[1]);
                setorg(tri2, sortarray[1]);
                setdest(tri3, sortarray[1]);
            }
            /* The topology does not depend on how the vertices are ordered. */
            bond(midtri, tri1);
            lnextself(midtri);
            bond(midtri, tri2);
            lnextself(midtri);
            bond(midtri, tri3);
            lprevself(tri1);
            lnextself(tri2);
            bond(tri1, tri2);
            lprevself(tri1);
            lprevself(tri3);
            bond(tri1, tri3);
            lnextself(tri2);
            lprevself(tri3);
            bond(tri2, tri3);
            /* Ensure that the origin of `farleft' is sortarray[0]. */
            otricopy(tri1, *farleft);
            /* Ensure that the destination of `farright' is sortarray[2]. */
            if (area > 0.0)
            {
                otricopy(tri2, *farright);
            }
            else
            {
                lnext(*farleft, *farright);
            }
        }

        return;
    }
    else
    {
        /* Split the vertices in half. */
        divider = vertices >> 1;
        /* Recursively triangulate each half. */
        divconqrecurse(m, b, sortarray, divider, 1 - axis, farleft, &innerleft);
        divconqrecurse(m,
                       b,
                       &sortarray[divider],
                       vertices - divider,
                       1 - axis,
                       &innerright,
                       farright);

        /* Merge the two triangulations into one. */
        mergehulls(m, b, farleft, &innerleft, &innerright, farright, axis);
    }
}

long DelaunayTriangle::removeghosts(struct mesh *m, struct behavior *b, struct otri *startghost)
{
    struct otri searchedge;
    struct otri dissolveedge;
    struct otri deadtriangle;
    vertex markorg;
    long hullsize;
    triangle ptr; /* Temporary variable used by sym(). */


    /* Find an edge on the convex hull to start point location from. */
    lprev(*startghost, searchedge);
    symself(searchedge);
    m->dummytri[0] = encode(searchedge);
    /* Remove the bounding box and count the convex hull edges. */
    otricopy(*startghost, dissolveedge);
    hullsize = 0;
    do
    {
        hullsize++;
        lnext(dissolveedge, deadtriangle);
        lprevself(dissolveedge);
        symself(dissolveedge);
        /* If no PSLG is involved, set the boundary markers of all the vertices
         */
        /*   on the convex hull.  If a PSLG is used, this step is done later. */
        if (!b->poly)
        {
            /* Watch out for the case where all the input vertices are
             * collinear. */
            if (dissolveedge.tri != m->dummytri)
            {
                org(dissolveedge, markorg);
                if (vertexmark(markorg) == 0)
                {
                    setvertexmark(markorg, 1);
                }
            }
        }
        /* Remove a bounding triangle from a convex hull triangle. */
        dissolve(dissolveedge);
        /* Find the next bounding triangle. */
        sym(deadtriangle, dissolveedge);
        /* Delete the bounding triangle. */
        triangledealloc(m, deadtriangle.tri);
    } while (!otriequal(dissolveedge, *startghost));
    return hullsize;
}

/*****************************************************************************/
/*                                                                           */
/*  divconqdelaunay()   Form a Delaunay triangulation by the divide-and-     */
/*                      conquer method.                                      */
/*                                                                           */
/*  Sorts the vertices, calls a recursive procedure to triangulate them, and */
/*  removes the bounding box, setting boundary markers as appropriate.       */
/*                                                                           */
/*****************************************************************************/

long DelaunayTriangle::divconqdelaunay(struct mesh *m, struct behavior *b)
{
    vertex *sortarray;
    struct otri hullleft, hullright;
    int divider;
    int i, j;

    /* Allocate an array of pointers to vertices for sorting. */
    sortarray = (vertex *)trimalloc(m->invertices * (int)sizeof(vertex));
    traversalinit(&m->vertices);
    for (i = 0; i < m->invertices; i++)
    {
        sortarray[i] = vertextraverse(m);
    }
    /* Sort the vertices. */
    vertexsort(sortarray, m->invertices);
    /* Discard duplicate vertices, which can really mess up the algorithm. */
    i = 0;
    for (j = 1; j < m->invertices; j++)
    {
        if ((sortarray[i][0] == sortarray[j][0]) &&
            (sortarray[i][1] == sortarray[j][1]))
        {
            setvertextype(sortarray[j], UNDEADVERTEX);
            m->undeads++;
        }
        else
        {
            i++;
            sortarray[i] = sortarray[j];
        }
    }
    i++;

    /* Re-sort the array of vertices to accommodate alternating cuts. */
    divider = i >> 1;
    if (i - divider >= 2)
    {
        if (divider >= 2)
        {
            alternateaxes(sortarray, divider, 1);
        }
        alternateaxes(&sortarray[divider], i - divider, 1);
    }

    /* Form the Delaunay triangulation. */
    divconqrecurse(m, b, sortarray, i, 0, &hullleft, &hullright);
    trifree((void *)sortarray);

    return removeghosts(m, b, &hullleft);
}

/**                                                                         **/
/**                                                                         **/
/********* Divide-and-conquer Delaunay triangulation ends here       *********/

/********* General mesh construction routines begin here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  delaunay()   Form a Delaunay triangulation.                              */
/*                                                                           */
/*****************************************************************************/

long DelaunayTriangle::delaunay(struct mesh *m, struct behavior *b)
{
    long hulledges;

    m->eextras = 0;
    initializetrisubpools(m, b);

    hulledges = divconqdelaunay(m, b);

    if (m->triangles.items == 0)
    {
        /* The input vertices were all collinear, so there are no triangles. */
        return 0l;
    }
    else
    {
        return hulledges;
    }
}

/**                                                                         **/
/**                                                                         **/
/********* General mesh construction routines end here               *********/

/********* Segment insertion begins here                             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  finddirection()   Find the first triangle on the path from one point     */
/*                    to another.                                            */
/*                                                                           */
/*  Finds the triangle that intersects a line segment drawn from the         */
/*  origin of `searchtri' to the point `searchpoint', and returns the result */
/*  in `searchtri'.  The origin of `searchtri' does not change, even though  */
/*  the triangle returned may differ from the one passed in.  This routine   */
/*  is used to find the direction to move in to get from one point to        */
/*  another.                                                                 */
/*                                                                           */
/*  The return value notes whether the destination or apex of the found      */
/*  triangle is collinear with the two points in question.                   */
/*                                                                           */
/*****************************************************************************/

enum finddirectionresult DelaunayTriangle::finddirection(struct mesh *m,
                                       struct behavior *b,
                                       struct otri *searchtri,
                                       vertex searchpoint)
{
    struct otri checktri;
    vertex startvertex;
    vertex leftvertex, rightvertex;
    double leftccw, rightccw;
    int leftflag, rightflag;
    triangle ptr; /* Temporary variable used by onext() and oprev(). */

    org(*searchtri, startvertex);
    dest(*searchtri, rightvertex);
    apex(*searchtri, leftvertex);
    /* Is `searchpoint' to the left? */
    leftccw  = counterclockwise(m, b, searchpoint, startvertex, leftvertex);
    leftflag = leftccw > 0.0;
    /* Is `searchpoint' to the right? */
    rightccw  = counterclockwise(m, b, startvertex, searchpoint, rightvertex);
    rightflag = rightccw > 0.0;
    if (leftflag && rightflag)
    {
        /* `searchtri' faces directly away from `searchpoint'.  We could go left
         */
        /*   or right.  Ask whether it's a triangle or a boundary on the left.
         */
        onext(*searchtri, checktri);
        if (checktri.tri == m->dummytri)
        {
            leftflag = 0;
        }
        else
        {
            rightflag = 0;
        }
    }
    while (leftflag)
    {
        /* Turn left until satisfied. */
        onextself(*searchtri);
        if (searchtri->tri == m->dummytri)
        {
            printf("Internal error in finddirection():  Unable to find a\n");
            printf("  triangle leading from (%.12g, %.12g) to",
                   startvertex[0],
                   startvertex[1]);
            printf("  (%.12g, %.12g).\n", searchpoint[0], searchpoint[1]);
            internalerror();
        }
        apex(*searchtri, leftvertex);
        rightccw = leftccw;
        leftccw  = counterclockwise(m, b, searchpoint, startvertex, leftvertex);
        leftflag = leftccw > 0.0;
    }
    while (rightflag)
    {
        /* Turn right until satisfied. */
        oprevself(*searchtri);
        if (searchtri->tri == m->dummytri)
        {
            printf("Internal error in finddirection():  Unable to find a\n");
            printf("  triangle leading from (%.12g, %.12g) to",
                   startvertex[0],
                   startvertex[1]);
            printf("  (%.12g, %.12g).\n", searchpoint[0], searchpoint[1]);
            internalerror();
        }
        dest(*searchtri, rightvertex);
        leftccw = rightccw;
        rightccw =
            counterclockwise(m, b, startvertex, searchpoint, rightvertex);
        rightflag = rightccw > 0.0;
    }
    if (leftccw == 0.0)
    {
        return LEFTCOLLINEAR;
    }
    else if (rightccw == 0.0)
    {
        return RIGHTCOLLINEAR;
    }
    else
    {
        return WITHIN;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  segmentintersection()   Find the intersection of an existing segment     */
/*                          and a segment that is being inserted.  Insert    */
/*                          a vertex at the intersection, splitting an       */
/*                          existing subsegment.                             */
/*                                                                           */
/*  The segment being inserted connects the apex of splittri to endpoint2.   */
/*  splitsubseg is the subsegment being split, and MUST adjoin splittri.     */
/*  Hence, endpoints of the subsegment being split are the origin and        */
/*  destination of splittri.                                                 */
/*                                                                           */
/*  On completion, splittri is a handle having the newly inserted            */
/*  intersection point as its origin, and endpoint1 as its destination.      */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::segmentintersection(struct mesh *m,
                         struct behavior *b,
                         struct otri *splittri,
                         struct osub *splitsubseg,
                         vertex endpoint2)
{
    struct osub opposubseg;
    vertex endpoint1;
    vertex torg, tdest;
    vertex leftvertex, rightvertex;
    vertex newvertex;
    enum insertvertexresult success;
    double ex, ey;
    double tx, ty;
    double etx, ety;
    double split, denom;
    int i;
    triangle ptr; /* Temporary variable used by onext(). */
    subseg sptr;  /* Temporary variable used by snext(). */

    /* Find the other three segment endpoints. */
    apex(*splittri, endpoint1);
    org(*splittri, torg);
    dest(*splittri, tdest);
    /* Segment intersection formulae; see the Antonio reference. */
    tx    = tdest[0] - torg[0];
    ty    = tdest[1] - torg[1];
    ex    = endpoint2[0] - endpoint1[0];
    ey    = endpoint2[1] - endpoint1[1];
    etx   = torg[0] - endpoint2[0];
    ety   = torg[1] - endpoint2[1];
    denom = ty * ex - tx * ey;
    if (denom == 0.0)
    {
        printf("Internal error in segmentintersection():");
        printf("  Attempt to find intersection of parallel segments.\n");
        internalerror();
    }
    split = (ey * etx - ex * ety) / denom;
    /* Create the new vertex. */
    newvertex = (vertex)poolalloc(&m->vertices);
    /* Interpolate its coordinate and attributes. */
    for (i = 0; i < 2 + m->nextras; i++)
    {
        newvertex[i] = torg[i] + split * (tdest[i] - torg[i]);
    }
    setvertexmark(newvertex, mark(*splitsubseg));
    setvertextype(newvertex, INPUTVERTEX);

    /* Insert the intersection vertex.  This should always succeed. */
    success = insertvertex(m, b, newvertex, splittri, splitsubseg, 0, 0);
    if (success != SUCCESSFULVERTEX)
    {
        printf("Internal error in segmentintersection():\n");
        printf("  Failure to split a segment.\n");
        internalerror();
    }
    /* Record a triangle whose origin is the new vertex. */
    setvertex2tri(newvertex, encode(*splittri));
    if (m->steinerleft > 0)
    {
        m->steinerleft--;
    }

    /* Divide the segment into two, and correct the segment endpoints. */
    ssymself(*splitsubseg);
    spivot(*splitsubseg, opposubseg);
    sdissolve(*splitsubseg);
    sdissolve(opposubseg);
    do
    {
        setsegorg(*splitsubseg, newvertex);
        snextself(*splitsubseg);
    } while (splitsubseg->ss != m->dummysub);
    do
    {
        setsegorg(opposubseg, newvertex);
        snextself(opposubseg);
    } while (opposubseg.ss != m->dummysub);

    /* Inserting the vertex may have caused edge flips.  We wish to rediscover
     */
    /*   the edge connecting endpoint1 to the new intersection vertex. */
    finddirection(m, b, splittri, endpoint1);
    dest(*splittri, rightvertex);
    apex(*splittri, leftvertex);
    if ((leftvertex[0] == endpoint1[0]) && (leftvertex[1] == endpoint1[1]))
    {
        onextself(*splittri);
    }
    else if ((rightvertex[0] != endpoint1[0]) ||
             (rightvertex[1] != endpoint1[1]))
    {
        printf("Internal error in segmentintersection():\n");
        printf("  Topological inconsistency after splitting a segment.\n");
        internalerror();
    }
    /* `splittri' should have destination endpoint1. */
}

/*****************************************************************************/
/*                                                                           */
/*  scoutsegment()   Scout the first triangle on the path from one endpoint  */
/*                   to another, and check for completion (reaching the      */
/*                   second endpoint), a collinear vertex, or the            */
/*                   intersection of two segments.                           */
/*                                                                           */
/*  Returns one if the entire segment is successfully inserted, and zero if  */
/*  the job must be finished by conformingedge() or constrainededge().       */
/*                                                                           */
/*  If the first triangle on the path has the second endpoint as its         */
/*  destination or apex, a subsegment is inserted and the job is done.       */
/*                                                                           */
/*  If the first triangle on the path has a destination or apex that lies on */
/*  the segment, a subsegment is inserted connecting the first endpoint to   */
/*  the collinear vertex, and the search is continued from the collinear     */
/*  vertex.                                                                  */
/*                                                                           */
/*  If the first triangle on the path has a subsegment opposite its origin,  */
/*  then there is a segment that intersects the segment being inserted.      */
/*  Their intersection vertex is inserted, splitting the subsegment.         */
/*                                                                           */
/*****************************************************************************/

int DelaunayTriangle::scoutsegment(struct mesh *m,
                 struct behavior *b,
                 struct otri *searchtri,
                 vertex endpoint2,
                 int newmark)
{
    struct otri crosstri;
    struct osub crosssubseg;
    vertex leftvertex, rightvertex;
    enum finddirectionresult collinear;
    subseg sptr; /* Temporary variable used by tspivot(). */

    collinear = finddirection(m, b, searchtri, endpoint2);
    dest(*searchtri, rightvertex);
    apex(*searchtri, leftvertex);
    if (((leftvertex[0] == endpoint2[0]) && (leftvertex[1] == endpoint2[1])) ||
        ((rightvertex[0] == endpoint2[0]) && (rightvertex[1] == endpoint2[1])))
    {
        /* The segment is already an edge in the mesh. */
        if ((leftvertex[0] == endpoint2[0]) && (leftvertex[1] == endpoint2[1]))
        {
            lprevself(*searchtri);
        }
        /* Insert a subsegment, if there isn't already one there. */
        insertsubseg(m, b, searchtri, newmark);
        return 1;
    }
    else if (collinear == LEFTCOLLINEAR)
    {
        /* We've collided with a vertex between the segment's endpoints. */
        /* Make the collinear vertex be the triangle's origin. */
        lprevself(*searchtri);
        insertsubseg(m, b, searchtri, newmark);
        /* Insert the remainder of the segment. */
        return scoutsegment(m, b, searchtri, endpoint2, newmark);
    }
    else if (collinear == RIGHTCOLLINEAR)
    {
        /* We've collided with a vertex between the segment's endpoints. */
        insertsubseg(m, b, searchtri, newmark);
        /* Make the collinear vertex be the triangle's origin. */
        lnextself(*searchtri);
        /* Insert the remainder of the segment. */
        return scoutsegment(m, b, searchtri, endpoint2, newmark);
    }
    else
    {
        lnext(*searchtri, crosstri);
        tspivot(crosstri, crosssubseg);
        /* Check for a crossing segment. */
        if (crosssubseg.ss == m->dummysub)
        {
            return 0;
        }
        else
        {
            /* Insert a vertex at the intersection. */
            segmentintersection(m, b, &crosstri, &crosssubseg, endpoint2);
            otricopy(crosstri, *searchtri);
            insertsubseg(m, b, searchtri, newmark);
            /* Insert the remainder of the segment. */
            return scoutsegment(m, b, searchtri, endpoint2, newmark);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  conformingedge()   Force a segment into a conforming Delaunay            */
/*                     triangulation by inserting a vertex at its midpoint,  */
/*                     and recursively forcing in the two half-segments if   */
/*                     necessary.                                            */
/*                                                                           */
/*  Generates a sequence of subsegments connecting `endpoint1' to            */
/*  `endpoint2'.  `newmark' is the boundary marker of the segment, assigned  */
/*  to each new splitting vertex and subsegment.                             */
/*                                                                           */
/*  Note that conformingedge() does not always maintain the conforming       */
/*  Delaunay property.  Once inserted, segments are locked into place;       */
/*  vertices inserted later (to force other segments in) may render these    */
/*  fixed segments non-Delaunay.  The conforming Delaunay property will be   */
/*  restored by enforcequality() by splitting encroached subsegments.        */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::conformingedge(struct mesh *m,
                    struct behavior *b,
                    vertex endpoint1,
                    vertex endpoint2,
                    int newmark)
{
    struct otri searchtri1, searchtri2;
    struct osub brokensubseg;
    vertex newvertex;
    vertex midvertex1, midvertex2;
    enum insertvertexresult success;
    int i;
    subseg sptr; /* Temporary variable used by tspivot(). */

    /* Create a new vertex to insert in the middle of the segment. */
    newvertex = (vertex)poolalloc(&m->vertices);
    /* Interpolate coordinates and attributes. */
    for (i = 0; i < 2 + m->nextras; i++)
    {
        newvertex[i] = 0.5 * (endpoint1[i] + endpoint2[i]);
    }
    setvertexmark(newvertex, newmark);
    setvertextype(newvertex, SEGMENTVERTEX);
    /* No known triangle to search from. */
    searchtri1.tri = m->dummytri;
    /* Attempt to insert the new vertex. */
    success =
        insertvertex(m, b, newvertex, &searchtri1, (struct osub *)NULL, 0, 0);
    if (success == DUPLICATEVERTEX)
    {

        /* Use the vertex that's already there. */
        vertexdealloc(m, newvertex);
        org(searchtri1, newvertex);
    }
    else
    {
        if (success == VIOLATINGVERTEX)
        {

            /* By fluke, we've landed right on another segment.  Split it. */
            tspivot(searchtri1, brokensubseg);
            success =
                insertvertex(m, b, newvertex, &searchtri1, &brokensubseg, 0, 0);
            if (success != SUCCESSFULVERTEX)
            {
                printf("Internal error in conformingedge():\n");
                printf("  Failure to split a segment.\n");
                internalerror();
            }
        }
        /* The vertex has been inserted successfully. */
        if (m->steinerleft > 0)
        {
            m->steinerleft--;
        }
    }
    otricopy(searchtri1, searchtri2);
    /* `searchtri1' and `searchtri2' are fastened at their origins to         */
    /*   `newvertex', and will be directed toward `endpoint1' and `endpoint2' */
    /*   respectively.  First, we must get `searchtri2' out of the way so it  */
    /*   won't be invalidated during the insertion of the first half of the   */
    /*   segment.                                                             */
    finddirection(m, b, &searchtri2, endpoint2);
    if (!scoutsegment(m, b, &searchtri1, endpoint1, newmark))
    {
        /* The origin of searchtri1 may have changed if a collision with an */
        /*   intervening vertex on the segment occurred.                    */
        org(searchtri1, midvertex1);
        conformingedge(m, b, midvertex1, endpoint1, newmark);
    }
    if (!scoutsegment(m, b, &searchtri2, endpoint2, newmark))
    {
        /* The origin of searchtri2 may have changed if a collision with an */
        /*   intervening vertex on the segment occurred.                    */
        org(searchtri2, midvertex2);
        conformingedge(m, b, midvertex2, endpoint2, newmark);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  delaunayfixup()   Enforce the Delaunay condition at an edge, fanning out */
/*                    recursively from an existing vertex.  Pay special      */
/*                    attention to stacking inverted triangles.              */
/*                                                                           */
/*  This is a support routine for inserting segments into a constrained      */
/*  Delaunay triangulation.                                                  */
/*                                                                           */
/*  The origin of fixuptri is treated as if it has just been inserted, and   */
/*  the local Delaunay condition needs to be enforced.  It is only enforced  */
/*  in one sector, however, that being the angular range defined by          */
/*  fixuptri.                                                                */
/*                                                                           */
/*  This routine also needs to make decisions regarding the "stacking" of    */
/*  triangles.  (Read the description of constrainededge() below before      */
/*  reading on here, so you understand the algorithm.)  If the position of   */
/*  the new vertex (the origin of fixuptri) indicates that the vertex before */
/*  it on the polygon is a reflex vertex, then "stack" the triangle by       */
/*  doing nothing.  (fixuptri is an inverted triangle, which is how stacked  */
/*  triangles are identified.)                                               */
/*                                                                           */
/*  Otherwise, check whether the vertex before that was a reflex vertex.     */
/*  If so, perform an edge flip, thereby eliminating an inverted triangle    */
/*  (popping it off the stack).  The edge flip may result in the creation    */
/*  of a new inverted triangle, depending on whether or not the new vertex   */
/*  is visible to the vertex three edges behind on the polygon.              */
/*                                                                           */
/*  If neither of the two vertices behind the new vertex are reflex          */
/*  vertices, fixuptri and fartri, the triangle opposite it, are not         */
/*  inverted; hence, ensure that the edge between them is locally Delaunay.  */
/*                                                                           */
/*  `leftside' indicates whether or not fixuptri is to the left of the       */
/*  segment being inserted.  (Imagine that the segment is pointing up from   */
/*  endpoint1 to endpoint2.)                                                 */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::delaunayfixup(struct mesh *m,
                   struct behavior *b,
                   struct otri *fixuptri,
                   int leftside)
{
    struct otri neartri;
    struct otri fartri;
    struct osub faredge;
    vertex nearvertex, leftvertex, rightvertex, farvertex;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    lnext(*fixuptri, neartri);
    sym(neartri, fartri);
    /* Check if the edge opposite the origin of fixuptri can be flipped. */
    if (fartri.tri == m->dummytri)
    {
        return;
    }
    tspivot(neartri, faredge);
    if (faredge.ss != m->dummysub)
    {
        return;
    }
    /* Find all the relevant vertices. */
    apex(neartri, nearvertex);
    org(neartri, leftvertex);
    dest(neartri, rightvertex);
    apex(fartri, farvertex);
    /* Check whether the previous polygon vertex is a reflex vertex. */
    if (leftside)
    {
        if (counterclockwise(m, b, nearvertex, leftvertex, farvertex) <= 0.0)
        {
            /* leftvertex is a reflex vertex too.  Nothing can */
            /*   be done until a convex section is found.      */
            return;
        }
    }
    else
    {
        if (counterclockwise(m, b, farvertex, rightvertex, nearvertex) <= 0.0)
        {
            /* rightvertex is a reflex vertex too.  Nothing can */
            /*   be done until a convex section is found.       */
            return;
        }
    }
    if (counterclockwise(m, b, rightvertex, leftvertex, farvertex) > 0.0)
    {
        /* fartri is not an inverted triangle, and farvertex is not a reflex */
        /*   vertex.  As there are no reflex vertices, fixuptri isn't an     */
        /*   inverted triangle, either.  Hence, test the edge between the    */
        /*   triangles to ensure it is locally Delaunay.                     */
        if (incircle(m, b, leftvertex, farvertex, rightvertex, nearvertex) <=
            0.0)
        {
            return;
        }
        /* Not locally Delaunay; go on to an edge flip. */
    } /* else fartri is inverted; remove it from the stack by flipping. */
    flip(m, b, &neartri);
    lprevself(*fixuptri); /* Restore the origin of fixuptri after the flip. */
    /* Recursively process the two triangles that result from the flip. */
    delaunayfixup(m, b, fixuptri, leftside);
    delaunayfixup(m, b, &fartri, leftside);
}

/*****************************************************************************/
/*                                                                           */
/*  constrainededge()   Force a segment into a constrained Delaunay          */
/*                      triangulation by deleting the triangles it           */
/*                      intersects, and triangulating the polygons that      */
/*                      form on each side of it.                             */
/*                                                                           */
/*  Generates a single subsegment connecting `endpoint1' to `endpoint2'.     */
/*  The triangle `starttri' has `endpoint1' as its origin.  `newmark' is the */
/*  boundary marker of the segment.                                          */
/*                                                                           */
/*  To insert a segment, every triangle whose interior intersects the        */
/*  segment is deleted.  The union of these deleted triangles is a polygon   */
/*  (which is not necessarily monotone, but is close enough), which is       */
/*  divided into two polygons by the new segment.  This routine's task is    */
/*  to generate the Delaunay triangulation of these two polygons.            */
/*                                                                           */
/*  You might think of this routine's behavior as a two-step process.  The   */
/*  first step is to walk from endpoint1 to endpoint2, flipping each edge    */
/*  encountered.  This step creates a fan of edges connected to endpoint1,   */
/*  including the desired edge to endpoint2.  The second step enforces the   */
/*  Delaunay condition on each side of the segment in an incremental manner: */
/*  proceeding along the polygon from endpoint1 to endpoint2 (this is done   */
/*  independently on each side of the segment), each vertex is "enforced"    */
/*  as if it had just been inserted, but affecting only the previous         */
/*  vertices.  The result is the same as if the vertices had been inserted   */
/*  in the order they appear on the polygon, so the result is Delaunay.      */
/*                                                                           */
/*  In truth, constrainededge() interleaves these two steps.  The procedure  */
/*  walks from endpoint1 to endpoint2, and each time an edge is encountered  */
/*  and flipped, the newly exposed vertex (at the far end of the flipped     */
/*  edge) is "enforced" upon the previously flipped edges, usually affecting */
/*  only one side of the polygon (depending upon which side of the segment   */
/*  the vertex falls on).                                                    */
/*                                                                           */
/*  The algorithm is complicated by the need to handle polygons that are not */
/*  convex.  Although the polygon is not necessarily monotone, it can be     */
/*  triangulated in a manner similar to the stack-based algorithms for       */
/*  monotone polygons.  For each reflex vertex (local concavity) of the      */
/*  polygon, there will be an inverted triangle formed by one of the edge    */
/*  flips.  (An inverted triangle is one with negative area - that is, its   */
/*  vertices are arranged in clockwise order - and is best thought of as a   */
/*  wrinkle in the fabric of the mesh.)  Each inverted triangle can be       */
/*  thought of as a reflex vertex pushed on the stack, waiting to be fixed   */
/*  later.                                                                   */
/*                                                                           */
/*  A reflex vertex is popped from the stack when a vertex is inserted that  */
/*  is visible to the reflex vertex.  (However, if the vertex behind the     */
/*  reflex vertex is not visible to the reflex vertex, a new inverted        */
/*  triangle will take its place on the stack.)  These details are handled   */
/*  by the delaunayfixup() routine above.                                    */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::constrainededge(struct mesh *m,
                     struct behavior *b,
                     struct otri *starttri,
                     vertex endpoint2,
                     int newmark)
{
    struct otri fixuptri, fixuptri2;
    struct osub crosssubseg;
    vertex endpoint1;
    vertex farvertex;
    double area;
    int collision;
    int done;
    triangle ptr; /* Temporary variable used by sym() and oprev(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    org(*starttri, endpoint1);
    lnext(*starttri, fixuptri);
    flip(m, b, &fixuptri);
    /* `collision' indicates whether we have found a vertex directly */
    /*   between endpoint1 and endpoint2.                            */
    collision = 0;
    done      = 0;
    do
    {
        org(fixuptri, farvertex);
        /* `farvertex' is the extreme point of the polygon we are "digging" */
        /*   to get from endpoint1 to endpoint2.                           */
        if ((farvertex[0] == endpoint2[0]) && (farvertex[1] == endpoint2[1]))
        {
            oprev(fixuptri, fixuptri2);
            /* Enforce the Delaunay condition around endpoint2. */
            delaunayfixup(m, b, &fixuptri, 0);
            delaunayfixup(m, b, &fixuptri2, 1);
            done = 1;
        }
        else
        {
            /* Check whether farvertex is to the left or right of the segment */
            /*   being inserted, to decide which edge of fixuptri to dig      */
            /*   through next.                                                */
            area = counterclockwise(m, b, endpoint1, endpoint2, farvertex);
            if (area == 0.0)
            {
                /* We've collided with a vertex between endpoint1 and endpoint2.
                 */
                collision = 1;
                oprev(fixuptri, fixuptri2);
                /* Enforce the Delaunay condition around farvertex. */
                delaunayfixup(m, b, &fixuptri, 0);
                delaunayfixup(m, b, &fixuptri2, 1);
                done = 1;
            }
            else
            {
                if (area > 0.0)
                { /* farvertex is to the left of the segment. */
                    oprev(fixuptri, fixuptri2);
                    /* Enforce the Delaunay condition around farvertex, on the
                     */
                    /*   left side of the segment only. */
                    delaunayfixup(m, b, &fixuptri2, 1);
                    /* Flip the edge that crosses the segment.  After the edge
                     * is */
                    /*   flipped, one of its endpoints is the fan vertex, and
                     * the */
                    /*   destination of fixuptri is the fan vertex. */
                    lprevself(fixuptri);
                }
                else
                { /* farvertex is to the right of the segment. */
                    delaunayfixup(m, b, &fixuptri, 0);
                    /* Flip the edge that crosses the segment.  After the edge
                     * is */
                    /*   flipped, one of its endpoints is the fan vertex, and
                     * the */
                    /*   destination of fixuptri is the fan vertex. */
                    oprevself(fixuptri);
                }
                /* Check for two intersecting segments. */
                tspivot(fixuptri, crosssubseg);
                if (crosssubseg.ss == m->dummysub)
                {
                    flip(m,
                         b,
                         &fixuptri); /* May create inverted triangle at left. */
                }
                else
                {
                    /* We've collided with a segment between endpoint1 and
                     * endpoint2. */
                    collision = 1;
                    /* Insert a vertex at the intersection. */
                    segmentintersection(
                        m, b, &fixuptri, &crosssubseg, endpoint2);
                    done = 1;
                }
            }
        }
    } while (!done);
    /* Insert a subsegment to make the segment permanent. */
    insertsubseg(m, b, &fixuptri, newmark);
    /* If there was a collision with an interceding vertex, install another */
    /*   segment connecting that vertex with endpoint2.                     */
    if (collision)
    {
        /* Insert the remainder of the segment. */
        if (!scoutsegment(m, b, &fixuptri, endpoint2, newmark))
        {
            constrainededge(m, b, &fixuptri, endpoint2, newmark);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  insertsegment()   Insert a PSLG segment into a triangulation.            */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::insertsegment(struct mesh *m,
                   struct behavior *b,
                   vertex endpoint1,
                   vertex endpoint2,
                   int newmark)
{
    struct otri searchtri1, searchtri2;
    triangle encodedtri;
    vertex checkvertex;
    triangle ptr; /* Temporary variable used by sym(). */



    /* Find a triangle whose origin is the segment's first endpoint. */
    checkvertex = (vertex)NULL;
    encodedtri  = vertex2tri(endpoint1);
    if (encodedtri != (triangle)NULL)
    {
        decode(encodedtri, searchtri1);
        org(searchtri1, checkvertex);
    }
    if (checkvertex != endpoint1)
    {
        /* Find a boundary triangle to search from. */
        searchtri1.tri    = m->dummytri;
        searchtri1.orient = 0;
        symself(searchtri1);
        /* Search for the segment's first endpoint by point location. */
        if (locate(m, b, endpoint1, &searchtri1) != ONVERTEX)
        {
            printf("Internal error in insertsegment():  Unable to locate PSLG "
                   "vertex\n");
            printf("  (%.12g, %.12g) in triangulation.\n",
                   endpoint1[0],
                   endpoint1[1]);
            internalerror();
        }
    }
    /* Remember this triangle to improve subsequent point location. */
    otricopy(searchtri1, m->recenttri);
    /* Scout the beginnings of a path from the first endpoint */
    /*   toward the second.                                   */
    if (scoutsegment(m, b, &searchtri1, endpoint2, newmark))
    {
        /* The segment was easily inserted. */
        return;
    }
    /* The first endpoint may have changed if a collision with an intervening */
    /*   vertex on the segment occurred.                                      */
    org(searchtri1, endpoint1);

    /* Find a triangle whose origin is the segment's second endpoint. */
    checkvertex = (vertex)NULL;
    encodedtri  = vertex2tri(endpoint2);
    if (encodedtri != (triangle)NULL)
    {
        decode(encodedtri, searchtri2);
        org(searchtri2, checkvertex);
    }
    if (checkvertex != endpoint2)
    {
        /* Find a boundary triangle to search from. */
        searchtri2.tri    = m->dummytri;
        searchtri2.orient = 0;
        symself(searchtri2);
        /* Search for the segment's second endpoint by point location. */
        if (locate(m, b, endpoint2, &searchtri2) != ONVERTEX)
        {
            printf("Internal error in insertsegment():  Unable to locate PSLG "
                   "vertex\n");
            printf("  (%.12g, %.12g) in triangulation.\n",
                   endpoint2[0],
                   endpoint2[1]);
            internalerror();
        }
    }
    /* Remember this triangle to improve subsequent point location. */
    otricopy(searchtri2, m->recenttri);
    /* Scout the beginnings of a path from the second endpoint */
    /*   toward the first.                                     */
    if (scoutsegment(m, b, &searchtri2, endpoint1, newmark))
    {
        /* The segment was easily inserted. */
        return;
    }
    /* The second endpoint may have changed if a collision with an intervening
     */
    /*   vertex on the segment occurred. */
    org(searchtri2, endpoint2);


    /* Insert the segment directly into the triangulation. */
    constrainededge(m, b, &searchtri1, endpoint2, newmark);
}

/*****************************************************************************/
/*                                                                           */
/*  markhull()   Cover the convex hull of a triangulation with subsegments.  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::markhull(struct mesh *m, struct behavior *b)
{
    struct otri hulltri;
    struct otri nexttri;
    struct otri starttri;
    triangle ptr; /* Temporary variable used by sym() and oprev(). */

    /* Find a triangle handle on the hull. */
    hulltri.tri    = m->dummytri;
    hulltri.orient = 0;
    symself(hulltri);
    /* Remember where we started so we know when to stop. */
    otricopy(hulltri, starttri);
    /* Go once counterclockwise around the convex hull. */
    do
    {
        /* Create a subsegment if there isn't already one here. */
        insertsubseg(m, b, &hulltri, 1);
        /* To find the next hull edge, go clockwise around the next vertex. */
        lnextself(hulltri);
        oprev(hulltri, nexttri);
        while (nexttri.tri != m->dummytri)
        {
            otricopy(nexttri, hulltri);
            oprev(hulltri, nexttri);
        }
    } while (!otriequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  formskeleton()   Create the segments of a triangulation, including PSLG  */
/*                   segments and edges on the convex hull.                  */
/*                                                                           */
/*  The PSLG segments are read from a .poly file.  The return value is the   */
/*  number of segments in the file.                                          */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::formskeleton(struct mesh *m,
                  struct behavior *b,
                  int *segmentlist,
                  int *segmentmarkerlist,
                  int numberofsegments)
{
    char polyfilename[6];
    int index;
    vertex endpoint1, endpoint2;
    int segmentmarkers;
    int end1, end2;
    int boundmarker;
    int i;

    if (b->poly)
    {
        strcpy(polyfilename, "input");
        m->insegments  = numberofsegments;
        segmentmarkers = segmentmarkerlist != (int *)NULL;
        index          = 0;
        /* If the input vertices are collinear, there is no triangulation, */
        /*   so don't try to insert segments.                              */
        if (m->triangles.items == 0)
        {
            return;
        }

        /* If segments are to be inserted, compute a mapping */
        /*   from vertices to triangles.                     */
        if (m->insegments > 0)
        {
            makevertexmap(m, b);

        }

        boundmarker = 0;
        /* Read and insert the segments. */
        for (i = 0; i < m->insegments; i++)
        {
            end1 = segmentlist[index++];
            end2 = segmentlist[index++];
            if (segmentmarkers)
            {
                boundmarker = segmentmarkerlist[i];
            }
            if ((end1 < 0) ||
                (end1 >= 0 + m->invertices))
            {

            }
            else if ((end2 < 0) ||
                     (end2 >= 0 + m->invertices))
            {

            }
            else
            {
                /* Find the vertices numbered `end1' and `end2'. */
                endpoint1 = getvertex(m, b, end1);
                endpoint2 = getvertex(m, b, end2);
                if ((endpoint1[0] == endpoint2[0]) &&
                    (endpoint1[1] == endpoint2[1]))
                {

                }
                else
                {
                    insertsegment(m, b, endpoint1, endpoint2, boundmarker);
                }
            }
        }
    }
    else
    {
        m->insegments = 0;
    }
    if (!b->poly)
    {
        /* Enclose the convex hull with subsegments. */

        markhull(m, b);
    }
}

/**                                                                         **/
/**                                                                         **/
/********* Segment insertion ends here                               *********/

/********* Carving out holes and concavities begins here             *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  infecthull()   Virally infect all of the triangles of the convex hull    */
/*                 that are not protected by subsegments.  Where there are   */
/*                 subsegments, set boundary markers as appropriate.         */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::infecthull(struct mesh *m, struct behavior *b)
{
    struct otri hulltri;
    struct otri nexttri;
    struct otri starttri;
    struct osub hullsubseg;
    triangle **deadtriangle;
    vertex horg, hdest;
    triangle ptr; /* Temporary variable used by sym(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Find a triangle handle on the hull. */
    hulltri.tri    = m->dummytri;
    hulltri.orient = 0;
    symself(hulltri);
    /* Remember where we started so we know when to stop. */
    otricopy(hulltri, starttri);
    /* Go once counterclockwise around the convex hull. */
    do
    {
        /* Ignore triangles that are already infected. */
        if (!infected(hulltri))
        {
            /* Is the triangle protected by a subsegment? */
            tspivot(hulltri, hullsubseg);
            if (hullsubseg.ss == m->dummysub)
            {
                /* The triangle is not protected; infect it. */
                if (!infected(hulltri))
                {
                    infect(hulltri);
                    deadtriangle  = (triangle **)poolalloc(&m->viri);
                    *deadtriangle = hulltri.tri;
                }
            }
            else
            {
                /* The triangle is protected; set boundary markers if
                 * appropriate. */
                if (mark(hullsubseg) == 0)
                {
                    setmark(hullsubseg, 1);
                    org(hulltri, horg);
                    dest(hulltri, hdest);
                    if (vertexmark(horg) == 0)
                    {
                        setvertexmark(horg, 1);
                    }
                    if (vertexmark(hdest) == 0)
                    {
                        setvertexmark(hdest, 1);
                    }
                }
            }
        }
        /* To find the next hull edge, go clockwise around the next vertex. */
        lnextself(hulltri);
        oprev(hulltri, nexttri);
        while (nexttri.tri != m->dummytri)
        {
            otricopy(nexttri, hulltri);
            oprev(hulltri, nexttri);
        }
    } while (!otriequal(hulltri, starttri));
}

/*****************************************************************************/
/*                                                                           */
/*  plague()   Spread the virus from all infected triangles to any neighbors */
/*             not protected by subsegments.  Delete all infected triangles. */
/*                                                                           */
/*  This is the procedure that actually creates holes and concavities.       */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase identifies all   */
/*  the triangles that will die, and marks them as infected.  They are       */
/*  marked to ensure that each triangle is added to the virus pool only      */
/*  once, so the procedure will terminate.                                   */
/*                                                                           */
/*  The second phase actually eliminates the infected triangles.  It also    */
/*  eliminates orphaned vertices.                                            */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::plague(struct mesh *m, struct behavior *b)
{
    struct otri testtri;
    struct otri neighbor;
    triangle **virusloop;
    triangle **deadtriangle;
    struct osub neighborsubseg;
    vertex testvertex;
    vertex norg, ndest;
    int killorg;
    triangle ptr; /* Temporary variable used by sym() and onext(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */

    /* Loop through all the infected triangles, spreading the virus to */
    /*   their neighbors, then to their neighbors' neighbors.          */
    traversalinit(&m->viri);
    virusloop = (triangle **)traverse(&m->viri);
    while (virusloop != (triangle **)NULL)
    {
        testtri.tri = *virusloop;
        /* A triangle is marked as infected by messing with one of its pointers
         */
        /*   to subsegments, setting it to an illegal value.  Hence, we have to
         */
        /*   temporarily uninfect this triangle so that we can examine its */
        /*   adjacent subsegments. */
        uninfect(testtri);

        /* Check each of the triangle's three neighbors. */
        for (testtri.orient = 0; testtri.orient < 3; testtri.orient++)
        {
            /* Find the neighbor. */
            sym(testtri, neighbor);
            /* Check for a subsegment between the triangle and its neighbor. */
            tspivot(testtri, neighborsubseg);
            /* Check if the neighbor is nonexistent or already infected. */
            if ((neighbor.tri == m->dummytri) || infected(neighbor))
            {
                if (neighborsubseg.ss != m->dummysub)
                {
                    /* There is a subsegment separating the triangle from its */
                    /*   neighbor, but both triangles are dying, so the
                     * subsegment */
                    /*   dies too. */
                    subsegdealloc(m, neighborsubseg.ss);
                    if (neighbor.tri != m->dummytri)
                    {
                        /* Make sure the subsegment doesn't get deallocated
                         * again */
                        /*   later when the infected neighbor is visited. */
                        uninfect(neighbor);
                        tsdissolve(neighbor);
                        infect(neighbor);
                    }
                }
            }
            else
            { /* The neighbor exists and is not infected. */
                if (neighborsubseg.ss == m->dummysub)
                {
                    /* There is no subsegment protecting the neighbor, so */
                    /*   the neighbor becomes infected.                   */

                    infect(neighbor);
                    /* Ensure that the neighbor's neighbors will be infected. */
                    deadtriangle  = (triangle **)poolalloc(&m->viri);
                    *deadtriangle = neighbor.tri;
                }
                else
                { /* The neighbor is protected by a subsegment. */
                    /* Remove this triangle from the subsegment. */
                    stdissolve(neighborsubseg);
                    /* The subsegment becomes a boundary.  Set markers
                     * accordingly. */
                    if (mark(neighborsubseg) == 0)
                    {
                        setmark(neighborsubseg, 1);
                    }
                    org(neighbor, norg);
                    dest(neighbor, ndest);
                    if (vertexmark(norg) == 0)
                    {
                        setvertexmark(norg, 1);
                    }
                    if (vertexmark(ndest) == 0)
                    {
                        setvertexmark(ndest, 1);
                    }
                }
            }
        }
        /* Remark the triangle as infected, so it doesn't get added to the */
        /*   virus pool again.                                             */
        infect(testtri);
        virusloop = (triangle **)traverse(&m->viri);
    }


    traversalinit(&m->viri);
    virusloop = (triangle **)traverse(&m->viri);
    while (virusloop != (triangle **)NULL)
    {
        testtri.tri = *virusloop;

        /* Check each of the three corners of the triangle for elimination. */
        /*   This is done by walking around each vertex, checking if it is  */
        /*   still connected to at least one live triangle.                 */
        for (testtri.orient = 0; testtri.orient < 3; testtri.orient++)
        {
            org(testtri, testvertex);
            /* Check if the vertex has already been tested. */
            if (testvertex != (vertex)NULL)
            {
                killorg = 1;
                /* Mark the corner of the triangle as having been tested. */
                setorg(testtri, NULL);
                /* Walk counterclockwise about the vertex. */
                onext(testtri, neighbor);
                /* Stop upon reaching a boundary or the starting triangle. */
                while ((neighbor.tri != m->dummytri) &&
                       (!otriequal(neighbor, testtri)))
                {
                    if (infected(neighbor))
                    {
                        /* Mark the corner of this triangle as having been
                         * tested. */
                        setorg(neighbor, NULL);
                    }
                    else
                    {
                        /* A live triangle.  The vertex survives. */
                        killorg = 0;
                    }
                    /* Walk counterclockwise about the vertex. */
                    onextself(neighbor);
                }
                /* If we reached a boundary, we must walk clockwise as well. */
                if (neighbor.tri == m->dummytri)
                {
                    /* Walk clockwise about the vertex. */
                    oprev(testtri, neighbor);
                    /* Stop upon reaching a boundary. */
                    while (neighbor.tri != m->dummytri)
                    {
                        if (infected(neighbor))
                        {
                            /* Mark the corner of this triangle as having been
                             * tested. */
                            setorg(neighbor, NULL);
                        }
                        else
                        {
                            /* A live triangle.  The vertex survives. */
                            killorg = 0;
                        }
                        /* Walk clockwise about the vertex. */
                        oprevself(neighbor);
                    }
                }
                if (killorg)
                {

                    setvertextype(testvertex, UNDEADVERTEX);
                    m->undeads++;
                }
            }
        }

        /* Record changes in the number of boundary edges, and disconnect */
        /*   dead triangles from their neighbors.                         */
        for (testtri.orient = 0; testtri.orient < 3; testtri.orient++)
        {
            sym(testtri, neighbor);
            if (neighbor.tri == m->dummytri)
            {
                /* There is no neighboring triangle on this edge, so this edge
                 */
                /*   is a boundary edge.  This triangle is being deleted, so
                 * this */
                /*   boundary edge is deleted. */
                m->hullsize--;
            }
            else
            {
                /* Disconnect the triangle from its neighbor. */
                dissolve(neighbor);
                /* There is a neighboring triangle on this edge, so this edge */
                /*   becomes a boundary edge when this triangle is deleted.   */
                m->hullsize++;
            }
        }
        /* Return the dead triangle to the pool of triangles. */
        triangledealloc(m, testtri.tri);
        virusloop = (triangle **)traverse(&m->viri);
    }
    /* Empty the virus pool. */
    poolrestart(&m->viri);
}

/*****************************************************************************/
/*                                                                           */
/*  regionplague()   Spread regional attributes and/or area constraints      */
/*                   (from a .poly file) throughout the mesh.                */
/*                                                                           */
/*  This procedure operates in two phases.  The first phase spreads an       */
/*  attribute and/or an area constraint through a (segment-bounded) region.  */
/*  The triangles are marked to ensure that each triangle is added to the    */
/*  virus pool only once, so the procedure will terminate.                   */
/*                                                                           */
/*  The second phase uninfects all infected triangles, returning them to     */
/*  normal.                                                                  */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::regionplague(struct mesh *m,
                  struct behavior *b,
                  double attribute,
                  double area)
{
    struct otri testtri;
    struct otri neighbor;
    triangle **virusloop;
    triangle **regiontri;
    struct osub neighborsubseg;
    triangle ptr; /* Temporary variable used by sym() and onext(). */
    subseg sptr;  /* Temporary variable used by tspivot(). */


    /* Loop through all the infected triangles, spreading the attribute      */
    /*   and/or area constraint to their neighbors, then to their neighbors' */
    /*   neighbors.                                                          */
    traversalinit(&m->viri);
    virusloop = (triangle **)traverse(&m->viri);
    while (virusloop != (triangle **)NULL)
    {
        testtri.tri = *virusloop;
        /* A triangle is marked as infected by messing with one of its pointers
         */
        /*   to subsegments, setting it to an illegal value.  Hence, we have to
         */
        /*   temporarily uninfect this triangle so that we can examine its */
        /*   adjacent subsegments. */
        uninfect(testtri);

        /* Check each of the triangle's three neighbors. */
        for (testtri.orient = 0; testtri.orient < 3; testtri.orient++)
        {
            /* Find the neighbor. */
            sym(testtri, neighbor);
            /* Check for a subsegment between the triangle and its neighbor. */
            tspivot(testtri, neighborsubseg);
            /* Make sure the neighbor exists, is not already infected, and */
            /*   isn't protected by a subsegment.                          */
            if ((neighbor.tri != m->dummytri) && !infected(neighbor) &&
                (neighborsubseg.ss == m->dummysub))
            {

                /* Infect the neighbor. */
                infect(neighbor);
                /* Ensure that the neighbor's neighbors will be infected. */
                regiontri  = (triangle **)poolalloc(&m->viri);
                *regiontri = neighbor.tri;
            }
        }
        /* Remark the triangle as infected, so it doesn't get added to the */
        /*   virus pool again.                                             */
        infect(testtri);
        virusloop = (triangle **)traverse(&m->viri);
    }

    /* Uninfect all triangles. */

    traversalinit(&m->viri);
    virusloop = (triangle **)traverse(&m->viri);
    while (virusloop != (triangle **)NULL)
    {
        testtri.tri = *virusloop;
        uninfect(testtri);
        virusloop = (triangle **)traverse(&m->viri);
    }
    /* Empty the virus pool. */
    poolrestart(&m->viri);
}

/*****************************************************************************/
/*                                                                           */
/*  carveholes()   Find the holes and infect them.  Find the area            */
/*                 constraints and infect them.  Infect the convex hull.     */
/*                 Spread the infection and kill triangles.  Spread the      */
/*                 area constraints.                                         */
/*                                                                           */
/*  This routine mainly calls other routines to carry out all these          */
/*  functions.                                                               */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::carveholes(struct mesh *m,
                struct behavior *b,
                double *holelist,
                int holes,
                double *regionlist,
                int regions)
{
    struct otri searchtri;
    struct otri *regiontris;
    triangle **holetri;
    triangle **regiontri;
    vertex searchorg, searchdest;
    enum locateresult intersect;
    int i;
    triangle ptr; /* Temporary variable used by sym(). */

    if (regions > 0)
    {
        /* Allocate storage for the triangles in which region points fall. */
        regiontris =
            (struct otri *)trimalloc(regions * (int)sizeof(struct otri));
    }
    else
    {
        regiontris = (struct otri *)NULL;
    }


    /* Initialize a pool of viri to be used for holes, concavities, */
    /*   regional attributes, and/or regional area constraints.     */
    poolinit(&m->viri, sizeof(triangle *), VIRUSPERBLOCK, VIRUSPERBLOCK, 0);

    /* Mark as infected any unprotected triangles on the boundary. */
    /*   This is one way by which concavities are created.         */
    infecthull(m, b);

    if (holes > 0)
    {
        /* Infect each triangle in which a hole lies. */
        for (i = 0; i < 2 * holes; i += 2)
        {
            /* Ignore holes that aren't within the bounds of the mesh. */
            if ((holelist[i] >= m->xmin) && (holelist[i] <= m->xmax) &&
                (holelist[i + 1] >= m->ymin) && (holelist[i + 1] <= m->ymax))
            {
                /* Start searching from some triangle on the outer boundary. */
                searchtri.tri    = m->dummytri;
                searchtri.orient = 0;
                symself(searchtri);
                /* Ensure that the hole is to the left of this boundary edge; */
                /*   otherwise, locate() will falsely report that the hole    */
                /*   falls within the starting triangle.                      */
                org(searchtri, searchorg);
                dest(searchtri, searchdest);
                if (counterclockwise(
                        m, b, searchorg, searchdest, &holelist[i]) > 0.0)
                {
                    /* Find a triangle that contains the hole. */
                    intersect = locate(m, b, &holelist[i], &searchtri);
                    if ((intersect != OUTSIDE) && (!infected(searchtri)))
                    {
                        /* Infect the triangle.  This is done by marking the
                         * triangle  */
                        /*   as infected and including the triangle in the virus
                         * pool. */
                        infect(searchtri);
                        holetri  = (triangle **)poolalloc(&m->viri);
                        *holetri = searchtri.tri;
                    }
                }
            }
        }
    }

    /* Now, we have to find all the regions BEFORE we carve the holes, because
     */
    /*   locate() won't work when the triangulation is no longer convex. */
    /*   (Incidentally, this is the reason why regional attributes and area */
    /*   constraints can't be used when refining a preexisting mesh, which */
    /*   might not be convex; they can only be used with a freshly */
    /*   triangulated PSLG.) */
    if (regions > 0)
    {
        /* Find the starting triangle for each region. */
        for (i = 0; i < regions; i++)
        {
            regiontris[i].tri = m->dummytri;
            /* Ignore region points that aren't within the bounds of the mesh.
             */
            if ((regionlist[4 * i] >= m->xmin) &&
                (regionlist[4 * i] <= m->xmax) &&
                (regionlist[4 * i + 1] >= m->ymin) &&
                (regionlist[4 * i + 1] <= m->ymax))
            {
                /* Start searching from some triangle on the outer boundary. */
                searchtri.tri    = m->dummytri;
                searchtri.orient = 0;
                symself(searchtri);
                /* Ensure that the region point is to the left of this boundary
                 */
                /*   edge; otherwise, locate() will falsely report that the */
                /*   region point falls within the starting triangle. */
                org(searchtri, searchorg);
                dest(searchtri, searchdest);
                if (counterclockwise(
                        m, b, searchorg, searchdest, &regionlist[4 * i]) > 0.0)
                {
                    /* Find a triangle that contains the region point. */
                    intersect = locate(m, b, &regionlist[4 * i], &searchtri);
                    if ((intersect != OUTSIDE) && (!infected(searchtri)))
                    {
                        /* Record the triangle for processing after the */
                        /*   holes have been carved.                    */
                        otricopy(searchtri, regiontris[i]);
                    }
                }
            }
        }
    }

    if (m->viri.items > 0)
    {
        /* Carve the holes and concavities. */
        plague(m, b);
    }
    /* The virus pool should be empty now. */

    if (regions > 0)
    {
        for (i = 0; i < regions; i++)
        {
            if (regiontris[i].tri != m->dummytri)
            {
                /* Make sure the triangle under consideration still exists. */
                /*   It may have been eaten by the virus.                   */
                if (!deadtri(regiontris[i].tri))
                {
                    /* Put one triangle in the virus pool. */
                    infect(regiontris[i]);
                    regiontri  = (triangle **)poolalloc(&m->viri);
                    *regiontri = regiontris[i].tri;
                    /* Apply one region's attribute and/or area constraint. */
                    regionplague(
                        m, b, regionlist[4 * i + 2], regionlist[4 * i + 3]);
                    /* The virus pool should be empty now. */
                }
            }
        }
    }


    pooldeinit(&m->viri);

    if (regions > 0)
    {
        trifree((void *)regiontris);
    }
}

/**                                                                         **/
/**                                                                         **/
/********* Carving out holes and concavities ends here               *********/

/********* Mesh quality maintenance begins here                      *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  tallyencs()   Traverse the entire list of subsegments, and check each    */
/*                to see if it is encroached.  If so, add it to the list.    */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::tallyencs(struct mesh *m, struct behavior *b)
{
    struct osub subsegloop;

    traversalinit(&m->subsegs);
    subsegloop.ssorient = 0;
    subsegloop.ss       = subsegtraverse(m);
    while (subsegloop.ss != (subseg *)NULL)
    {
        /* If the segment is encroached, add it to the list. */
        checkseg4encroach(m, b, &subsegloop);
        subsegloop.ss = subsegtraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  precisionerror()  Print an error message for precision problems.         */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::precisionerror()
{
    printf("Try increasing the area criterion and/or reducing the minimum\n");
    printf("  allowable angle so that tiny triangles are not created.\n");
#ifdef SINGLE
    printf("Alternatively, try recompiling me with double precision\n");
    printf("  arithmetic (by removing \"#define SINGLE\" from the\n");
    printf("  source file or \"-DSINGLE\" from the makefile).\n");
#endif /* SINGLE */
}

/*****************************************************************************/
/*                                                                           */
/*  splitencsegs()   Split all the encroached subsegments.                   */
/*                                                                           */
/*  Each encroached subsegment is repaired by splitting it - inserting a     */
/*  vertex at or near its midpoint.  Newly inserted vertices may encroach    */
/*  upon other subsegments; these are also repaired.                         */
/*                                                                           */
/*  `triflaws' is a flag that specifies whether one should take note of new  */
/*  bad triangles that result from inserting vertices to repair encroached   */
/*  subsegments.                                                             */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::splitencsegs(struct mesh *m, struct behavior *b, int triflaws)
{
    struct otri enctri;
    struct otri testtri;
    struct osub testsh;
    struct osub currentenc;
    struct badsubseg *encloop;
    vertex eorg, edest, eapex;
    vertex newvertex;
    enum insertvertexresult success;
    double segmentlength, nearestpoweroftwo;
    double split;
    double multiplier, divisor;
    int acuteorg, acuteorg2, acutedest, acutedest2;
    int i;
    triangle ptr; /* Temporary variable used by stpivot(). */
    subseg sptr;  /* Temporary variable used by snext(). */

    /* Note that steinerleft == -1 if an unlimited number */
    /*   of Steiner points is allowed.                    */
    while ((m->badsubsegs.items > 0) && (m->steinerleft != 0))
    {
        traversalinit(&m->badsubsegs);
        encloop = badsubsegtraverse(m);
        while ((encloop != (struct badsubseg *)NULL) && (m->steinerleft != 0))
        {
            sdecode(encloop->encsubseg, currentenc);
            sorg(currentenc, eorg);
            sdest(currentenc, edest);
            /* Make sure that this segment is still the same segment it was   */
            /*   when it was determined to be encroached.  If the segment was */
            /*   enqueued multiple times (because several newly inserted      */
            /*   vertices encroached it), it may have already been split.     */
            if (!deadsubseg(currentenc.ss) && (eorg == encloop->subsegorg) &&
                (edest == encloop->subsegdest))
            {
                /* To decide where to split a segment, we need to know if the */
                /*   segment shares an endpoint with an adjacent segment. */
                /*   The concern is that, if we simply split every encroached */
                /*   segment in its center, two adjacent segments with a small
                 */
                /*   angle between them might lead to an infinite loop; each */
                /*   vertex added to split one segment will encroach upon the */
                /*   other segment, which must then be split with a vertex that
                 */
                /*   will encroach upon the first segment, and so on forever. */
                /* To avoid this, imagine a set of concentric circles, whose */
                /*   radii are powers of two, about each segment endpoint. */
                /*   These concentric circles determine where the segment is */
                /*   split.  (If both endpoints are shared with adjacent */
                /*   segments, split the segment in the middle, and apply the */
                /*   concentric circles for later splittings.) */

                /* Is the origin shared with another segment? */
                stpivot(currentenc, enctri);
                lnext(enctri, testtri);
                tspivot(testtri, testsh);
                acuteorg = testsh.ss != m->dummysub;
                /* Is the destination shared with another segment? */
                lnextself(testtri);
                tspivot(testtri, testsh);
                acutedest = testsh.ss != m->dummysub;

                /* If we're using Chew's algorithm (rather than Ruppert's) */
                /*   to define encroachment, delete free vertices from the */
                /*   subsegment's diametral circle.                        */
                if (!acuteorg && !acutedest)
                {
                    apex(enctri, eapex);
                    while ((vertextype(eapex) == FREEVERTEX) &&
                           ((eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                                (eorg[1] - eapex[1]) * (edest[1] - eapex[1]) <
                            0.0))
                    {
                        deletevertex(m, b, &testtri);
                        stpivot(currentenc, enctri);
                        apex(enctri, eapex);
                        lprev(enctri, testtri);
                    }
                }

                /* Now, check the other side of the segment, if there's a
                 * triangle */
                /*   there. */
                sym(enctri, testtri);
                if (testtri.tri != m->dummytri)
                {
                    /* Is the destination shared with another segment? */
                    lnextself(testtri);
                    tspivot(testtri, testsh);
                    acutedest2 = testsh.ss != m->dummysub;
                    acutedest  = acutedest || acutedest2;
                    /* Is the origin shared with another segment? */
                    lnextself(testtri);
                    tspivot(testtri, testsh);
                    acuteorg2 = testsh.ss != m->dummysub;
                    acuteorg  = acuteorg || acuteorg2;

                    /* Delete free vertices from the subsegment's diametral
                     * circle. */
                    if (!acuteorg2 && !acutedest2)
                    {
                        org(testtri, eapex);
                        while (
                            (vertextype(eapex) == FREEVERTEX) &&
                            ((eorg[0] - eapex[0]) * (edest[0] - eapex[0]) +
                                 (eorg[1] - eapex[1]) * (edest[1] - eapex[1]) <
                             0.0))
                        {
                            deletevertex(m, b, &testtri);
                            sym(enctri, testtri);
                            apex(testtri, eapex);
                            lprevself(testtri);
                        }
                    }
                }

                /* Use the concentric circles if exactly one endpoint is shared
                 */
                /*   with another adjacent segment. */
                if (acuteorg || acutedest)
                {
                    segmentlength =
                        sqrt((edest[0] - eorg[0]) * (edest[0] - eorg[0]) +
                             (edest[1] - eorg[1]) * (edest[1] - eorg[1]));
                    /* Find the power of two that most evenly splits the
                     * segment.  */
                    /*   The worst case is a 2:1 ratio between subsegment
                     * lengths. */
                    nearestpoweroftwo = 1.0;
                    while (segmentlength > 3.0 * nearestpoweroftwo)
                    {
                        nearestpoweroftwo *= 2.0;
                    }
                    while (segmentlength < 1.5 * nearestpoweroftwo)
                    {
                        nearestpoweroftwo *= 0.5;
                    }
                    /* Where do we split the segment? */
                    split = nearestpoweroftwo / segmentlength;
                    if (acutedest)
                    {
                        split = 1.0 - split;
                    }
                }
                else
                {
                    /* If we're not worried about adjacent segments, split */
                    /*   this segment in the middle.                       */
                    split = 0.5;
                }

                /* Create the new vertex. */
                newvertex = (vertex)poolalloc(&m->vertices);
                /* Interpolate its coordinate and attributes. */
                for (i = 0; i < 2 + m->nextras; i++)
                {
                    newvertex[i] = eorg[i] + split * (edest[i] - eorg[i]);
                }

                /* Roundoff in the above calculation may yield a `newvertex'
                 */
                /*   that is not precisely collinear with `eorg' and
                 * `edest'.  */
                /*   Improve collinearity by one step of iterative
                 * refinement. */
                multiplier = counterclockwise(m, b, eorg, edest, newvertex);
                divisor    = ((eorg[0] - edest[0]) * (eorg[0] - edest[0]) +
                           (eorg[1] - edest[1]) * (eorg[1] - edest[1]));
                if ((multiplier != 0.0) && (divisor != 0.0))
                {
                    multiplier = multiplier / divisor;
                    /* Watch out for NANs. */
                    if (multiplier == multiplier)
                    {
                        newvertex[0] += multiplier * (edest[1] - eorg[1]);
                        newvertex[1] += multiplier * (eorg[0] - edest[0]);
                    }
                }

                setvertexmark(newvertex, mark(currentenc));
                setvertextype(newvertex, SEGMENTVERTEX);

                /* Check whether the new vertex lies on an endpoint. */
                if (((newvertex[0] == eorg[0]) && (newvertex[1] == eorg[1])) ||
                    ((newvertex[0] == edest[0]) && (newvertex[1] == edest[1])))
                {
                    printf("Error:  Ran out of precision at (%.12g, %.12g).\n",
                           newvertex[0],
                           newvertex[1]);
                    printf("I attempted to split a segment to a smaller size "
                           "than\n");
                    printf(
                        "  can be accommodated by the finite precision of\n");
                    printf("  floating point arithmetic.\n");
                    precisionerror();
                    triexit(1);
                }
                /* Insert the splitting vertex.  This should always succeed. */
                success = insertvertex(
                    m, b, newvertex, &enctri, &currentenc, 1, triflaws);
                if ((success != SUCCESSFULVERTEX) &&
                    (success != ENCROACHINGVERTEX))
                {
                    printf("Internal error in splitencsegs():\n");
                    printf("  Failure to split a segment.\n");
                    internalerror();
                }
                if (m->steinerleft > 0)
                {
                    m->steinerleft--;
                }
                /* Check the two new subsegments to see if they're encroached.
                 */
                checkseg4encroach(m, b, &currentenc);
                snextself(currentenc);
                checkseg4encroach(m, b, &currentenc);
            }

            badsubsegdealloc(m, encloop);
            encloop = badsubsegtraverse(m);
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  tallyfaces()   Test every triangle in the mesh for quality measures.     */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::tallyfaces(struct mesh *m, struct behavior *b)
{
    struct otri triangleloop;


    traversalinit(&m->triangles);
    triangleloop.orient = 0;
    triangleloop.tri    = triangletraverse(m);
    while (triangleloop.tri != (triangle *)NULL)
    {
        /* If the triangle is bad, enqueue it. */
        testtriangle(m, b, &triangleloop);
        triangleloop.tri = triangletraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  splittriangle()   Inserts a vertex at the circumcenter of a triangle.    */
/*                    Deletes the newly inserted vertex if it encroaches     */
/*                    upon a segment.                                        */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::splittriangle(struct mesh *m, struct behavior *b, struct badtriang *badtri)
{
    struct otri badotri;
    vertex borg, bdest, bapex;
    vertex newvertex;
    double xi, eta;
    enum insertvertexresult success;
    int errorflag;
    int i;

    decode(badtri->poortri, badotri);
    org(badotri, borg);
    dest(badotri, bdest);
    apex(badotri, bapex);
    /* Make sure that this triangle is still the same triangle it was      */
    /*   when it was tested and determined to be of bad quality.           */
    /*   Subsequent transformations may have made it a different triangle. */
    if (!deadtri(badotri.tri) && (borg == badtri->triangorg) &&
        (bdest == badtri->triangdest) && (bapex == badtri->triangapex))
    {


        errorflag = 0;
        /* Create a new vertex at the triangle's circumcenter. */
        newvertex = (vertex)poolalloc(&m->vertices);
        findcircumcenter(m, b, borg, bdest, bapex, newvertex, &xi, &eta, 1);

        /* Check whether the new vertex lies on a triangle vertex. */
        if (((newvertex[0] == borg[0]) && (newvertex[1] == borg[1])) ||
            ((newvertex[0] == bdest[0]) && (newvertex[1] == bdest[1])) ||
            ((newvertex[0] == bapex[0]) && (newvertex[1] == bapex[1])))
        {
            vertexdealloc(m, newvertex);
        }
        else
        {
            for (i = 2; i < 2 + m->nextras; i++)
            {
                /* Interpolate the vertex attributes at the circumcenter. */
                newvertex[i] = borg[i] + xi * (bdest[i] - borg[i]) +
                               eta * (bapex[i] - borg[i]);
            }
            /* The new vertex must be in the interior, and therefore is a */
            /*   free vertex with a marker of zero.                       */
            setvertexmark(newvertex, 0);
            setvertextype(newvertex, FREEVERTEX);

            /* Ensure that the handle `badotri' does not represent the longest
             */
            /*   edge of the triangle.  This ensures that the circumcenter must
             */
            /*   fall to the left of this edge, so point location will work. */
            /*   (If the angle org-apex-dest exceeds 90 degrees, then the */
            /*   circumcenter lies outside the org-dest edge, and eta is */
            /*   negative.  Roundoff error might prevent eta from being */
            /*   negative when it should be, so I test eta against xi.) */
            if (eta < xi)
            {
                lprevself(badotri);
            }

            /* Insert the circumcenter, searching from the edge of the triangle,
             */
            /*   and maintain the Delaunay property of the triangulation. */
            success = insertvertex(
                m, b, newvertex, &badotri, (struct osub *)NULL, 1, 1);
            if (success == SUCCESSFULVERTEX)
            {
                if (m->steinerleft > 0)
                {
                    m->steinerleft--;
                }
            }
            else if (success == ENCROACHINGVERTEX)
            {
                /* If the newly inserted vertex encroaches upon a subsegment, */
                /*   delete the new vertex.                                   */
                undovertex(m, b);

                vertexdealloc(m, newvertex);
            }
            else if (success == VIOLATINGVERTEX)
            {
                /* Failed to insert the new vertex, but some subsegment was */
                /*   marked as being encroached.                            */
                vertexdealloc(m, newvertex);
            }
            else
            {   /* success == DUPLICATEVERTEX */
                /* Couldn't insert the new vertex because a vertex is already
                 * there. */

                vertexdealloc(m, newvertex);
            }
        }
        if (errorflag)
        {

            printf(
                "This probably means that I am trying to refine triangles\n");
            printf(
                "  to a smaller size than can be accommodated by the finite\n");
            printf("  precision of floating point arithmetic.  (You can be\n");
            printf("  sure of this if I fail to terminate.)\n");
            precisionerror();
        }
    }
}

/*****************************************************************************/
/*                                                                           */
/*  enforcequality()   Remove all the encroached subsegments and bad         */
/*                     triangles from the triangulation.                     */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::enforcequality(struct mesh *m, struct behavior *b)
{
    struct badtriang *badtri;
    int i;

    /* Initialize the pool of encroached subsegments. */
    poolinit(&m->badsubsegs,
             sizeof(struct badsubseg),
             BADSUBSEGPERBLOCK,
             BADSUBSEGPERBLOCK,
             0);

    /* Test all segments to see if they're encroached. */
    tallyencs(m, b);

    /* Fix encroached subsegments without noting bad triangles. */
    splitencsegs(m, b, 0);
    /* At this point, if we haven't run out of Steiner points, the */
    /*   triangulation should be (conforming) Delaunay.            */

    /* Next, we worry about enforcing triangle quality. */
    if ((b->minangle > 0.0) || b->usertest)
    {
        /* Initialize the pool of bad triangles. */
        poolinit(&m->badtriangles,
                 sizeof(struct badtriang),
                 BADTRIPERBLOCK,
                 BADTRIPERBLOCK,
                 0);
        /* Initialize the queues of bad triangles. */
        for (i = 0; i < 4096; i++)
        {
            m->queuefront[i] = (struct badtriang *)NULL;
        }
        m->firstnonemptyq = -1;
        /* Test all triangles to see if they're bad. */
        tallyfaces(m, b);
        /* Initialize the pool of recently flipped triangles. */
        poolinit(&m->flipstackers,
                 sizeof(struct flipstacker),
                 FLIPSTACKERPERBLOCK,
                 FLIPSTACKERPERBLOCK,
                 0);
        m->checkquality = 1;

        while ((m->badtriangles.items > 0) && (m->steinerleft != 0))
        {
            /* Fix one bad triangle by inserting a vertex at its circumcenter.
             */
            badtri = dequeuebadtriang(m);
            splittriangle(m, b, badtri);
            if (m->badsubsegs.items > 0)
            {
                /* Put bad triangle back in queue for another try later. */
                enqueuebadtriang(m, b, badtri);
                /* Fix any encroached subsegments that resulted. */
                /*   Record any new bad triangles that result.   */
                splitencsegs(m, b, 1);
            }
            else
            {
                /* Return the bad triangle to the pool. */
                pooldealloc(&m->badtriangles, (void *)badtri);
            }
        }
    }
    /* At this point, if the "-D" switch was selected and we haven't run out  */
    /*   of Steiner points, the triangulation should be (conforming) Delaunay */
    /*   and have no low-quality triangles.                                   */

}

/**                                                                         **/
/**                                                                         **/
/********* Mesh quality maintenance ends here                        *********/


/********* File I/O routines begin here                              *********/
/**                                                                         **/
/**                                                                         **/

/*****************************************************************************/
/*                                                                           */
/*  transfernodes()   Read the vertices from memory.                         */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::transfernodes(struct mesh *m,
                   struct behavior *b,
                   double *pointlist,
                   double *pointattriblist,
                   int *pointmarkerlist,
                   int numberofpoints,
                   int numberofpointattribs)
{
    vertex vertexloop;
    double x, y;
    int i, j;
    int coordindex;
    int attribindex;

    m->invertices   = numberofpoints;
    m->mesh_dim     = 2;
    m->nextras      = numberofpointattribs;
    m->readnodefile = 0;
    if (m->invertices < 3)
    {
        printf("Error:  Input must have at least three input vertices.\n");
        triexit(1);
    }
    if (m->nextras == 0)
    {
        b->weighted = 0;
    }

    initializevertexpool(m, b);

    /* Read the vertices. */
    coordindex  = 0;
    attribindex = 0;
    for (i = 0; i < m->invertices; i++)
    {
        vertexloop = (vertex)poolalloc(&m->vertices);
        /* Read the vertex coordinates. */
        x = vertexloop[0] = pointlist[coordindex++];
        y = vertexloop[1] = pointlist[coordindex++];
        /* Read the vertex attributes. */
        for (j = 0; j < numberofpointattribs; j++)
        {
            vertexloop[2 + j] = pointattriblist[attribindex++];
        }
        if (pointmarkerlist != (int *)NULL)
        {
            /* Read a vertex marker. */
            setvertexmark(vertexloop, pointmarkerlist[i]);
        }
        else
        {
            /* If no markers are specified, they default to zero. */
            setvertexmark(vertexloop, 0);
        }
        setvertextype(vertexloop, INPUTVERTEX);
        /* Determine the smallest and largest x and y coordinates. */
        if (i == 0)
        {
            m->xmin = m->xmax = x;
            m->ymin = m->ymax = y;
        }
        else
        {
            m->xmin = (x < m->xmin) ? x : m->xmin;
            m->xmax = (x > m->xmax) ? x : m->xmax;
            m->ymin = (y < m->ymin) ? y : m->ymin;
            m->ymax = (y > m->ymax) ? y : m->ymax;
        }
    }

    /* Nonexistent x value used as a flag to mark circle events in sweepline */
    /*   Delaunay algorithm.                                                 */
    m->xminextreme = 10 * m->xmin - 9 * m->xmax;
}

/*****************************************************************************/
/*                                                                           */
/*  writenodes()   Number the vertices and write them to a .node file.       */
/*                                                                           */
/*  To save memory, the vertex numbers are written over the boundary markers */
/*  after the vertices are written to a file.                                */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::writenodes(struct mesh *m,
                struct behavior *b,
                double **pointlist,
                double **pointattriblist,
                int **pointmarkerlist)
{
    double *plist;
    double *palist;
    int *pmlist;
    int coordindex;
    int attribindex;
    vertex vertexloop;
    long outvertices;
    int vertexnumber;
    int i;

    if (b->jettison)
    {
        outvertices = m->vertices.items - m->undeads;
    }
    else
    {
        outvertices = m->vertices.items;
    }

    /* Allocate memory for output vertices if necessary. */
    if (*pointlist == (double *)NULL)
    {
        *pointlist =
            (double *)trimalloc((int)(outvertices * 2 * sizeof(double)));
    }
    /* Allocate memory for output vertex attributes if necessary. */
    if ((m->nextras > 0) && (*pointattriblist == (double *)NULL))
    {
        *pointattriblist = (double *)trimalloc(
            (int)(outvertices * m->nextras * sizeof(double)));
    }
    /* Allocate memory for output vertex markers if necessary. */
    if (*pointmarkerlist == (int *)NULL)
    {
        *pointmarkerlist = (int *)trimalloc((int)(outvertices * sizeof(int)));
    }
    plist       = *pointlist;
    palist      = *pointattriblist;
    pmlist      = *pointmarkerlist;
    coordindex  = 0;
    attribindex = 0;

    traversalinit(&m->vertices);
    vertexnumber = 0;
    vertexloop   = vertextraverse(m);
    while (vertexloop != (vertex)NULL)
    {
        if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX))
        {
            /* X and y coordinates. */
            plist[coordindex++] = vertexloop[0];
            plist[coordindex++] = vertexloop[1];
            /* Vertex attributes. */
            for (i = 0; i < m->nextras; i++)
            {
                palist[attribindex++] = vertexloop[2 + i];
            }

            /* Copy the boundary marker. */
            pmlist[vertexnumber] = vertexmark(vertexloop);

            setvertexmark(vertexloop, vertexnumber);
            vertexnumber++;
        }
        vertexloop = vertextraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  numbernodes()   Number the vertices.                                     */
/*                                                                           */
/*  Each vertex is assigned a marker equal to its number.                    */
/*                                                                           */
/*  Used when writenodes() is not called because no .node file is written.   */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::numbernodes(struct mesh *m, struct behavior *b)
{
    vertex vertexloop;
    int vertexnumber;

    traversalinit(&m->vertices);
    vertexnumber = 0;
    vertexloop   = vertextraverse(m);
    while (vertexloop != (vertex)NULL)
    {
        setvertexmark(vertexloop, vertexnumber);
        if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX))
        {
            vertexnumber++;
        }
        vertexloop = vertextraverse(m);
    }
}

/*****************************************************************************/
/*                                                                           */
/*  writeelements()   Write the triangles to an .ele file.                   */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::writeelements(struct mesh *m,
                   struct behavior *b,
                   int **trianglelist,
                   double **triangleattriblist)
{
    int *tlist;
    double *talist;
    int vertexindex;
    int attribindex;
    struct otri triangleloop;
    vertex p1, p2, p3;
    long elementnumber;
    int i;

    /* Allocate memory for output triangles if necessary. */
    if (*trianglelist == (int *)NULL)
    {
        *trianglelist = (int *)trimalloc(
            (int)(m->triangles.items * 3 *
                  sizeof(int)));
    }
    /* Allocate memory for output triangle attributes if necessary. */
    if ((m->eextras > 0) && (*triangleattriblist == (double *)NULL))
    {
        *triangleattriblist = (double *)trimalloc(
            (int)(m->triangles.items * m->eextras * sizeof(double)));
    }
    tlist       = *trianglelist;
    talist      = *triangleattriblist;
    vertexindex = 0;
    attribindex = 0;

    traversalinit(&m->triangles);
    triangleloop.tri    = triangletraverse(m);
    triangleloop.orient = 0;
    elementnumber       = 0;
    while (triangleloop.tri != (triangle *)NULL)
    {
        org(triangleloop, p1);
        dest(triangleloop, p2);
        apex(triangleloop, p3);

            tlist[vertexindex++] = vertexmark(p1);
            tlist[vertexindex++] = vertexmark(p2);
            tlist[vertexindex++] = vertexmark(p3);



        for (i = 0; i < m->eextras; i++)
        {
            talist[attribindex++] = elemattribute(triangleloop, i);
        }

        triangleloop.tri = triangletraverse(m);
        elementnumber++;
    }
}

/*****************************************************************************/
/*                                                                           */
/*  writepoly()   Write the segments and holes to a .poly file.              */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::writepoly(struct mesh *m,
               struct behavior *b,
               int **segmentlist,
               int **segmentmarkerlist)
{
    int *slist;
    int *smlist;
    int index;
    struct osub subsegloop;
    vertex endpoint1, endpoint2;
    long subsegnumber;

    /* Allocate memory for output segments if necessary. */
    if (*segmentlist == (int *)NULL)
    {
        *segmentlist =
            (int *)trimalloc((int)(m->subsegs.items * 2 * sizeof(int)));
    }
    /* Allocate memory for output segment markers if necessary. */
    if (*segmentmarkerlist == (int *)NULL)
    {
        *segmentmarkerlist =
            (int *)trimalloc((int)(m->subsegs.items * sizeof(int)));
    }
    slist  = *segmentlist;
    smlist = *segmentmarkerlist;
    index  = 0;

    traversalinit(&m->subsegs);
    subsegloop.ss       = subsegtraverse(m);
    subsegloop.ssorient = 0;
    subsegnumber        = 0;
    while (subsegloop.ss != (subseg *)NULL)
    {
        sorg(subsegloop, endpoint1);
        sdest(subsegloop, endpoint2);
        /* Copy indices of the segment's two endpoints. */
        slist[index++] = vertexmark(endpoint1);
        slist[index++] = vertexmark(endpoint2);

        /* Copy the boundary marker. */
        smlist[subsegnumber] = mark(subsegloop);

        subsegloop.ss = subsegtraverse(m);
        subsegnumber++;
    }
}

/**                                                                         **/
/**                                                                         **/
/********* File I/O routines end here                                *********/



/*****************************************************************************/
/*                                                                           */
/*  main() or triangulate()   Gosh, do everything.                           */
/*                                                                           */
/*  The sequence is roughly as follows.  Many of these steps can be skipped, */
/*  depending on the command line switches.                                  */
/*                                                                           */
/*  - Initialize constants and parse the command line.                       */
/*  - Read the vertices from a file and either                               */
/*    - triangulate them (no -r), or                                         */
/*    - read an old mesh from files and reconstruct it (-r).                 */
/*  - Insert the PSLG segments (-p), and possibly segments on the convex     */
/*      hull (-c).                                                           */
/*  - Read the holes (-p), regional attributes (-pA), and regional area      */
/*      constraints (-pa).  Carve the holes and concavities, and spread the  */
/*      regional attributes and area constraints.                            */
/*  - Enforce the constraints on minimum angle (-q) and maximum area (-a).   */
/*      Also enforce the conforming Delaunay property (-q and -a).           */
/*  - Compute the number of edges in the resulting mesh.                     */
/*  - Promote the mesh's linear triangles to higher order elements (-o).     */
/*  - Write the output files and print the statistics.                       */
/*  - Check the consistency and Delaunay property of the mesh (-C).          */
/*                                                                           */
/*****************************************************************************/

void DelaunayTriangle::triangulate(char *triswitches)
{
    struct mesh m;
    struct behavior b;
    double *holearray; /* Array of holes. */
    double
        *regionarray; /* Array of regional attributes and area constraints. */

    plus1mod3[0]= 1;
    plus1mod3[1]= 2;
    plus1mod3[2]= 0;
    minus1mod3[0] = 2;
    minus1mod3[1] = 0;
    minus1mod3[2] = 1;

    triangleinit(&m);

    parsecommandline(1, &triswitches, &b);

    m.steinerleft = -1;

    transfernodes(&m,
                  &b,
                  in.pointlist,
                  in.pointattributelist,
                  in.pointmarkerlist,
                  in.numberofpoints,
                  in.numberofpointattributes);


    m.hullsize = delaunay(&m, &b); /* Triangulate the vertices. */

    /* Ensure that no vertex can be mistaken for a triangular bounding */
    /*   box vertex in insertvertex().                                 */
    m.infvertex1 = (vertex)NULL;
    m.infvertex2 = (vertex)NULL;
    m.infvertex3 = (vertex)NULL;

    if (b.usesegments)
    {
        m.checksegments = 1; /* Segments will be introduced next. */

        formskeleton(&m,
                     &b,
                     in.segmentlist,
                     in.segmentmarkerlist,
                     in.numberofsegments);
    }

    if (b.poly && (m.triangles.items > 0))
    {
        holearray   = in.holelist;
        m.holes     = in.numberofholes;
        regionarray = in.regionlist;
        m.regions   = in.numberofregions;

        /* Carve out holes and concavities. */
        carveholes(&m, &b, holearray, m.holes, regionarray, m.regions);
    }
    else
    {
        /* Without a PSLG, there can be no holes or regional attributes   */
        /*   or area constraints.  The following are set to zero to avoid */
        /*   an accidental free() later.                                  */
        m.holes   = 0;
        m.regions = 0;
    }

    if (b.quality && (m.triangles.items > 0))
    {
        enforcequality(&m, &b); /* Enforce angle and area constraints. */
    }

    /* Calculate the number of edges. */
    m.edges = (3l * m.triangles.items + m.hullsize) / 2l;

    if (b.jettison)
    {
        out.numberofpoints = m.vertices.items - m.undeads;
    }
    else
    {
        out.numberofpoints = m.vertices.items;
    }
    out.numberofpointattributes    = m.nextras;
    out.numberoftriangles          = m.triangles.items;
    out.numberofcorners            = 3;
    out.numberoftriangleattributes = m.eextras;
    out.numberofedges              = m.edges;
    if (b.usesegments)
    {
        out.numberofsegments = m.subsegs.items;
    }
    else
    {
        out.numberofsegments = m.hullsize;
    }

    /* writenodes() numbers the vertices too. */
    writenodes(&m,
               &b,
               &out.pointlist,
               &out.pointattributelist,
               &out.pointmarkerlist);

    writeelements(&m, &b, &out.trianglelist, &out.triangleattributelist);

    /* The -c switch (convex switch) causes a PSLG to be written */
    /*   even if none was read.                                  */
    if (b.poly)
    {

        writepoly(&m, &b, &out.segmentlist, &out.segmentmarkerlist);
        out.numberofholes   = m.holes;
        out.numberofregions = m.regions;
        if (b.poly)
        {
            out.holelist   = in.holelist;
            out.regionlist = in.regionlist;
        }
        else
        {
            out.holelist   = (double *)NULL;
            out.regionlist = (double *)NULL;
        }
    }

    triangledeinit(&m, &b);
}

}
}
