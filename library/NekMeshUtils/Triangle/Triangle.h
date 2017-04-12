////////////////////////////////////////////////////////////////////////////////
//
//  File: Triangle.h
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
/*  (triangle.h)                                                             */
/*                                                                           */
/*  Include file for programs that call Triangle.                            */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 28, 2005                                                            */
/*                                                                           */
/*  Copyright 1996, 2005                                                     */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/
#ifndef NEKTAR_MESHUTILS_TRIANGLE_DT_H
#define NEKTAR_MESHUTILS_TRIANGLE_DT_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
namespace NekMeshUtils
{

/* For efficiency, a variety of data structures are allocated in bulk.  The  */
/*   following constants determine how many of each structure is allocated   */
/*   at once.                                                                */

#define TRIPERBLOCK 4092    /* Number of triangles allocated at once. */
#define SUBSEGPERBLOCK 508  /* Number of subsegments allocated at once. */
#define VERTEXPERBLOCK 4092 /* Number of vertices allocated at once. */
#define VIRUSPERBLOCK 1020  /* Number of virus triangles allocated at once. */
/* Number of encroached subsegments allocated at once. */
#define BADSUBSEGPERBLOCK 252
/* Number of skinny triangles allocated at once. */
#define BADTRIPERBLOCK 4092
/* Number of flipped triangles allocated at once. */
#define FLIPSTACKERPERBLOCK 252
/* Number of splay tree nodes allocated at once. */
#define SPLAYNODEPERBLOCK 508

/* The vertex types.   A DEADVERTEX has been deleted entirely.  An           */
/*   UNDEADVERTEX is not part of the mesh, but is written to the output      */
/*   .node file and affects the node indexing in the other output files.     */

#define INPUTVERTEX 0
#define SEGMENTVERTEX 1
#define FREEVERTEX 2
#define DEADVERTEX -32768
#define UNDEADVERTEX -32767

/* Two constants for algorithms based on random sampling.  Both constants    */
/*   have been chosen empirically to optimize their respective algorithms.   */

/* Used for the point location scheme of Mucke, Saias, and Zhu, to decide    */
/*   how large a random sample of triangles to inspect.                      */

#define SAMPLEFACTOR 11

/* Used in Fortune's sweepline Delaunay algorithm to determine what fraction */
/*   of boundary edges should be maintained in the splay tree for point      */
/*   location on the front.                                                  */

#define SAMPLERATE 10

/* A number that speaks for itself, every kissable digit.                    */

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

/* Another fave.                                                             */

#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732

/* And here's one for those of you who are intimidated by math.              */

#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333

/********* Mesh manipulation primitives begin here                   *********/
/**                                                                         **/
/**                                                                         **/

/********* Primitives for triangles                                  *********/
/*                                                                           */
/*                                                                           */

/* decode() converts a pointer to an oriented triangle.  The orientation is  */
/*   extracted from the two least significant bits of the pointer.           */

#define decode(ptr, otri)                                                      \
    (otri).orient = (int)((unsigned long)(ptr) & (unsigned long)3l);           \
    (otri).tri =                                                               \
        (triangle *)((unsigned long)(ptr) ^ (unsigned long)(otri).orient)

/* encode() compresses an oriented triangle into a single pointer.  It       */
/*   relies on the assumption that all triangles are aligned to four-byte    */
/*   boundaries, so the two least significant bits of (otri).tri are zero.   */

#define encode(otri)                                                           \
    (triangle)((unsigned long)(otri).tri | (unsigned long)(otri).orient)

/* The following handle manipulation primitives are all described by Guibas  */
/*   and Stolfi.  However, Guibas and Stolfi use an edge-based data          */
/*   structure, whereas I use a triangle-based data structure.               */

/* sym() finds the abutting triangle, on the same edge.  Note that the edge  */
/*   direction is necessarily reversed, because the handle specified by an   */
/*   oriented triangle is directed counterclockwise around the triangle.     */

#define sym(otri1, otri2)                                                      \
    ptr = (otri1).tri[(otri1).orient];                                         \
    decode(ptr, otri2);

#define symself(otri)                                                          \
    ptr = (otri).tri[(otri).orient];                                           \
    decode(ptr, otri);

/* lnext() finds the next edge (counterclockwise) of a triangle.             */

#define lnext(otri1, otri2)                                                    \
    (otri2).tri    = (otri1).tri;                                              \
    (otri2).orient = plus1mod3[(otri1).orient]

#define lnextself(otri) (otri).orient = plus1mod3[(otri).orient]

/* lprev() finds the previous edge (clockwise) of a triangle.                */

#define lprev(otri1, otri2)                                                    \
    (otri2).tri    = (otri1).tri;                                              \
    (otri2).orient = minus1mod3[(otri1).orient]

#define lprevself(otri) (otri).orient = minus1mod3[(otri).orient]

/* onext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same origin in the counterclockwise direction.  This */
/*   edge is part of a different triangle.                                   */

#define onext(otri1, otri2)                                                    \
    lprev(otri1, otri2);                                                       \
    symself(otri2);

#define onextself(otri)                                                        \
    lprevself(otri);                                                           \
    symself(otri);

/* oprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same origin in the clockwise direction.  This edge is part of  */
/*   a different triangle.                                                   */

#define oprev(otri1, otri2)                                                    \
    sym(otri1, otri2);                                                         \
    lnextself(otri2);

#define oprevself(otri)                                                        \
    symself(otri);                                                             \
    lnextself(otri);

/* dnext() spins counterclockwise around a vertex; that is, it finds the     */
/*   next edge with the same destination in the counterclockwise direction.  */
/*   This edge is part of a different triangle.                              */

#define dnext(otri1, otri2)                                                    \
    sym(otri1, otri2);                                                         \
    lprevself(otri2);

#define dnextself(otri)                                                        \
    symself(otri);                                                             \
    lprevself(otri);

/* dprev() spins clockwise around a vertex; that is, it finds the next edge  */
/*   with the same destination in the clockwise direction.  This edge is     */
/*   part of a different triangle.                                           */

#define dprev(otri1, otri2)                                                    \
    lnext(otri1, otri2);                                                       \
    symself(otri2);

#define dprevself(otri)                                                        \
    lnextself(otri);                                                           \
    symself(otri);

/* rnext() moves one edge counterclockwise about the adjacent triangle.      */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rnext(otri1, otri2)                                                    \
    sym(otri1, otri2);                                                         \
    lnextself(otri2);                                                          \
    symself(otri2);

#define rnextself(otri)                                                        \
    symself(otri);                                                             \
    lnextself(otri);                                                           \
    symself(otri);

/* rprev() moves one edge clockwise about the adjacent triangle.             */
/*   (It's best understood by reading Guibas and Stolfi.  It involves        */
/*   changing triangles twice.)                                              */

#define rprev(otri1, otri2)                                                    \
    sym(otri1, otri2);                                                         \
    lprevself(otri2);                                                          \
    symself(otri2);

#define rprevself(otri)                                                        \
    symself(otri);                                                             \
    lprevself(otri);                                                           \
    symself(otri);

/* These primitives determine or set the origin, destination, or apex of a   */
/* triangle.                                                                 */

#define org(otri, vertexptr)                                                   \
    vertexptr = (vertex)(otri).tri[plus1mod3[(otri).orient] + 3]

#define dest(otri, vertexptr)                                                  \
    vertexptr = (vertex)(otri).tri[minus1mod3[(otri).orient] + 3]

#define apex(otri, vertexptr) vertexptr = (vertex)(otri).tri[(otri).orient + 3]

#define setorg(otri, vertexptr)                                                \
    (otri).tri[plus1mod3[(otri).orient] + 3] = (triangle)vertexptr

#define setdest(otri, vertexptr)                                               \
    (otri).tri[minus1mod3[(otri).orient] + 3] = (triangle)vertexptr

#define setapex(otri, vertexptr)                                               \
    (otri).tri[(otri).orient + 3] = (triangle)vertexptr

/* Bond two triangles together.                                              */

#define bond(otri1, otri2)                                                     \
    (otri1).tri[(otri1).orient] = encode(otri2);                               \
    (otri2).tri[(otri2).orient] = encode(otri1)

/* Dissolve a bond (from one side).  Note that the other triangle will still */
/*   think it's connected to this triangle.  Usually, however, the other     */
/*   triangle is being deleted entirely, or bonded to another triangle, so   */
/*   it doesn't matter.                                                      */

#define dissolve(otri) (otri).tri[(otri).orient] = (triangle)m->dummytri

/* Copy an oriented triangle.                                                */

#define otricopy(otri1, otri2)                                                 \
    (otri2).tri    = (otri1).tri;                                              \
    (otri2).orient = (otri1).orient

/* Test for equality of oriented triangles.                                  */

#define otriequal(otri1, otri2)                                                \
    (((otri1).tri == (otri2).tri) && ((otri1).orient == (otri2).orient))

/* Primitives to infect or cure a triangle with the virus.  These rely on    */
/*   the assumption that all subsegments are aligned to four-byte boundaries.*/

#define infect(otri)                                                           \
    (otri).tri[6] = (triangle)((unsigned long)(otri).tri[6] | (unsigned long)2l)

#define uninfect(otri)                                                         \
    (otri).tri[6] =                                                            \
        (triangle)((unsigned long)(otri).tri[6] & ~(unsigned long)2l)

/* Test a triangle for viral infection.                                      */

#define infected(otri)                                                         \
    (((unsigned long)(otri).tri[6] & (unsigned long)2l) != 0l)

/* Check or set a triangle's attributes.                                     */

#define elemattribute(otri, attnum)                                            \
    ((double *)(otri).tri)[m->elemattribindex + (attnum)]

#define setelemattribute(otri, attnum, value)                                  \
    ((double *)(otri).tri)[m->elemattribindex + (attnum)] = value

/* Check or set a triangle's maximum area bound.                             */

#define areabound(otri) ((double *)(otri).tri)[m->areaboundindex]

#define setareabound(otri, value)                                              \
    ((double *)(otri).tri)[m->areaboundindex] = value

/* Check or set a triangle's deallocation.  Its second pointer is set to     */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its fourth pointer (its first vertex)    */
/*   is set to NULL in case a `badtriang' structure points to it.            */

#define deadtri(tria) ((tria)[1] == (triangle)NULL)

#define killtri(tria)                                                          \
    (tria)[1] = (triangle)NULL;                                                \
    (tria)[3] = (triangle)NULL

/********* Primitives for subsegments                                *********/
/*                                                                           */
/*                                                                           */

/* sdecode() converts a pointer to an oriented subsegment.  The orientation  */
/*   is extracted from the least significant bit of the pointer.  The two    */
/*   least significant bits (one for orientation, one for viral infection)   */
/*   are masked out to produce the double pointer. */

#define sdecode(sptr, osub)                                                    \
    (osub).ssorient = (int)((unsigned long)(sptr) & (unsigned long)1l);        \
    (osub).ss       = (subseg *)((unsigned long)(sptr) & ~(unsigned long)3l)

/* sencode() compresses an oriented subsegment into a single pointer.  It    */
/*   relies on the assumption that all subsegments are aligned to two-byte   */
/*   boundaries, so the least significant bit of (osub).ss is zero.          */

#define sencode(osub)                                                          \
    (subseg)((unsigned long)(osub).ss | (unsigned long)(osub).ssorient)

/* ssym() toggles the orientation of a subsegment.                           */

#define ssym(osub1, osub2)                                                     \
    (osub2).ss       = (osub1).ss;                                             \
    (osub2).ssorient = 1 - (osub1).ssorient

#define ssymself(osub) (osub).ssorient = 1 - (osub).ssorient

/* spivot() finds the other subsegment (from the same segment) that shares   */
/*   the same origin.                                                        */

#define spivot(osub1, osub2)                                                   \
    sptr = (osub1).ss[(osub1).ssorient];                                       \
    sdecode(sptr, osub2)

#define spivotself(osub)                                                       \
    sptr = (osub).ss[(osub).ssorient];                                         \
    sdecode(sptr, osub)

/* snext() finds the next subsegment (from the same segment) in sequence;    */
/*   one whose origin is the input subsegment's destination.                 */

#define snext(osub1, osub2)                                                    \
    sptr = (osub1).ss[1 - (osub1).ssorient];                                   \
    sdecode(sptr, osub2)

#define snextself(osub)                                                        \
    sptr = (osub).ss[1 - (osub).ssorient];                                     \
    sdecode(sptr, osub)

/* These primitives determine or set the origin or destination of a          */
/*   subsegment or the segment that includes it.                             */

#define sorg(osub, vertexptr) vertexptr = (vertex)(osub).ss[2 + (osub).ssorient]

#define sdest(osub, vertexptr)                                                 \
    vertexptr = (vertex)(osub).ss[3 - (osub).ssorient]

#define setsorg(osub, vertexptr)                                               \
    (osub).ss[2 + (osub).ssorient] = (subseg)vertexptr

#define setsdest(osub, vertexptr)                                              \
    (osub).ss[3 - (osub).ssorient] = (subseg)vertexptr

#define segorg(osub, vertexptr)                                                \
    vertexptr = (vertex)(osub).ss[4 + (osub).ssorient]

#define segdest(osub, vertexptr)                                               \
    vertexptr = (vertex)(osub).ss[5 - (osub).ssorient]

#define setsegorg(osub, vertexptr)                                             \
    (osub).ss[4 + (osub).ssorient] = (subseg)vertexptr

#define setsegdest(osub, vertexptr)                                            \
    (osub).ss[5 - (osub).ssorient] = (subseg)vertexptr

/* These primitives read or set a boundary marker.  Boundary markers are     */
/*   used to hold user-defined tags for setting boundary conditions in       */
/*   finite element solvers.                                                 */

#define mark(osub) (*(int *)((osub).ss + 8))

#define setmark(osub, value) *(int *)((osub).ss + 8) = value

/* Bond two subsegments together.                                            */

#define sbond(osub1, osub2)                                                    \
    (osub1).ss[(osub1).ssorient] = sencode(osub2);                             \
    (osub2).ss[(osub2).ssorient] = sencode(osub1)

/* Dissolve a subsegment bond (from one side).  Note that the other          */
/*   subsegment will still think it's connected to this subsegment.          */

#define sdissolve(osub) (osub).ss[(osub).ssorient] = (subseg)m->dummysub

/* Copy a subsegment.                                                        */

#define subsegcopy(osub1, osub2)                                               \
    (osub2).ss       = (osub1).ss;                                             \
    (osub2).ssorient = (osub1).ssorient

/* Test for equality of subsegments.                                         */

#define subsegequal(osub1, osub2)                                              \
    (((osub1).ss == (osub2).ss) && ((osub1).ssorient == (osub2).ssorient))

/* Check or set a subsegment's deallocation.  Its second pointer is set to   */
/*   NULL to indicate that it is not allocated.  (Its first pointer is used  */
/*   for the stack of dead items.)  Its third pointer (its first vertex)     */
/*   is set to NULL in case a `badsubseg' structure points to it.            */

#define deadsubseg(sub) ((sub)[1] == (subseg)NULL)

#define killsubseg(sub)                                                        \
    (sub)[1] = (subseg)NULL;                                                   \
    (sub)[2] = (subseg)NULL

/********* Primitives for interacting triangles and subsegments      *********/
/*                                                                           */
/*                                                                           */

/* tspivot() finds a subsegment abutting a triangle.                         */

#define tspivot(otri, osub)                                                    \
    sptr = (subseg)(otri).tri[6 + (otri).orient];                              \
    sdecode(sptr, osub)

/* stpivot() finds a triangle abutting a subsegment.  It requires that the   */
/*   variable `ptr' of type `triangle' be defined.                           */

#define stpivot(osub, otri)                                                    \
    ptr = (triangle)(osub).ss[6 + (osub).ssorient];                            \
    decode(ptr, otri)

/* Bond a triangle to a subsegment.                                          */

#define tsbond(otri, osub)                                                     \
    (otri).tri[6 + (otri).orient]  = (triangle)sencode(osub);                  \
    (osub).ss[6 + (osub).ssorient] = (subseg)encode(otri)

/* Dissolve a bond (from the triangle side).                                 */

#define tsdissolve(otri) (otri).tri[6 + (otri).orient] = (triangle)m->dummysub

/* Dissolve a bond (from the subsegment side).                               */

#define stdissolve(osub) (osub).ss[6 + (osub).ssorient] = (subseg)m->dummytri

/********* Primitives for vertices                                   *********/
/*                                                                           */
/*                                                                           */

#define vertexmark(vx) ((int *)(vx))[m->vertexmarkindex]

#define setvertexmark(vx, value) ((int *)(vx))[m->vertexmarkindex] = value

#define vertextype(vx) ((int *)(vx))[m->vertexmarkindex + 1]

#define setvertextype(vx, value) ((int *)(vx))[m->vertexmarkindex + 1] = value

#define vertex2tri(vx) ((triangle *)(vx))[m->vertex2triindex]

#define setvertex2tri(vx, value) ((triangle *)(vx))[m->vertex2triindex] = value

/**                                                                         **/
/**                                                                         **/
/********* Mesh manipulation primitives end here                     *********/

struct triangulateio
{
    double *pointlist;             /* In / out */
    double *pointattributelist;    /* In / out */
    int *pointmarkerlist;        /* In / out */
    int numberofpoints;          /* In / out */
    int numberofpointattributes; /* In / out */

    int *trianglelist;              /* In / out */
    double *triangleattributelist;    /* In / out */
    double *trianglearealist;         /* In only */
    int *neighborlist;              /* Out only */
    int numberoftriangles;          /* In / out */
    int numberofcorners;            /* In / out */
    int numberoftriangleattributes; /* In / out */

    int *segmentlist;       /* In / out */
    int *segmentmarkerlist; /* In / out */
    int numberofsegments;   /* In / out */

    double *holelist;    /* In / pointer to array copied out */
    int numberofholes; /* In / copied out */

    double *regionlist;    /* In / pointer to array copied out */
    int numberofregions; /* In / copied out */

    int *edgelist;       /* Out only */
    int *edgemarkerlist; /* Not used with Voronoi diagram; out only */
    double *normlist;      /* Used only with Voronoi diagram; out only */
    int numberofedges;   /* Out only */
};

/* Labels that signify the result of point location.  The result of a        */
/*   search indicates that the point falls in the interior of a triangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */

enum locateresult
{
    INTRIANGLE,
    ONEDGE,
    ONVERTEX,
    OUTSIDE
};

/* Labels that signify the result of vertex insertion.  The result indicates */
/*   that the vertex was inserted with complete success, was inserted but    */
/*   encroaches upon a subsegment, was not inserted because it lies on a     */
/*   segment, or was not inserted because another vertex occupies the same   */
/*   location.                                                               */

enum insertvertexresult
{
    SUCCESSFULVERTEX,
    ENCROACHINGVERTEX,
    VIOLATINGVERTEX,
    DUPLICATEVERTEX
};

/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction triangle, along the left edge of the direction triangle,  */
/*   or along the right edge of the direction triangle.                      */

enum finddirectionresult
{
    WITHIN,
    LEFTCOLLINEAR,
    RIGHTCOLLINEAR
};

/*****************************************************************************/
/*                                                                           */
/*  The basic mesh data structures                                           */
/*                                                                           */
/*  There are three:  vertices, triangles, and subsegments (abbreviated      */
/*  `subseg').  These three data structures, linked by pointers, comprise    */
/*  the mesh.  A vertex simply represents a mesh vertex and its properties.  */
/*  A triangle is a triangle.  A subsegment is a special data structure used */
/*  to represent an impenetrable edge of the mesh (perhaps on the outer      */
/*  boundary, on the boundary of a hole, or part of an internal boundary     */
/*  separating two triangulated regions).  Subsegments represent boundaries, */
/*  defined by the user, that triangles may not lie across.                  */
/*                                                                           */
/*  A triangle consists of a list of three vertices, a list of three         */
/*  adjoining triangles, a list of three adjoining subsegments (when         */
/*  segments exist), an arbitrary number of optional user-defined            */
/*  floating-point attributes, and an optional area constraint.  The latter  */
/*  is an upper bound on the permissible area of each triangle in a region,  */
/*  used for mesh refinement.                                                */
/*                                                                           */
/*  For a triangle on a boundary of the mesh, some or all of the neighboring */
/*  triangles may not be present.  For a triangle in the interior of the     */
/*  mesh, often no neighboring subsegments are present.  Such absent         */
/*  triangles and subsegments are never represented by NULL pointers; they   */
/*  are represented by two special records:  `dummytri', the triangle that   */
/*  fills "outer space", and `dummysub', the omnipresent subsegment.         */
/*  `dummytri' and `dummysub' are used for several reasons; for instance,    */
/*  they can be dereferenced and their contents examined without violating   */
/*  protected memory.                                                        */
/*                                                                           */
/*  However, it is important to understand that a triangle includes other    */
/*  information as well.  The pointers to adjoining vertices, triangles, and */
/*  subsegments are ordered in a way that indicates their geometric relation */
/*  to each other.  Furthermore, each of these pointers contains orientation */
/*  information.  Each pointer to an adjoining triangle indicates which face */
/*  of that triangle is contacted.  Similarly, each pointer to an adjoining  */
/*  subsegment indicates which side of that subsegment is contacted, and how */
/*  the subsegment is oriented relative to the triangle.                     */
/*                                                                           */
/*  The data structure representing a subsegment may be thought to be        */
/*  abutting the edge of one or two triangle data structures:  either        */
/*  sandwiched between two triangles, or resting against one triangle on an  */
/*  exterior boundary or hole boundary.                                      */
/*                                                                           */
/*  A subsegment consists of a list of four vertices--the vertices of the    */
/*  subsegment, and the vertices of the segment it is a part of--a list of   */
/*  two adjoining subsegments, and a list of two adjoining triangles.  One   */
/*  of the two adjoining triangles may not be present (though there should   */
/*  always be one), and neighboring subsegments might not be present.        */
/*  Subsegments also store a user-defined integer "boundary marker".         */
/*  Typically, this integer is used to indicate what boundary conditions are */
/*  to be applied at that location in a finite element simulation.           */
/*                                                                           */
/*  Like triangles, subsegments maintain information about the relative      */
/*  orientation of neighboring objects.                                      */
/*                                                                           */
/*  Vertices are relatively simple.  A vertex is a list of floating-point    */
/*  numbers, starting with the x, and y coordinates, followed by an          */
/*  arbitrary number of optional user-defined floating-point attributes,     */
/*  followed by an integer boundary marker.  During the segment insertion    */
/*  phase, there is also a pointer from each vertex to a triangle that may   */
/*  contain it.  Each pointer is not always correct, but when one is, it     */
/*  speeds up segment insertion.  These pointers are assigned values once    */
/*  at the beginning of the segment insertion phase, and are not used or     */
/*  updated except during this phase.  Edge flipping during segment          */
/*  insertion will render some of them incorrect.  Hence, don't rely upon    */
/*  them for anything.                                                       */
/*                                                                           */
/*  Other than the exception mentioned above, vertices have no information   */
/*  about what triangles, subfacets, or subsegments they are linked to.      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  Handles                                                                  */
/*                                                                           */
/*  The oriented triangle (`otri') and oriented subsegment (`osub') data     */
/*  structures defined below do not themselves store any part of the mesh.   */
/*  The mesh itself is made of `triangle's, `subseg's, and `vertex's.        */
/*                                                                           */
/*  Oriented triangles and oriented subsegments will usually be referred to  */
/*  as "handles."  A handle is essentially a pointer into the mesh; it       */
/*  allows you to "hold" one particular part of the mesh.  Handles are used  */
/*  to specify the regions in which one is traversing and modifying the mesh.*/
/*  A single `triangle' may be held by many handles, or none at all.  (The   */
/*  latter case is not a memory leak, because the triangle is still          */
/*  connected to other triangles in the mesh.)                               */
/*                                                                           */
/*  An `otri' is a handle that holds a triangle.  It holds a specific edge   */
/*  of the triangle.  An `osub' is a handle that holds a subsegment.  It     */
/*  holds either the left or right side of the subsegment.                   */
/*                                                                           */
/*  Navigation about the mesh is accomplished through a set of mesh          */
/*  manipulation primitives, further below.  Many of these primitives take   */
/*  a handle and produce a new handle that holds the mesh near the first     */
/*  handle.  Other primitives take two handles and glue the corresponding    */
/*  parts of the mesh together.  The orientation of the handles is           */
/*  important.  For instance, when two triangles are glued together by the   */
/*  bond() primitive, they are glued at the edges on which the handles lie.  */
/*                                                                           */
/*  Because vertices have no information about which triangles they are      */
/*  attached to, I commonly represent a vertex by use of a handle whose      */
/*  origin is the vertex.  A single handle can simultaneously represent a    */
/*  triangle, an edge, and a vertex.                                         */
/*                                                                           */
/*****************************************************************************/

/* The triangle data structure.  Each triangle contains three pointers to    */
/*   adjoining triangles, plus three pointers to vertices, plus three        */
/*   pointers to subsegments (declared below; these pointers are usually     */
/*   `dummysub').  It may or may not also contain user-defined attributes    */
/*   and/or a floating-point "area constraint."  It may also contain extra   */
/*   pointers for nodes, when the user asks for high-order elements.         */
/*   Because the size and structure of a `triangle' is not decided until     */
/*   runtime, I haven't simply declared the type `triangle' as a struct.     */

typedef double **triangle; /* Really:  typedef triangle *triangle   */

/* An oriented triangle:  includes a pointer to a triangle and orientation.  */
/*   The orientation denotes an edge of the triangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge always points    */
/*   counterclockwise about the corresponding triangle.                      */

struct otri
{
    triangle *tri;
    int orient; /* Ranges from 0 to 2. */
};

/* The subsegment data structure.  Each subsegment contains two pointers to  */
/*   adjoining subsegments, plus four pointers to vertices, plus two         */
/*   pointers to adjoining triangles, plus one boundary marker, plus one     */
/*   segment number.                                                         */

typedef double **subseg; /* Really:  typedef subseg *subseg   */

/* An oriented subsegment:  includes a pointer to a subsegment and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */

struct osub
{
    subseg *ss;
    int ssorient; /* Ranges from 0 to 1. */
};

/* The vertex data structure.  Each vertex is actually an array of doubles. */
/*   The number of doubles is unknown until runtime.  An integer boundary */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   doubles. */

typedef double *vertex;

/* A queue used to store encroached subsegments.  Each subsegment's vertices */
/*   are stored so that we can check whether a subsegment is still the same. */

struct badsubseg
{
    subseg encsubseg;             /* An encroached subsegment. */
    vertex subsegorg, subsegdest; /* Its two vertices. */
};

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the triangle.  Each triangle's vertices are    */
/*   stored so that one can check whether a triangle is still the same.      */

struct badtriang
{
    triangle poortri; /* A skinny or too-large triangle. */
    double key;       /* cos^2 of smallest (apical) angle. */
    vertex triangorg, triangdest, triangapex; /* Its three vertices. */
    struct badtriang *nexttriang; /* Pointer to next bad triangle. */
};

/* A stack of triangles flipped during the most recent vertex insertion.     */
/*   The stack is used to undo the vertex insertion if the vertex encroaches */
/*   upon a subsegment.                                                      */

struct flipstacker
{
    triangle flippedtri;          /* A recently flipped triangle. */
    struct flipstacker *prevflip; /* Previous flip in the stack. */
};

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not point directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `vertex' or to a triangle (in encoded format, so     */
/*   that an orientation is included).  In the latter case, the origin of    */
/*   the oriented triangle is the apex of a "circle event" of the sweepline  */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */

struct event
{
    double xkey, ykey; /* Coordinates of the event. */
    void *eventptr;    /* Can be a vertex or the location of a circle event. */
    int heapposition;  /* Marks this event's position in the heap. */
};

/* A node in the splay tree.  Each node holds an oriented ghost triangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a triangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a triangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */

struct splaynode
{
    struct otri keyedge; /* Lprev of an edge on the front. */
    vertex keydest;      /* Used to verify that splay node is still live. */
    struct splaynode *lchild, *rchild; /* Children in splay tree. */
};

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* alignbytes determines how new records should be aligned in memory.        */
/*   itembytes is the length of a record in bytes (after rounding up).       */
/*   itemsperblock is the number of items allocated at once in a single      */
/*   block.  itemsfirstblock is the number of items in the first block,      */
/*   which can vary from the others.  items is the number of currently       */
/*   allocated items.  maxitems is the maximum number of items that have     */
/*   been allocated at once; it is the current number of items plus the      */
/*   number of records kept on deaditemstack.                                */

struct memorypool
{
    void **firstblock, **nowblock;
    void *nextitem;
    void *deaditemstack;
    void **pathblock;
    void *pathitem;
    int alignbytes;
    int itembytes;
    int itemsperblock;
    int itemsfirstblock;
    long items, maxitems;
    int unallocateditems;
    int pathitemsleft;
};

/* Mesh data structure.  Triangle operates on only one mesh, but the mesh    */
/*   structure is used (instead of global variables) to allow reentrancy.    */

struct mesh
{

    /* Variables used to allocate memory for triangles, subsegments, vertices,
     */
    /*   viri (triangles being eaten), encroached segments, bad (skinny or too
     */
    /*   large) triangles, and splay tree nodes. */

    struct memorypool triangles;
    struct memorypool subsegs;
    struct memorypool vertices;
    struct memorypool viri;
    struct memorypool badsubsegs;
    struct memorypool badtriangles;
    struct memorypool flipstackers;
    struct memorypool splaynodes;

    /* Variables that maintain the bad triangle queues.  The queues are */
    /*   ordered from 4095 (highest priority) to 0 (lowest priority). */

    struct badtriang *queuefront[4096];
    struct badtriang *queuetail[4096];
    int nextnonemptyq[4096];
    int firstnonemptyq;

    /* Variable that maintains the stack of recently flipped triangles. */

    struct flipstacker *lastflip;

    /* Other variables. */

    double xmin, xmax, ymin, ymax; /* x and y bounds. */
    double xminextreme; /* Nonexistent x value used as a flag in sweepline. */
    int invertices;     /* Number of input vertices. */
    int inelements;     /* Number of input triangles. */
    int insegments;     /* Number of input segments. */
    int holes;          /* Number of input holes. */
    int regions;        /* Number of input regions. */
    int undeads;   /* Number of input vertices that don't appear in the mesh. */
    long edges;    /* Number of output edges. */
    int mesh_dim;  /* Dimension (ought to be 2). */
    int nextras;   /* Number of attributes per vertex. */
    int eextras;   /* Number of attributes per triangle. */
    long hullsize; /* Number of edges in convex hull. */
    int steinerleft;     /* Number of Steiner points not yet used. */
    int vertexmarkindex; /* Index to find boundary marker of a vertex. */
    int vertex2triindex; /* Index to find a triangle adjacent to a vertex. */
    int highorderindex; /* Index to find extra nodes for high-order elements. */
    int elemattribindex; /* Index to find attributes of a triangle. */
    int areaboundindex;  /* Index to find area bound of a triangle. */
    int checksegments;   /* Are there segments in the triangulation yet? */
    int checkquality;    /* Has quality triangulation begun yet? */
    int readnodefile;    /* Has a .node file been read? */
    long samples;        /* Number of random samples for point location. */

    long incirclecount;     /* Number of incircle tests performed. */
    long counterclockcount; /* Number of counterclockwise tests performed. */
    long orient3dcount;     /* Number of 3D orientation tests performed. */
    long hyperbolacount;    /* Number of right-of-hyperbola tests performed. */
    long circumcentercount; /* Number of circumcenter calculations performed. */
    long circletopcount;    /* Number of circle top calculations performed. */

    /* Triangular bounding box vertices. */

    vertex infvertex1, infvertex2, infvertex3;

    /* Pointer to the `triangle' that occupies all of "outer space." */

    triangle *dummytri;
    triangle *dummytribase; /* Keep base address so we can free() it later. */

    /* Pointer to the omnipresent subsegment.  Referenced by any triangle or */
    /*   subsegment that isn't really connected to a subsegment at that */
    /*   location. */

    subseg *dummysub;
    subseg *dummysubbase; /* Keep base address so we can free() it later. */

    /* Pointer to a recently visited triangle.  Improves point location if */
    /*   proximate vertices are inserted sequentially. */

    struct otri recenttri;

}; /* End of `struct mesh'. */

/* Data structure for command line switches and file names.  This structure  */
/*   is used (instead of global variables) to allow reentrancy.              */

struct behavior
{

    /* Switches for the triangulator. */
    /*   poly: -p switch. */
    /*   quality: -q switch. */
    /*     minangle: minimum angle bound, specified after -q switch. */
    /*     goodangle: cosine squared of minangle. */
    /*     offconstant: constant used to place off-center Steiner points. */
    /*   usertest: -u switch. */ //pretend to use this one to make triangle do its own calcs
    /*   weighted: 1 for -w switch, 2 for -W switch.  jettison: -j switch */
    /*   nobisect: count of how often -Y switch is selected. */
    /*   usesegments: -p, -r, -q, or -c switch; determines whether segments are
     */
    /*     used at all. */
    /*                                                                           */
    /* Read the instructions to find out the meaning of these switches. */

    int poly, quality, usertest;
    int weighted, jettison;
    int usesegments;
    int nobisect;
    double minangle, goodangle, offconstant;

}; /* End of `struct behavior'. */

class DelaunayTriangle
{
public:
    friend class MemoryManager<DelaunayTriangle>;

    void triangulate(char *);
    void trifree(void *memptr);

    struct triangulateio in,out;

private:
    int triunsuitable(vertex triorg, vertex tridest, vertex triapex, double area);
    void triexit(int status);
    void *trimalloc(int size);
    void internalerror();
    void parsecommandline(int argc, char **argv, struct behavior *b);
    void poolzero(struct memorypool *pool);
    void poolrestart(struct memorypool *pool);
    void pooldeinit(struct memorypool *pool);
    void *poolalloc(struct memorypool *pool);
    void pooldealloc(struct memorypool *pool, void *dyingitem);
    void traversalinit(struct memorypool *pool);
    void *traverse(struct memorypool *pool);
    void initializevertexpool(struct mesh *m, struct behavior *b);
    void initializetrisubpools(struct mesh *m, struct behavior *b);
    void triangledealloc(struct mesh *m, triangle *dyingtriangle);
    triangle *triangletraverse(struct mesh *m);
    void subsegdealloc(struct mesh *m, subseg *dyingsubseg);
    subseg *subsegtraverse(struct mesh *m);
    void vertexdealloc(struct mesh *m, vertex dyingvertex);
    vertex vertextraverse(struct mesh *m);
    void badsubsegdealloc(struct mesh *m, struct badsubseg *dyingseg);
    struct badsubseg *badsubsegtraverse(struct mesh *m);
    vertex getvertex(struct mesh *m, struct behavior *b, int number);
    void triangledeinit(struct mesh *m, struct behavior *b);
    void maketriangle(struct mesh *m, struct behavior *b, struct otri *newotri);
    void makesubseg(struct mesh *m, struct osub *newsubseg);
    void exactinit();
    int scale_expansion_zeroelimTRI(int elen, double *e, double b, double *h);
    double estimateTRI(int elen, double *e);
    double counterclockwise(struct mesh *m, struct behavior *b, vertex pa, vertex pb, vertex pc);
    double counterclockwiseadapt(vertex pa, vertex pb, vertex pc, double detsum);
    void poolinit(struct memorypool *pool,  int bytecount,  int itemcount, int firstitemcount, int alignment);
    void dummyinit(struct mesh *m, struct behavior *b, int trianglebytes, int subsegbytes);
    double incircleadaptTRI(vertex pa, vertex pb, vertex pc, vertex pd, double permanent);
    double incircle(struct mesh *m, struct behavior *b, vertex pa, vertex pb, vertex pc, vertex pd);
    void triangleinit(struct mesh *m);
    unsigned long randomnation(unsigned int choices);
    struct badtriang *dequeuebadtriang(struct mesh *m);
    void testtriangle(struct mesh *m, struct behavior *b, struct otri *testtri);
    void makevertexmap(struct mesh *m, struct behavior *b);
    void flip(struct mesh *m, struct behavior *b, struct otri *flipedge);
    void unflip(struct mesh *m, struct behavior *b, struct otri *flipedge);
    int fast_expansion_sum_zeroelimTRI(int elen, double *e, int flen, double *f, double *h);
    double orient3dadapt(vertex pa, vertex pb, vertex pc, vertex pd, double aheight, double bheight, double cheight, double dheight, double permanent);
    double orient3d(struct mesh *m, struct behavior *b, vertex pa, vertex pb, vertex pc, vertex pd, double aheight, double bheight, double cheight, double dheight);
    double nonregular(struct mesh *m, struct behavior *b, vertex pa, vertex pb, vertex pc, vertex pd);
    void findcircumcenter(struct mesh *m, struct behavior *b, vertex torg, vertex tdest, vertex tapex, vertex circumcenter, double *xi, double *eta, int offcenter);
    void enqueuebadtriang(struct mesh *m, struct behavior *b, struct badtriang *badtri);
    void enqueuebadtri(struct mesh *m, struct behavior *b, struct otri *enqtri, double minedge, vertex enqapex, vertex enqorg, vertex enqdest);
    int checkseg4encroach(struct mesh *m, struct behavior *b, struct osub *testsubseg);
    enum locateresult preciselocate(struct mesh *m, struct behavior *b, vertex searchpoint, struct otri *searchtri, int stopatsubsegment);
    enum locateresult locate(struct mesh *m, struct behavior *b, vertex searchpoint, struct otri *searchtri);
    void insertsubseg(struct mesh *m, struct behavior *b, struct otri *tri, int subsegmark);
    enum insertvertexresult insertvertex(struct mesh *m, struct behavior *b, vertex newvertex, struct otri *searchtri, struct osub *splitseg, int segmentflaws, int triflaws);
    void triangulatepolygon(struct mesh *m, struct behavior *b, struct otri *firstedge, struct otri *lastedge, int edgecount, int doflip, int triflaws);
    void deletevertex(struct mesh *m, struct behavior *b, struct otri *deltri);
    void undovertex(struct mesh *m, struct behavior *b);
    void vertexsort(vertex *sortarray, unsigned int arraysize);
    void vertexmedian(vertex *sortarray, int arraysize, int median, int axis);
    void alternateaxes(vertex *sortarray, int arraysize, int axis);
    long removeghosts(struct mesh *m, struct behavior *b, struct otri *startghost);
    long divconqdelaunay(struct mesh *m, struct behavior *b);
    long delaunay(struct mesh *m, struct behavior *b);
    void mergehulls(struct mesh *m, struct behavior *b, struct otri *farleft, struct otri *innerleft, struct otri *innerright, struct otri *farright, int axis);
    void divconqrecurse(struct mesh *m, struct behavior *b, vertex *sortarray, int vertices, int axis, struct otri *farleft, struct otri *farright);
    enum finddirectionresult finddirection(struct mesh *m, struct behavior *b, struct otri *searchtri, vertex searchpoint);
    void segmentintersection(struct mesh *m,  struct behavior *b, struct otri *splittri, struct osub *splitsubseg, vertex endpoint2);
    int scoutsegment(struct mesh *m, struct behavior *b, struct otri *searchtri, vertex endpoint2, int newmark);
    void conformingedge(struct mesh *m, struct behavior *b, vertex endpoint1, vertex endpoint2, int newmark);
    void delaunayfixup(struct mesh *m, struct behavior *b, struct otri *fixuptri, int leftside);
    void constrainededge(struct mesh *m, struct behavior *b, struct otri *starttri, vertex endpoint2, int newmark);
    void insertsegment(struct mesh *m, struct behavior *b, vertex endpoint1, vertex endpoint2, int newmark);
    void markhull(struct mesh *m, struct behavior *b);
    void formskeleton(struct mesh *m, struct behavior *b, int *segmentlist, int *segmentmarkerlist, int numberofsegments);
    void infecthull(struct mesh *m, struct behavior *b);
    void plague(struct mesh *m, struct behavior *b);
    void regionplague(struct mesh *m, struct behavior *b, double attribute, double area);
    void carveholes(struct mesh *m, struct behavior *b, double *holelist, int holes, double *regionlist, int regions);
    void tallyencs(struct mesh *m, struct behavior *b);
    void precisionerror();
    void splitencsegs(struct mesh *m, struct behavior *b, int triflaws);
    void tallyfaces(struct mesh *m, struct behavior *b);
    void splittriangle(struct mesh *m, struct behavior *b, struct badtriang *badtri);
    void enforcequality(struct mesh *m, struct behavior *b);
    void transfernodes(struct mesh *m, struct behavior *b, double *pointlist, double *pointattriblist, int *pointmarkerlist, int numberofpoints, int numberofpointattribs);
    void writenodes(struct mesh *m, struct behavior *b, double **pointlist, double **pointattriblist, int **pointmarkerlist);
    void numbernodes(struct mesh *m, struct behavior *b);
    void writeelements(struct mesh *m, struct behavior *b, int **trianglelist, double **triangleattriblist);
    void writepoly(struct mesh *m, struct behavior *b, int **segmentlist, int **segmentmarkerlist);

    /* Global constants.                                                         */

    double splitter; /* Used to split double factors for exact multiplication. */
    double epsilon;  /* Floating-point machine epsilon. */
    double resulterrbound;
    double ccwerrboundA, ccwerrboundB, ccwerrboundC;
    double iccerrboundA, iccerrboundB, iccerrboundC;
    double o3derrboundA, o3derrboundB, o3derrboundC;

    /* Random number seed is not constant, but I've made it global anyway.       */

    unsigned long randomseed; /* Current random number seed. */

    /* Fast lookup arrays to speed some of the mesh manipulation primitives.     */

    int plus1mod3[3];
    int minus1mod3[3];
};

}
}
#endif
