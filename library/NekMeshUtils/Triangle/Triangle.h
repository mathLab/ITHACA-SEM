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

void triangulate(char *,
                 struct triangulateio *,
                 struct triangulateio *);
void trifree(void *memptr);
