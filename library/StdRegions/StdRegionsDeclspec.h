///////////////////////////////////////////////////////////////////////////////
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_STD_REGIONS_DECLSPEC_H
#define NEKTAR_STD_REGIONS_DECLSPEC_H

#if defined(_MSC_VER)
    #ifdef STD_REGIONS_EXPORTS
        #define STD_REGIONS_EXPORT _declspec(dllexport)
    #else
        #define STD_REGIONS_EXPORT _declspec(dllimport)
    #endif
#else
    #define STD_REGIONS_EXPORT
#endif

#endif //NEKTAR_STD_REGIONS_DECLSPEC_H

