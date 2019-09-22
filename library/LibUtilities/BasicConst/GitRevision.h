///////////////////////////////////////////////////////////////////////////////
//
// File GitRevision.h
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
// Description: Constants for git SHA1 hash and branch name
//
///////////////////////////////////////////////////////////////////////////////

#ifndef  GITREVISION_H
#define  GITREVISION_H

#include <string>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    namespace NekConstants
    {
        LIB_UTILITIES_EXPORT extern const std::string kGitSha1;
        LIB_UTILITIES_EXPORT extern const std::string kGitBranch;
    }

    //This class is a workaround for a windows quirk which means that it cant
    //figure out how the extern works when a library other than LibUtilities
    //wants accsess to the information. This class wraps the consts with a class
    //so that they can be used elsewhere (such as nekmesh)
    namespace LibUtilities
    {
        class GitConsts
        {
        public:
            LIB_UTILITIES_EXPORT GitConsts(){}
            LIB_UTILITIES_EXPORT ~GitConsts(){}
            LIB_UTILITIES_EXPORT std::string GetSha1(){return NekConstants::kGitSha1;}
            LIB_UTILITIES_EXPORT std::string GetBranch(){return NekConstants::kGitBranch;}
        };
    }
}

#endif
