///////////////////////////////////////////////////////////////////////////////
//
// File: Misc.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Miscellaneous Time Integration routines.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Core/TimeIntegrationUtils.h>

namespace Nektar {
namespace SolverUtils {

void ParseTimeIntegrationParameters( std::string &sMethod,
                                     std::string &sVariant,
                                     std::string sParameter,
                                     std::vector< NekDouble > &params )
{
    // The Fractional-in-time uses a shorted class name.
    if( sMethod.find( "FractionalInTime" ) == 0 ) {
      sMethod = "FractionalIn";
    }

    // At this time nothing is done with the vairant.
    boost::ignore_unused(sVariant);
    
    // Parse the free parameters which are in a string.
    while( sParameter.size() )
    {
        size_t found = sParameter.find(" ");
        
        if( found == 0 )
        {
            sParameter = sParameter.substr( found+1 );
        }
        else if( found == std::string::npos )
        {
            if( sParameter.size() )
            {
                params.push_back( stoi( sParameter ) );
            }
          
            break;
        }
        else if( found != std::string::npos )
        {
            int fp = stoi( sParameter.substr(0, found) );
            
            params.push_back(fp);
            
            sParameter = sParameter.substr(found+1);
        }
    }
}

} // end namespace SolverUtils
} // end namespace Nektar
