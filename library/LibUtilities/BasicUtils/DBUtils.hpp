///////////////////////////////////////////////////////////////////////////////
//
// File DBUtils.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: output function for use in debugging
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_DBUTILS_HPP
#define NEKTAR_LIB_LIBUTILITIES_DBUTILS_HPP

#include <LibUtilities/BasicUtils/VmathArray.hpp> 

namespace DBUtils
{
    using namespace std;
    using namespace Nektar;

    const int StopDefault = -99;

    template<class T> void Output1DArray(const Array<OneD, const T> &in, const int start = 0,
                                         const int stop = StopDefault)
    {
        int i;
        
        ASSERTL1(start < in.num_elements(), "Start value is outside array range ");

        if(stop == StopDefault)
        {

            for(i = start; i < in.num_elements(); ++i)
            {
                cout << in[i] << endl;
            }
        }
        else
        {
            ASSERTL1(stop <= in.num_elements(), "Stop value is outside array range ");
            
            for(i = start; i < stop; ++i)
            {
                cout << in[i] << endl;
            }
        }
    }

    template<class T> void Output1DArray(const Array<OneD, const T> &in, std::string outfile, 
                                         const int start = 0,
                                         const int stop = StopDefault)
    {
        int i;
        
        ASSERTL1(start < in.num_elements(), "Start value is outside array range ");

        ofstream ofile(outfile.c_str());

        if(stop == StopDefault)
        {

            for(i = start; i < in.num_elements(); ++i)
            {
                ofile << in[i] << endl;
            }
        }
        else
        {
            ASSERTL1(stop <= in.num_elements(), "Stop value is outside array range ");
            
            for(i = start; i < stop; ++i)
            {
                ofile << in[i] << endl;
            }
        }

    }

    template<class T> void Output1DArray(const Array<OneD, const T> &in, ofstream &ofile, 
                                         const int start = 0,
                                         const int stop = StopDefault)
    {
        int i;
        
        ASSERTL1(start < in.num_elements(), "Start value is outside array range ");

        if(stop == StopDefault)
        {

            for(i = start; i < in.num_elements(); ++i)
            {
                ofile << in[i] << endl;
            }
        }
        else
        {
            ASSERTL1(stop <= in.num_elements(), "Stop value is outside array range ");
            
            for(i = start; i < stop; ++i)
            {
                ofile << in[i] << endl;
            }
        }

    }

    template<class T> void Output1DArray(const NekVector<T> &in, ofstream &ofile, 
                                         const int start = 0,
                                         const int stop = StopDefault)
    {
        int i;
        
        ASSERTL1(start < in.GetDimension(), "Start value is outside array range ");

        if(stop == StopDefault)
        {

            for(i = start; i < in.GetDimension(); ++i)
            {
                ofile << in[i] << endl;
            }
        }
        else
        {
            ASSERTL1(stop <= in.GetDimension(), "Stop value is outside array range ");
            
            for(i = start; i < stop; ++i)
            {
                ofile << in[i] << endl;
            }
        }
    }

    template<class T> void NormGlobalVector(
            const int n,
            Array<OneD, const T> &in,
            std::ostream &out,
            MultiRegions::AssemblyMapCGSharedPtr &map)
    {
        Array<OneD, NekDouble> vExchange(1);
        Array<OneD, int> m_map = map->GetGlobalToUniversalMapUnique();
        vExchange[0] = Vmath::Dot2(n, in, in, m_map);
        map->GetComm()->AllReduce(vExchange, Nektar::LibUtilities::ReduceSum);
        out << "Norm: " << vExchange[0] << endl;
    }
}
#endif
