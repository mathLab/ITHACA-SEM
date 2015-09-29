////////////////////////////////////////////////////////////////////////////////
//
//  File: CompressData.h
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
//  Description: Routines for compressing and inflating data
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_COMPRESSDATA_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_COMPRESSDATA_H

#include <LibUtilities/BasicUtils/SessionReader.h>

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/assign/list_of.hpp>

namespace Nektar
{
namespace LibUtilities
{

    namespace CompressData
    {
        /**
         * Compress a vector of NekDouble values into a string using zlib.
         */
        int Deflate(std::vector<NekDouble>& in, std::string& out);


        /**
         * Convert a string containing compressed binary (i.e. from
         * deflate) into a base 64 string
         */
        void BinaryStrToBase64Str(std::string &compressedDataString,
                                 std::string &base64string);

        /**
         * Compress a vector of NekDouble values into a base64 string.
         */
        int DeflateToBase64Str(std::vector<NekDouble>& in, std::string& out64);

        /**
         * Decompress a zlib-compressed string into a vector of NekDouble
         * values.
         */
        int Inflate(std::string& in, std::vector<NekDouble>& out);


        /**
         * Convert a string containing base 64 (i.e. from xml file)
         * into a binary string
         */
        void Base64StrToBinaryStr(std::string &base64string,
                                  std::string &compressedDataString);



        /**
         * Decompress a base 64 compressed binary string into a vector
         * of NekDouble values.
         */
        int InflateFromBase64Str(std::string& in64,
                                 std::vector<NekDouble>& out);

    }
}
}
#endif
