////////////////////////////////////////////////////////////////////////////////
//
//  File: CompressData.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
#define NOMINMAX

#include <LibUtilities/BasicUtils/CompressData.h>
#include <LibUtilities/BasicConst/GitRevision.h>

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/lexical_cast.hpp>

#include <set>
#include <cstdint>

#ifdef NEKTAR_USE_MPI
#include <mpi.h>
#endif


#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace Nektar
{
namespace LibUtilities
{

    /**
     * run time determination of endianness, returning an EndianType
     */
    EndianType Endianness(void)
    {
        union
        {
            std::uint32_t value;
            std::uint8_t  data[sizeof(std::uint32_t)];
        } number;

        number.data[0] = 0x00;
        number.data[1] = 0x01;
        number.data[2] = 0x02;
        number.data[3] = 0x03;

        switch (number.value)
        {
        case UINT32_C(0x00010203): return eEndianBig;
        case UINT32_C(0x03020100): return eEndianLittle;
        case UINT32_C(0x02030001): return eEndianBigWord;
        case UINT32_C(0x01000302): return eEndianLittleWord;
        default:                   return eEndianUnknown;
        }
    }

    namespace CompressData
    {

        /**
         * Return a string describing this compression and endianness
         */
        std::string GetCompressString(void)
        {
            return  "B64Z-"+ EndianTypeMap[Endianness()];
        }

        std::string GetBitSizeStr(void)
        {
            return boost::lexical_cast<std::string>(sizeof(void*)*8);
        }

        /**
         * Convert a binary string to Base 64 string
         */
        void BinaryStrToBase64Str(std::string &compressedDataString,
                                  std::string &base64string)
        {
            // If the string length is not divisible by 3,
            // pad it. There is a bug in transform_width
            // that will make it reference past the end
            // and crash.
            switch (compressedDataString.length() % 3)
            {
            case 1:
                compressedDataString += '\0';
                /* Falls through. */
            case 2:
                compressedDataString += '\0';
                break;
            }

            // Convert from binary to base64.
            typedef boost::archive::iterators::base64_from_binary<
                boost::archive::iterators::transform_width<
                    std::string::const_iterator, 6, 8> > base64_t;
            base64string = std::string(base64_t(compressedDataString.begin()),
                                       base64_t(compressedDataString.end()));
        }

        /**
         * Convert a Base 64 string into a binary string
         */
        void Base64StrToBinaryStr(std::string &base64string,
                                  std::string &compressedDataString)
        {
            // Convert from base64 to binary.
            typedef boost::archive::iterators::transform_width<
                boost::archive::iterators::binary_from_base64<
                    std::string::const_iterator>, 8, 6 > binary_t;
            compressedDataString = std::string(binary_t(base64string.begin()),
                                               binary_t(base64string.end()));
        }
    }
}
}
