////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldConvertComm.hpp
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
//  Description: Specialised FieldConvert communicator.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_FIELDCONVERTCOMM
#define FIELDUTILS_FIELDCONVERTCOMM

#include <LibUtilities/Communication/CommSerial.h>

namespace Nektar
{
namespace FieldUtils
{

class FieldConvertComm : public LibUtilities::CommSerial
{
public:
    FieldConvertComm(int argc, char *argv[], int size, int rank)
        : CommSerial(argc, argv)
    {
        m_size = size;
        m_rank = rank;
        m_type = "FieldConvert parallel";
    }
    FieldConvertComm(int size, int rank) : CommSerial(0, NULL)
    {
        m_size = size;
        m_rank = rank;
        m_type = "FieldConvert parallel";
    }
    virtual ~FieldConvertComm()
    {
    }
    void v_SplitComm(int pRows, int pColumns)
    {
        // Compute row and column in grid.
        m_commRow = std::shared_ptr<FieldConvertComm>(
            new FieldConvertComm(pColumns, m_rank));
        m_commColumn =
            std::shared_ptr<FieldConvertComm>(new FieldConvertComm(pRows, 0));
    }

protected:
    int v_GetRank(void)
    {
        return m_rank;
    }

    bool v_TreatAsRankZero(void)
    {
        return true;
    }

    bool v_IsSerial(void)
    {
        return true;
    }

    bool v_RemoveExistingFiles(void)
    {
        return false;
    }

private:
    int m_rank;
};
}
}
#endif
