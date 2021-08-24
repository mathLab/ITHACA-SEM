////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFileBase.h
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
//  Description: Base class for outputting to a file
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTFILEBASE
#define FIELDUTILS_OUTPUTFILEBASE

#include "../Module.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{

/// Converter from fld to vtk.
class OutputFileBase : public OutputModule
{
public:
    OutputFileBase(FieldSharedPtr f);
    virtual ~OutputFileBase();

    /// Write fld to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "OutputFileBase";
    }

    virtual std::string GetModuleDescription()
    {
        return "Writing file";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eOutput;
    }

protected:
    /// Write from pts to output file.
    virtual void OutputFromPts(po::variables_map &vm) = 0;

    /// Write from m_exp to output file.
    virtual void OutputFromExp(po::variables_map &vm) = 0;

    /// Write from data to output file.
    virtual void OutputFromData(po::variables_map &vm) = 0;

    virtual fs::path GetPath(std::string &filename,
                                    po::variables_map &vm) = 0;

    virtual fs::path GetFullOutName(std::string &filename,
                                    po::variables_map &vm) = 0;

    bool m_requireEquiSpaced;

private:
    bool WriteFile(std::string &filename, po::variables_map &vm);

    void ConvertExpToEquispaced(po::variables_map &vm);

    void PrintErrorFromPts();

    void PrintErrorFromExp();
};
}
}

#endif
