////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.h
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
//  Description: Vtk output module
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTVTK
#define FIELDUTILS_OUTPUTVTK

#include "OutputFileBase.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{

/// Converter from fld to vtk.
class OutputVtk : public OutputFileBase
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<OutputVtk>::AllocateSharedPtr(f);
    }
    static ModuleKey m_className;

    OutputVtk(FieldSharedPtr f);
    virtual ~OutputVtk();

    virtual std::string GetModuleName()
    {
        return "OutputVtk";
    }

protected:
    /// Write from pts to output file.
    virtual void OutputFromPts(po::variables_map &vm);

    /// Write from m_exp to output file.
    virtual void OutputFromExp(po::variables_map &vm);

    /// Write from data to output file.
    virtual void OutputFromData(po::variables_map &vm);

    virtual fs::path GetPath(std::string &filename,
                                    po::variables_map &vm);

    virtual fs::path GetFullOutName(std::string &filename,
                                    po::variables_map &vm);

private:
    void WriteVtkHeader(std::ostream &outfile);

    void WriteVtkFooter(std::ostream &outfile);

    void WriteEmptyVtkPiece(std::ofstream &outfile);

    void WritePVtu(po::variables_map &vm);

    std::string PrepareOutput(po::variables_map &vm);
};
}
}

#endif
