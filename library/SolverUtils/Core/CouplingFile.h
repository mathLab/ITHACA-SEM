///////////////////////////////////////////////////////////////////////////////
//
// File CouplingFile.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: File based Coupling class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COUPLINGFILE
#define NEKTAR_COUPLINGFILE

#include <SolverUtils/Core/Coupling.h>
#include <SolverUtils/Core/SessionFunction.h>

namespace Nektar
{
namespace SolverUtils
{

class CouplingFile;

typedef boost::shared_ptr<CouplingFile> CouplingFileSharedPointer;

class CouplingFile : public Coupling
{

public:
    static std::string className;

    /// Creates an instance of this class
    static CouplingSharedPointer create(MultiRegions::ExpListSharedPtr field)
    {
        CouplingSharedPointer p =
            MemoryManager<CouplingFile>::AllocateSharedPtr(field);
        return p;
    }

    CouplingFile(MultiRegions::ExpListSharedPtr field);

    ~CouplingFile();

protected:
    virtual void v_Send(const int step,
                        const NekDouble time,
                        const Array<OneD, const Array<OneD, NekDouble> > &field,
                        LibUtilities::FieldMetaDataMap &fieldMetaDataMap);

    virtual void v_Receive(const int step,
                           const NekDouble time,
                           Array<OneD, Array<OneD, NekDouble> > &field,
                           LibUtilities::FieldMetaDataMap &fieldMetaDataMap);

    virtual void v_Finalize(void)
    {};

private:

    SessionFunctionSharedPtr m_inputFunction;
};
}
}

#endif
