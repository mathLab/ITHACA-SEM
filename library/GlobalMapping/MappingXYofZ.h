///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXYofZ.h
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
// Description: Mapping of the type X = x + f(z), Y = y + g(z)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_GLOBALMAPPING_MAPPINGXYOFZ
#define NEKTAR_GLOBALMAPPING_MAPPINGXYOFZ

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <GlobalMapping/GlobalMappingDeclspec.h>
#include <GlobalMapping/Mapping.h>

namespace Nektar
{
namespace GlobalMapping
{
class MappingXYofZ: public Mapping
{
public:

    friend class MemoryManager<MappingXYofZ> ;

    /// Creates an instance of this class
    GLOBAL_MAPPING_EXPORT
    static MappingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr        &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement                                *pMapping)
    {
        MappingSharedPtr p =
                MemoryManager<MappingXYofZ>::AllocateSharedPtr(pSession,
                                                               pFields);
        p->InitObject(pFields, pMapping);
        return p;
    }

    ///Name of the class
    static std::string className;

protected:

    // Constructor
    MappingXYofZ(const LibUtilities::SessionReaderSharedPtr       &pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

    // Virtual functions
    GLOBAL_MAPPING_EXPORT
    virtual void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement                                *pMapping);

    GLOBAL_MAPPING_EXPORT virtual void v_ContravarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_CovarToCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray);            

    GLOBAL_MAPPING_EXPORT virtual void v_ContravarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_CovarFromCartesian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray); 

    GLOBAL_MAPPING_EXPORT virtual void v_GetJacobian(
        Array<OneD, NekDouble>               &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_DotGradJacobian(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, NekDouble>                            &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_GetMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_GetInvMetricTensor(
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_ApplyChristoffelContravar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_ApplyChristoffelCovar(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray);

    GLOBAL_MAPPING_EXPORT virtual void v_UpdateGeomInfo();

private:

};

}
}

#endif
