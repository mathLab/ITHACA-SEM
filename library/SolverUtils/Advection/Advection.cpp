///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.cpp
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
namespace SolverUtils
{

/**
 * @returns The advection factory.
 */
AdvectionFactory& GetAdvectionFactory()
{
    typedef Loki::SingletonHolder<AdvectionFactory,
    Loki::CreateUsingNew,
    Loki::NoDestroy > Type;
    return Type::Instance();
}


/**
 * @param   pSession            Session configuration data.
 * @param   pFields             Array of ExpList objects.
 */
void Advection::InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    v_InitObject(pSession, pFields);
}


/**
 * @param   nConvectiveFields   Number of velocity components.
 * @param   pFields             Expansion lists for scalar fields.
 * @param   pAdvVel             Advection velocity.
 * @param   pInarray            Scalar data to advect.
 * @param   pOutarray           Advected scalar data.
 * @param   pTime               Simulation time.
 */
void Advection::Advect(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble> >        &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble> >        &pInarray,
    Array<OneD, Array<OneD, NekDouble> >              &pOutarray,
    const NekDouble                                   &pTime)
{
    v_Advect(nConvectiveFields, pFields, pAdvVel, pInarray, pOutarray, pTime);
}


/**
 * This function should be overridden in derived classes to initialise the
 * specific advection data members. However, this base class function should
 * be called as the first statement of the overridden function to ensure the
 * base class is correctly initialised in order.
 *
 * @param   pSession            Session information.
 * @param   pFields             Expansion lists for scalar fields.
 */
void Advection::v_InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    m_spaceDim = pFields[0]->GetCoordim(0);

    if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = pSession->GetSolverInfo("HOMOGENEOUS");
        if (HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
            HomoStr == "1D"            || HomoStr == "Homo1D")
        {
            m_spaceDim++;
        }
        else
        {
            ASSERTL0(false, "Only 1D homogeneous dimension supported.");
        }
    }
}


/**
 *
 */
void Advection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >    &inarray)
{
    ASSERTL0(false,
            "A baseflow is not appropriate for this advection type.");
}

}
}
