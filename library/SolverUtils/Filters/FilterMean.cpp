///////////////////////////////////////////////////////////////////////////////
//
// File FilterMean.cpp
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
// Description: Output mean.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Filters/FilterMean.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{
std::string FilterMean::className = SolverUtils::GetFilterFactory().
    RegisterCreatorFunction("Mean", FilterMean::create);

FilterMean::FilterMean(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem>        &pEquation,
    const ParamMap &pParams)
    : Filter        (pSession, pEquation),
      m_index       (-1),
      m_homogeneous (false),
      m_planes      ()
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile = m_session->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");
        m_outputFile = it->second;
    }
    m_outputFile += ".avg";

    

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    ASSERTL0(it != pParams.end(), "Missing parameter 'OutputFrequency'.");
    LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
    m_outputFrequency = round(equ.Evaluate());
    
    pSession->LoadParameter("LZ", m_homogeneousLength, 0.0);
}

FilterMean::~FilterMean()
{

}

void FilterMean::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    MultiRegions::ExpListSharedPtr areaField;

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e1D,
             "1D expansion not supported for energy filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e2D,
             "2D expansion not supported for energy filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e3DH2D,
             "Homogeneous 2D expansion not supported for energy filter");

    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_homogeneous = true;
    }

    // Calculate area/volume of domain.
    if (m_homogeneous)
    {
        m_planes  = pFields[0]->GetZIDs();
        areaField = pFields[0]->GetPlane(0);
    }
    else
    {
        areaField = pFields[0];
    }

    Array<OneD, NekDouble> inarray(areaField->GetNpoints(), 1.0);
    m_area = areaField->Integral(inarray);

    if (m_homogeneous)
    {
        m_area *= m_homogeneousLength;
    }

    // Open OutputFile
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    if (vComm->GetRank() == 0)
    {
        m_outputStream.open(m_outputFile.c_str());
        ASSERTL0(m_outputStream.good(), "Unable to open: '" + m_outputFile + "'");
        m_outputStream.setf(ios::scientific, ios::floatfield);
        m_outputStream << "# Time                Ua"
                  << setw(22) << " Va"
                  << setw(22) << " Wa"
                  << setw(22) << " Area: " << m_area;
        m_outputStream << endl;
    }

    // Output values at initial time.
    m_index = 0;
    v_Update(pFields, time);
}

void FilterMean::v_Update(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    // Only output every m_outputFrequency
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }
    
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    int i, nPoints = pFields[0]->GetNpoints();
    // Lock equation system pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Mean filter is incompatible with this solver.");

    // Store physical values in an array
    Array<OneD, Array<OneD, NekDouble> > physfields(pFields.size());
    for(i = 0; i < pFields.size(); ++i)
    {
        physfields[i] = pFields[i]->GetPhys();
    }

    // Calculate average values.
    NekDouble ua, va, wa;
    ua = va = wa = 0.0;
    Array<OneD, Array<OneD, NekDouble> > u(3);
    for (i = 0; i < 3; ++i)
    {
        u[i] = Array<OneD, NekDouble>(nPoints);
    }
    fluidEqu->GetVelocity(physfields, u);

    if (m_homogeneous)
    {
        ua = pFields[0]->GetPlane(0)->Integral(u[0]) * m_homogeneousLength;
        va = pFields[0]->GetPlane(0)->Integral(u[1]) * m_homogeneousLength;
        wa = pFields[0]->GetPlane(0)->Integral(u[2]) * m_homogeneousLength;
    }
    else
    {
        ua = pFields[0]->Integral(u[0]);
        va = pFields[0]->Integral(u[1]);
        wa = pFields[0]->Integral(u[2]);
    }
    
    ua /= m_area;
    va /= m_area;
    wa /= m_area;

    if (vComm->GetRank() == 0)
    {
        m_outputStream << setw(17) << setprecision(8) << time
                  << setw(22) << setprecision(11) << ua
                  << setw(22) << setprecision(11) << va
                  << setw(22) << setprecision(11) << wa << endl;
    }
}

void FilterMean::v_Finalise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    boost::ignore_unused(pFields, time);

    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}

bool FilterMean::v_IsTimeDependent()
{
    return true;
}

}
}
