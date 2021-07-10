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
    int spacedim = 2;

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e1D,
             "1D expansion not supported for mean filter");

    ASSERTL0(pFields[0]->GetExpType() != MultiRegions::e3DH2D,
             "Homogeneous 2D expansion not supported for mean filter");
    
    // Lock equation system pointer
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer expired");

    auto fluidEqu = std::dynamic_pointer_cast<FluidInterface>(equ);
    ASSERTL0(fluidEqu, "Mean filter is incompatible with this solver.");
    
    if (pFields[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        m_homogeneous = true;
        spacedim = 3;
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
    std::string volname[3] = {"length", "area", "volume"};
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    if (vComm->GetRank() == 0)
    {
        m_outputStream.open(m_outputFile.c_str());
        ASSERTL0(m_outputStream.good(), "Unable to open: '" + m_outputFile + "'");
        m_outputStream.setf(ios::scientific, ios::floatfield);
        m_outputStream << "# Time";
        for(int i=0; i<pFields.size(); ++i)
            m_outputStream << setw(22) << equ->GetVariable(i);
        m_outputStream << setw(22) << volname[spacedim - 1] << " " << m_area;
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
    
    int i;
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    // Calculate average values.
    Array<OneD, NekDouble> avg(pFields.size());
    for (i = 0; i < pFields.size(); ++i)
    {
        avg[i] = 0.0;
    }
    
    if (m_homogeneous)
    {
        for (i = 0; i < pFields.size(); ++i)
            avg[i] = pFields[0]->GetPlane(0)->Integral(pFields[i]->GetPhys()) * m_homogeneousLength;
    }
    else
    {
        for (i = 0; i < pFields.size(); ++i)
            avg[i] = pFields[0]->Integral(pFields[i]->GetPhys());
    }
    
    for (i = 0; i < pFields.size(); ++i)
        avg[i] /= m_area;

    if (vComm->GetRank() == 0)
    {
        m_outputStream << setw(17) << setprecision(8) << time;
        for(int i=0; i<pFields.size(); ++i)
            m_outputStream << setw(22) << setprecision(11) << avg[i];
        m_outputStream << endl;
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
