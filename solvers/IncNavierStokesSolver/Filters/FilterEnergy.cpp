///////////////////////////////////////////////////////////////////////////////
//
// File FilterEnergy.cpp
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
// Description: Output kinetic energy and enstrophy.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include <IncNavierStokesSolver/Filters/FilterEnergy.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string FilterEnergy::className = GetFilterFactory().
            RegisterCreatorFunction("Energy", FilterEnergy::create);

        FilterEnergy::FilterEnergy(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession), m_planes(), m_index(0), m_homogeneous(false)
        {
            std::string outName;
            if (pParams.find("OutputFile") == pParams.end())
            {
                outName = m_session->GetSessionName();
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                         "Missing parameter 'OutputFile'.");
                outName = pParams.find("OutputFile")->second;
            }

            m_comm = pSession->GetComm();
            outName += ".eny";
            if (m_comm->GetRank() == 0)
            {
                m_outFile.open(outName.c_str());
                ASSERTL0(m_outFile.good(), "Unable to open: '" + outName + "'");
            }
            
            ASSERTL0(pParams.find("OutputFrequency") != pParams.end(),
                     "Missing parameter 'OutputFrequency'.");
            m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());

            if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
                if (HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
                    HomoStr == "1D"            || HomoStr == "Homo1D")
                {
                    m_homogeneous = true;
                    pSession->LoadParameter("LZ", m_homogeneousLength);
                }
                else
                {
                    ASSERTL0(false, "The energy filter only supports 1D "
                                    "homogeneous expansions.");
                }
            }
        }

        FilterEnergy::~FilterEnergy()
        {

        }

        void FilterEnergy::v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            m_index = 0;
            MultiRegions::ExpListSharedPtr areaField;

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

            // Output values at initial time.
            v_Update(pFields, time);
        }

        void FilterEnergy::v_Update(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            int i, n, nPoints = pFields[0]->GetNpoints(), nPlanePts = 0;

            if (m_homogeneous)
            {
                nPlanePts = pFields[0]->GetPlane(0)->GetNpoints();
            }
            
            m_index++;
            if (m_index % m_outputFrequency > 0)
            {
                return;
            }

            // Calculate kinetic energy.
            NekDouble Ek = 0.0;
            Array<OneD, NekDouble> tmp(nPoints);
            Array<OneD, Array<OneD, NekDouble> > u(3);
            for (i = 0; i < 3; ++i)
            {
                if (m_homogeneous)
                {
                    pFields[i]->HomogeneousBwdTrans(pFields[i]->GetPhys(), tmp);
                }
                else
                {
                    tmp = pFields[i]->GetPhys();
                }

                u[i] = Array<OneD, NekDouble>(nPoints);
                Vmath::Vcopy(nPoints, tmp, 1, u[i], 1);
                Vmath::Vmul (nPoints, tmp, 1, tmp, 1, tmp, 1);

                if (m_homogeneous)
                {
                    pFields[i]->HomogeneousFwdTrans(tmp, tmp);
                    Ek += pFields[i]->GetPlane(0)->Integral(tmp) * 2.0 * M_PI;
                }
                else
                {
                    Ek += pFields[i]->Integral(tmp);
                }
            }

            Ek /= 2.0*m_area;

            if (m_comm->GetRank() == 0)
            {
                m_outFile << setw(10) << time << setw(20) << Ek;
            }

            bool waveSpace[3] = {
                pFields[0]->GetWaveSpace(),
                pFields[1]->GetWaveSpace(),
                pFields[2]->GetWaveSpace()
            };

            if (m_homogeneous)
            {
                for (i = 0; i < 3; ++i)
                {
                    pFields[i]->SetWaveSpace(false);
                }
            }

            // First calculate vorticity field.
            Array<OneD, NekDouble> tmp2(nPoints), tmp3(nPoints);
            Vmath::Zero(nPoints, tmp, 1);
            for (i = 0; i < 3; ++i)
            {
                int f1 = (i+2) % 3, c2 = f1;
                int c1 = (i+1) % 3, f2 = c1;
                pFields[f1]->PhysDeriv(c1, u[f1], tmp2);
                pFields[f2]->PhysDeriv(c2, u[f2], tmp3);
                Vmath::Vsub (nPoints, tmp2, 1, tmp3, 1, tmp2, 1);
                Vmath::Vvtvp(nPoints, tmp2, 1, tmp2, 1, tmp, 1, tmp, 1);
            }

            if (m_homogeneous)
            {
                for (i = 0; i < 3; ++i)
                {
                    pFields[i]->SetWaveSpace(waveSpace[i]);
                }
                pFields[0]->HomogeneousFwdTrans(tmp, tmp);
                Ek = pFields[0]->GetPlane(0)->Integral(tmp) * 2 * M_PI;
            }
            else
            {
                Ek = pFields[0]->Integral(tmp);
            }

            Ek /= 2.0*m_area;

            if (m_comm->GetRank() == 0)
            {
                m_outFile << setw(20) << Ek << endl;
            }
        }

        void FilterEnergy::v_Finalise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time)
        {
            m_outFile.close();
        }

        bool FilterEnergy::v_IsTimeDependent()
        {
            return true;
        }
    }
}
