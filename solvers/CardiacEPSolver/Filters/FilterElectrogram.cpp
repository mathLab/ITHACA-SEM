///////////////////////////////////////////////////////////////////////////////
//
// File FilterElectrogram.cpp
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
// Description: Outputs virtual electrograms at specific locations.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <CardiacEPSolver/Filters/FilterElectrogram.h>

namespace Nektar
{
    std::string FilterElectrogram::className =
            GetFilterFactory().RegisterCreatorFunction(
                    "Electrogram",
                    FilterElectrogram::create);

    /**
     *
     */
    FilterElectrogram::FilterElectrogram(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::map<std::string, std::string> &pParams) :
        Filter(pSession)
    {
        if (pParams.find("OutputFile") == pParams.end())
        {
            m_outputFile = m_session->GetSessionName();
        }
        else
        {
            ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                     "Missing parameter 'OutputFile'.");
            m_outputFile = pParams.find("OutputFile")->second;
        }
        if (!(m_outputFile.length() >= 4
              && m_outputFile.substr(m_outputFile.length() - 4) == ".ecg"))
        {
            m_outputFile += ".ecg";
        }

        if (pParams.find("OutputFrequency") == pParams.end())
        {
            m_outputFrequency = 1;
        }
        else
        {
            m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());
        }

        ASSERTL0(pParams.find("Points") != pParams.end(),
                 "Missing parameter 'Points'.");
        m_electrogramStream.str(pParams.find("Points")->second);
        m_index = 0;
    }


    /**
     *
     */
    FilterElectrogram::~FilterElectrogram()
    {

    }


    /**
     *
     */
    void FilterElectrogram::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        ASSERTL0(!m_electrogramStream.fail(),
                 "No history points in stream.");

        m_index = 0;
        Array<OneD, NekDouble>  gloCoord(3,0.0);
        LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

        // Read electrogram points
        // Always use dim = 3 to allow electrode to be above surface
        const int dim = 3;
        const NekDouble tol = 1e-06;
        int i = 0;

        while (!m_electrogramStream.fail())
        {
            m_electrogramStream >> gloCoord[0] >> gloCoord[1] >> gloCoord[2];
            if (!m_electrogramStream.fail())
            {
                ASSERTL0(pFields[0]->GetExpIndex(gloCoord, tol) == -1,
                         "Electrogram point must lie outside of domain: ("
                        + boost::lexical_cast<string>(gloCoord[0]) + ","
                        + boost::lexical_cast<string>(gloCoord[1]) + ","
                        + boost::lexical_cast<string>(gloCoord[2]) + ")");

                SpatialDomains::VertexComponentSharedPtr vert
                    = MemoryManager<SpatialDomains::VertexComponent>
                    ::AllocateSharedPtr(dim, i, gloCoord[0],
                                        gloCoord[1], gloCoord[2]);

                m_electrogramPoints.push_back(vert);
                ++i;
            }
        }

        if (vComm->GetRank() == 0)
        {
            // Open output stream
            m_outputStream.open(m_outputFile.c_str());
            m_outputStream << "# Electrogram data for variables (:";

            for (i = 0; i < pFields.num_elements(); ++i)
            {
                m_outputStream << m_session->GetVariable(i) <<",";
            }

            m_outputStream << ") at points:" << endl;

            for (i = 0; i < m_electrogramPoints.size(); ++i)
            {
                m_electrogramPoints[i]->GetCoords(  gloCoord[0],
                                                gloCoord[1],
                                                gloCoord[2]);

                m_outputStream << "# \t" << i;
                m_outputStream.width(8);
                m_outputStream << gloCoord[0];
                m_outputStream.width(8);
                m_outputStream << gloCoord[1];
                m_outputStream.width(8);
                m_outputStream << gloCoord[2];
                m_outputStream << endl;
            }
        }

        // Compute the distance function for each electrogram point
        const unsigned int nq = pFields[0]->GetNpoints();
        NekDouble px, py, pz;
        m_grad_R_x = Array<OneD, Array<OneD, NekDouble> >(m_electrogramPoints.size());
        m_grad_R_y = Array<OneD, Array<OneD, NekDouble> >(m_electrogramPoints.size());
        m_grad_R_z = Array<OneD, Array<OneD, NekDouble> >(m_electrogramPoints.size());

        Array<OneD, NekDouble> x(nq);
        Array<OneD, NekDouble> y(nq);
        Array<OneD, NekDouble> z(nq);

        Array<OneD, NekDouble> oneOverR(nq);
        for (unsigned int i = 0; i < m_electrogramPoints.size(); ++i)
        {
            m_grad_R_x[i] = Array<OneD, NekDouble>(nq);
            m_grad_R_y[i] = Array<OneD, NekDouble>(nq);
            m_grad_R_z[i] = Array<OneD, NekDouble>(nq);

            // Compute 1/R
            m_electrogramPoints[i]->GetCoords(px,py,pz);

            pFields[0]->GetCoords(x,y,z);

            Vmath::Sadd   (nq, -px, x, 1, x, 1);
            Vmath::Sadd   (nq, -py, y, 1, y, 1);
            Vmath::Sadd   (nq, -pz, z, 1, z, 1);
            Vmath::Vvtvvtp(nq, x, 1, x, 1, y, 1, y, 1, oneOverR, 1);
            Vmath::Vvtvp  (nq, z, 1, z, 1, oneOverR, 1, oneOverR, 1);
            Vmath::Vsqrt  (nq, oneOverR, 1, oneOverR, 1);
            Vmath::Sdiv   (nq, 1.0, oneOverR, 1, oneOverR, 1);

            // Compute grad 1/R
            pFields[0]->PhysDeriv(oneOverR, m_grad_R_x[i], m_grad_R_y[i], m_grad_R_z[i]);
        }

        // Compute electrogram point for initial condition
        v_Update(pFields, time);
    }


    /**
     *
     */
    void FilterElectrogram::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        // Only output every m_outputFrequency.
        if ((m_index++) % m_outputFrequency)
        {
            return;
        }

        const unsigned int nq = pFields[0]->GetNpoints();
        const unsigned int npoints = m_electrogramPoints.size();
        LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

        unsigned int i = 0;
        Array<OneD, NekDouble> e(npoints);

        // Compute grad V
        Array<OneD, NekDouble> grad_V_x(nq), grad_V_y(nq), grad_V_z(nq);
        pFields[0]->PhysDeriv(pFields[0]->GetPhys(), grad_V_x, grad_V_y, grad_V_z);

        for (i = 0; i < npoints; ++i)
        {
            // Multiply together
            Array<OneD, NekDouble> output(nq);
            Vmath::Vvtvvtp(nq, m_grad_R_x[i], 1, grad_V_x, 1, m_grad_R_y[i], 1, grad_V_y, 1, output, 1);
            Vmath::Vvtvp  (nq, m_grad_R_z[i], 1, grad_V_z, 1, output, 1, output, 1);

            e[i] = pFields[0]->Integral(output);
        }

        // Exchange history data
        // This could be improved to reduce communication but works for now
        vComm->AllReduce(e, LibUtilities::ReduceSum);

        // Only the root process writes out electrogram data
        if (vComm->GetRank() == 0)
        {
            m_outputStream.width(8);
            m_outputStream << setprecision(6) << time;

            // Write data values point by point
            for (i = 0; i < m_electrogramPoints.size(); ++i)
            {
                m_outputStream.width(25);
                m_outputStream << setprecision(16) << e[i];
            }
            m_outputStream << endl;
        }
    }


    /**
     *
     */
    void FilterElectrogram::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        if (pFields[0]->GetComm()->GetRank() == 0)
        {
            m_outputStream.close();
        }
    }


    /**
     *
     */
    bool FilterElectrogram::v_IsTimeDependent()
    {
        return true;
    }
}
