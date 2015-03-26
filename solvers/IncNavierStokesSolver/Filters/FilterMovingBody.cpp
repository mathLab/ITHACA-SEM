///////////////////////////////////////////////////////////////////////////////
//
// File FilterMovingBody.cpp
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
// Description: Output moving body motion during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <iomanip>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <IncNavierStokesSolver/Filters/FilterMovingBody.h>

namespace Nektar
{

std::string FilterMovingBody::className = SolverUtils::GetFilterFactory().
        RegisterCreatorFunction("MovingBody",
                                FilterMovingBody::create,
                                "Moving Body Filter");
/**
 *
 */
FilterMovingBody::FilterMovingBody(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::map<std::string, std::string> &pParams)
    : Filter(pSession),
      m_session(pSession)
{
    // Determine output file name
    if (pParams.find("OutputFile") == pParams.end())
    {
        m_outputFile = pSession->GetSessionName();
    }
    else
    {
        ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                 "Missing parameter 'OutputFile'.");

        m_outputFile = pParams.find("OutputFile")->second;
    }

    if (!(m_outputFile.length() >= 4 &&
          m_outputFile.substr(m_outputFile.length() - 4) == ".mot"))
    {
        m_outputFile += ".mot";
    }
    // Load output frequency
    if (pParams.find("OutputFrequency") == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        m_outputFrequency = atoi(
                        pParams.find("OutputFrequency")->second.c_str());
    }
    //Check if we have a homogeneous case
    bool isHomogeneous1D;
    pSession->MatchSolverInfo("Homogeneous", "1D", isHomogeneous1D, false);

    ASSERTL0(isHomogeneous1D, "Moving Body implemented just for 3D "
                                "Homogeneous 1D discetisations.");
}


/**
 *
 */
FilterMovingBody::~FilterMovingBody()
{

}


/**
 *
 */
void FilterMovingBody::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index = 0;

    // Write header
    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();
    if (vComm->GetRank() == 0)
    {
        // Open output stream for cable motions
        m_outputStream.open(m_outputFile.c_str());
        m_outputStream << "#";
        m_outputStream.width(7);
        m_outputStream << "Time";
        m_outputStream.width(25);
        m_outputStream << "z";
        m_outputStream.width(25);
        m_outputStream << "Disp_x";
        m_outputStream.width(25);
        m_outputStream << "Vel_x";
        m_outputStream.width(25);
        m_outputStream << "Acel_x";
        m_outputStream.width(25);
        m_outputStream << "Disp_y";
        m_outputStream.width(25);
        m_outputStream << "Vel_y";
        m_outputStream.width(25);
        m_outputStream << "Acel_y";
        m_outputStream << endl;
    }
}

/**
 *
 */
void FilterMovingBody::UpdateMotion(
        const LibUtilities::SessionReaderSharedPtr              &pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble>                            &MotionVars,
        const NekDouble                                         &time)
{
    // Only output every m_outputFrequency.
    if ((m_index++) % m_outputFrequency)
    {
        return;
    }

    // Get the number of local planes on the process and their IDs
    // to properly locate the forces in the Fx, Fy etc. vectors.
    Array<OneD, unsigned int> ZIDs;
    ZIDs = pFields[0]->GetZIDs();
    int local_planes = ZIDs.num_elements();

    LibUtilities::CommSharedPtr vColComm
                            = pFields[0]->GetComm()->GetColumnComm();

    //
    if(!pSession->DefinesSolverInfo("HomoStrip"))
    {
        int Num_z_pos = pFields[0]->GetHomogeneousBasis()->GetNumModes();
        Array<OneD, NekDouble> z_coords(Num_z_pos,0.0);
        Array<OneD, const NekDouble> pts
                            = pFields[0]->GetHomogeneousBasis()->GetZ();

        NekDouble LZ;
        pSession->LoadParameter("LZ", LZ);
        Vmath::Smul(Num_z_pos,LZ/2.0,pts,1,z_coords,1);
        Vmath::Sadd(Num_z_pos,LZ/2.0,z_coords,1,z_coords,1);

        //get and output moving body variables
        int nStrVars = 3;
        Array<OneD, Array<OneD, NekDouble> > Motion_x(nStrVars);
        Array<OneD, Array<OneD, NekDouble> > Motion_y(nStrVars);

        for(int i = 0; i < nStrVars; i++)
        {
            Motion_x[i] = Array<OneD, NekDouble>(local_planes,0.0);
            Motion_y[i] = Array<OneD, NekDouble>(local_planes,0.0);
        }

        for(int plane = 0; plane < local_planes; plane++)
        {
            for (int var = 0; var < nStrVars; var++)
            {
                int xoffset = var*local_planes+plane;
                int yoffset = nStrVars*local_planes+xoffset;
                Motion_x[var][plane] = MotionVars[xoffset];
                Motion_y[var][plane] = MotionVars[yoffset];
            }
        }

        Array <OneD, NekDouble> CableAccelX;
        Array <OneD, NekDouble> CableVelocX;
        Array <OneD, NekDouble> CableDisplX;
        Array <OneD, NekDouble> CableAccelY;
        Array <OneD, NekDouble> CableVelocY;
        Array <OneD, NekDouble> CableDisplY;

        int npoints = Motion_x[0].num_elements();
        CableAccelX = Array <OneD, NekDouble>(npoints);
        CableVelocX = Array <OneD, NekDouble>(npoints);
        CableDisplX = Array <OneD, NekDouble>(npoints);
        CableAccelY = Array <OneD, NekDouble>(npoints);
        CableVelocY = Array <OneD, NekDouble>(npoints);
        CableDisplY = Array <OneD, NekDouble>(npoints);

        Vmath::Vcopy(npoints, Motion_x[0], 1, CableDisplX, 1);
        Vmath::Vcopy(npoints, Motion_x[1], 1, CableVelocX, 1);
        Vmath::Vcopy(npoints, Motion_x[2], 1, CableAccelX, 1);
        Vmath::Vcopy(npoints, Motion_y[0], 1, CableDisplY, 1);
        Vmath::Vcopy(npoints, Motion_y[1], 1, CableVelocY, 1);
        Vmath::Vcopy(npoints, Motion_y[2], 1, CableAccelY, 1);

        int colrank = vColComm->GetRank();
        int nproc   = vColComm->GetSize();
        // Send to root process.
        if (colrank == 0)
        {
            for (int j = 0; j <Motion_x[0].num_elements(); j++)
            {
                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                m_outputStream.width(25);
                m_outputStream << setprecision(6) << z_coords[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableDisplX[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableVelocX[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableAccelX[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableDisplY[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableVelocY[j];
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableAccelY[j];
                m_outputStream << endl;
            }

            for (int i = 1; i < nproc; ++i)
            {
                vColComm->Recv(i, CableAccelX);
                vColComm->Recv(i, CableVelocX);
                vColComm->Recv(i, CableDisplX);
                vColComm->Recv(i, CableAccelY);
                vColComm->Recv(i, CableVelocY);
                vColComm->Recv(i, CableDisplY);


                for (int j = 0; j < Motion_x[0].num_elements(); ++j)
                {
                    int n = Num_z_pos/nproc * i + j;
                    m_outputStream.width(8);
                    m_outputStream << setprecision(6) << time;
                    m_outputStream.width(25);
                    m_outputStream << setprecision(6) << z_coords[n];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableDisplX[j];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableVelocX[j];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableAccelX[j];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableDisplY[j];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableVelocY[j];
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableAccelY[j];
                    m_outputStream << endl;
                }
            }
        }
        else
        {
            vColComm->Send(0, CableAccelX);
            vColComm->Send(0, CableVelocX);
            vColComm->Send(0, CableDisplX);
            vColComm->Send(0, CableAccelY);
            vColComm->Send(0, CableVelocY);
            vColComm->Send(0, CableDisplY);
        }
    }
    else
    {
        int colrank = vColComm->GetRank();
        int nstrips;

        NekDouble DistStrip;

        pSession->LoadParameter("Strip_Z", nstrips);
        pSession->LoadParameter("DistStrip", DistStrip);

        Array<OneD, NekDouble> z_coords(nstrips);
        for(int i = 0; i < nstrips; i++)
        {
            z_coords[i] = i * DistStrip;
        }

        //get and output moving body variables
        int nStrVars = 3;
        Array<OneD, Array<OneD, NekDouble> > Motion_x(nStrVars);
        Array<OneD, Array<OneD, NekDouble> > Motion_y(nStrVars);

        for(int i = 0; i < nStrVars; i++)
        {
            Motion_x[i] = Array<OneD, NekDouble>(local_planes,0.0);
            Motion_y[i] = Array<OneD, NekDouble>(local_planes,0.0);
        }

        for(int plane = 0; plane < local_planes; plane++)
        {
            for (int var = 0; var < nStrVars; var++)
            {
                int xoffset = var*local_planes+plane;
                int yoffset = nStrVars*local_planes+xoffset;
                Motion_x[var][plane] = MotionVars[xoffset];
                Motion_y[var][plane] = MotionVars[yoffset];
            }
        }

        Array <OneD, NekDouble> CableMotions(6);

        for(int var = 0; var <nStrVars; var++)
        {
            CableMotions[var]   = Motion_x[var][0];
            CableMotions[3+var] = Motion_y[var][0];
        }
        // Send to root process.
        if (colrank == 0)
        {
            m_outputStream.width(8);
            m_outputStream << setprecision(6) << time;
            m_outputStream.width(25);
            m_outputStream << setprecision(6) << z_coords[0];
            for(int var = 0; var < 2*nStrVars; var++)
            {
                m_outputStream.width(25);
                m_outputStream << setprecision(8) << CableMotions[var];
            }
            m_outputStream << endl;

            for (int i = 1; i < nstrips; ++i)
            {
                vColComm->Recv(i, CableMotions);

                m_outputStream.width(8);
                m_outputStream << setprecision(6) << time;
                m_outputStream.width(25);
                m_outputStream << setprecision(6) << z_coords[i];
                for(int var = 0; var < 2*nStrVars; var++)
                {
                    m_outputStream.width(25);
                    m_outputStream << setprecision(8) << CableMotions[var];
                }
                m_outputStream << endl;
            }
        }
        else
        {
            for(int i = 1; i < nstrips; i++)
            {
                if(colrank == i)
                {
                    vColComm->Send(0, CableMotions);
                }
            }
        }
    }
}


/**
 *
 */
void FilterMovingBody::v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble                                         &time)
{
    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream.close();
    }
}


/**
 *
 */
bool FilterMovingBody::v_IsTimeDependent()
{
    return true;
}
}
