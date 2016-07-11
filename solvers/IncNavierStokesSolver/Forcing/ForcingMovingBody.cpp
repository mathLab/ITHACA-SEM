///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.cpp
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
// Description: Moving Body m_forcing (movement of a body in a domain is achieved
// via a m_forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
std::string ForcingMovingBody::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("MovingBody",
                                    ForcingMovingBody::create,
                                    "Moving Body Forcing");

ForcingMovingBody::ForcingMovingBody(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : Forcing(pSession)
{
}

void ForcingMovingBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    // Just 3D homogenous 1D problems can use this techinque
    ASSERTL0(pFields[0]->GetExpType()==MultiRegions::e3DH1D,
             "Moving body implemented just for 3D Homogenous 1D expansions.");

    // At this point we know in the xml file where those quantities
    // are declared (equation or file) - via a function name which is now
    // stored in funcNameD etc. We need now to fill in with this info the
    // m_zta and m_eta vectors (actuallythey are matrices) Array to control if
    // the motion is determined by an equation or is from a file.(not Nektar++)
    // check if we need to load a file or we have an equation
    CheckIsFromFile(pForce);

    // Initialise movingbody filter
    InitialiseFilter(m_session, pFields, pForce);

    // Initialise the cable model
    InitialiseCableModel(m_session, pFields);

    // Load mapping
    m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);
    m_mapping->SetTimeDependent( true );

    if(m_vdim > 0)
    {
        m_mapping->SetFromFunction( false );
    }

    m_zta = Array<OneD, Array< OneD, NekDouble> > (3);
    m_eta = Array<OneD, Array< OneD, NekDouble> > (3);
    // What are this bi-dimensional vectors ------------------------------------------
    // m_zta[0] = m_zta                     |  m_eta[0] = m_eta                      |
    // m_zta[1] = d(m_zta)/dt               |  m_eta[1] = d(m_eta)/dt                |
    // m_zta[2] = dd(m_zta)/ddtt            |  m_eta[2] = dd(m_eta)/ddtt             |
    //--------------------------------------------------------------------------------
    int phystot = pFields[0]->GetTotPoints();
    for(int i = 0; i < m_zta.num_elements(); i++)
    {
        m_zta[i] = Array<OneD, NekDouble>(phystot,0.0);
        m_eta[i] = Array<OneD, NekDouble>(phystot,0.0);
    }
}

void ForcingMovingBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    // Update the forces from the calculation of fluid field, which is
    // implemented in the movingbody filter
    Array<OneD, NekDouble> Hydroforces (2*m_np,0.0);
    m_MovBodyfilter->UpdateForce(m_session, pFields, Hydroforces, time);

    // for "free" type (m_vdim = 2), the cable vibrates both in streamwise and crossflow
    // dimections, for "constrained" type (m_vdim = 1), the cable only vibrates in crossflow
    // dimection, and for "forced" type (m_vdim = 0), the calbe vibrates specifically along
    // a given function or file
    if(m_vdim == 1 || m_vdim == 2)
    {
        // For free vibration case, displacements, velocities and acceleartions
        // are obtained through solving structure dynamic model
         EvaluateStructDynModel(pFields, Hydroforces, time);

        // Convert result to format required by mapping
        int physTot = pFields[0]->GetTotPoints();
        Array< OneD, Array< OneD, NekDouble> >  coords(3);
        Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);
        for(int i =0; i<3; i++)
        {
            coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
            coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
        }
        // Get original coordinates
        pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

        // Add displacement to coordinates
        Vmath::Vadd(physTot, coords[0], 1, m_zta[0], 1, coords[0], 1);
        Vmath::Vadd(physTot, coords[1], 1, m_eta[0], 1, coords[1], 1);
        // Copy velocities
        Vmath::Vcopy(physTot, m_zta[1], 1, coordsVel[0], 1);
        Vmath::Vcopy(physTot, m_eta[1], 1, coordsVel[1], 1);
        
        // Update mapping
        m_mapping->UpdateMapping(time, coords, coordsVel);
    }
    else if(m_vdim == 0)
    {
        // For forced vibration case, load from specified file or function
        int cnt = 0;
        for(int j = 0; j < m_funcName.num_elements(); j++)
        {
            if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
            {
                ASSERTL0(false, "Motion loading from file needs specific "
                                "implementation: Work in Progress!");
            }
            else
            {
                EvaluateFunction(pFields, m_session, m_motion[0], m_zta[j],
                                 m_funcName[j], time);
                EvaluateFunction(pFields, m_session, m_motion[1], m_eta[j],
                                 m_funcName[j], time);
                cnt = cnt + 2;
            }
        }
        
        // Update mapping
        m_mapping->UpdateMapping(time);

        // Convert result from mapping       
        int physTot = pFields[0]->GetTotPoints();
        Array< OneD, Array< OneD, NekDouble> >  coords(3);
        Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);
        for(int i =0; i<3; i++)
        {
            coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
            coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
        }
        // Get original coordinates
        pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

        // Get Coordinates and coord velocity from mapping
        m_mapping->GetCartesianCoordinates(m_zta[0], m_eta[0], coords[2]);
        m_mapping->GetCoordVelocity(coordsVel);

        // Calculate displacement
        Vmath::Vsub(physTot, m_zta[0], 1, coords[0], 1, m_zta[0], 1);
        Vmath::Vsub(physTot, m_eta[0], 1, coords[1], 1, m_eta[0], 1);

        Vmath::Vcopy(physTot, coordsVel[0], 1, m_zta[1], 1);
        Vmath::Vcopy(physTot, coordsVel[1], 1, m_eta[1], 1);

        for(int var = 0; var < 3; var++)
        {
            for(int plane = 0; plane < m_np; plane++)
            {
                int n = pFields[0]->GetPlane(plane)->GetTotPoints();
                int offset  = plane * n;
                int Offset = var * m_np+plane;

                m_MotionVars[0][Offset] = m_zta[var][offset];
                m_MotionVars[1][Offset] = m_eta[var][offset];
            }
        }
    }
    else
    {
        ASSERTL0(false, 
                 "Unrecogized vibration type for cable's dynamic model");
    }

    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    // Pass the variables of the cable's motion to the movingbody filter
    if(colrank == 0)
    {
        int n = m_MotionVars[0].num_elements();
        Array<OneD, NekDouble> tmpArray(2*n),tmp(n);
        Vmath::Vcopy(n,m_MotionVars[0],1,tmpArray,1);
        Vmath::Vcopy(n,m_MotionVars[1],1,tmp=tmpArray+n,1);
        m_MovBodyfilter->UpdateMotion(m_session, pFields, tmpArray, time);
    }
}



/**
 *
 */
void ForcingMovingBody::EvaluateStructDynModel(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &Hydroforces,
              NekDouble  time)
{
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();

    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);

    //number of structral modes and number of strips
    int npts, nstrips;
    if(!homostrip) //full resolutions
    {
        npts = m_session->GetParameter("HomModesZ");
    }
    else
    {
        m_session->LoadParameter("HomStructModesZ", npts);
        m_session->LoadParameter("Strip_Z", nstrips);
         
        //ASSERTL0(nstrips < = npts,
        //     "the number of struct. modes must be larger than that of the strips.");
    }

    //the hydrodynamic forces
    Array<OneD, Array <OneD, NekDouble> > fces(2);
    //forces in x-direction
    fces[0] = Array <OneD, NekDouble> (npts,0.0);
    //forces in y-direction
    fces[1] = Array <OneD, NekDouble> (npts,0.0);

    //fill the force vectors
    if(colrank == 0)
    {
        if(!homostrip) //full resolutions
        {
            Vmath::Vcopy(m_np, Hydroforces,      1, fces[0], 1);
            Vmath::Vcopy(m_np, Hydroforces+m_np, 1, fces[1], 1);
        }
        else //strip modelling
        {
            fces[0][0] = Hydroforces[0];
            fces[1][0] = Hydroforces[m_np];
        }

        if(!homostrip) //full resolutions
        { 
            Array<OneD, NekDouble> tmp(2*m_np);
            for (int i = 1; i < nproc; ++i)
            {
                vcomm->GetColumnComm()->Recv(i, tmp);
                for(int n = 0; n < m_np; n++)
                {
                    for(int j = 0; j < 2; ++j)
                    {
                        fces[j][i*m_np + n] = tmp[j*m_np + n];
                    }
                }
            }
        }
        else //strip modelling
        //if the body is submerged partly, the fces are filled partly
        //by the flow induced forces
        {
            Array<OneD, NekDouble> tmp(2);
            for(int i = 1; i < nstrips; ++i)
            {
                vcomm->GetColumnComm()->Recv(i, tmp);

                for(int j = 0 ; j < 2; ++j)
                {
                    fces[j][i] = tmp[j];
                }
            }
        }
    }
    else
    {
        if(!homostrip) //full resolutions
        {
            vcomm->GetColumnComm()->Send(0, Hydroforces);
        }
        else //strip modelling
        {
            for(int i = 1; i < nstrips; ++i)
            {
                if(colrank == i)
                {
                    Array<OneD, NekDouble> tmp(2);
                    tmp[0] = Hydroforces[0];
                    tmp[1] = Hydroforces[m_np];
                    vcomm->GetColumnComm()->Send(0, tmp);
                }
            }
        }
    }

    if(colrank == 0)
    {
        // Fictitious mass method used to stablize the explicit coupling at
        // relatively lower mass ratio
        bool fictmass;
        m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                    fictmass, false);
        if(fictmass)
        {
            NekDouble fictrho, fictdamp;
            m_session->LoadParameter("FictMass", fictrho);
            m_session->LoadParameter("FictDamp", fictdamp);

            static NekDouble Betaq_Coeffs[2][2] = 
                                {{1.0,  0.0},{2.0, -1.0}};

            // only consider second order approximation for fictitious variables
            int  intOrder= 2;
            int  nint    = min(m_movingBodyCalls+1,intOrder);
            int  nlevels = m_fV[0].num_elements();

            for(int i = 0; i < m_motion.num_elements(); ++i)
            {
                RollOver(m_fV[i]);
                RollOver(m_fA[i]);

                int Voffset = npts;
                int Aoffset = 2*npts;

                Vmath::Vcopy(npts, m_MotionVars[i]+Voffset, 1, m_fV[i][0], 1);
                Vmath::Vcopy(npts, m_MotionVars[i]+Aoffset, 1, m_fA[i][0], 1);

                // Extrapolate to n+1
                Vmath::Smul(npts, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fV[i][nint-1],    1,
                            m_fV[i][nlevels-1], 1);
                Vmath::Smul(npts, 
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fA[i][nint-1],    1,
                            m_fA[i][nlevels-1], 1);

                for(int n = 0; n < nint-1; ++n)
                {
                    Vmath::Svtvp(npts, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fV[i][n],1,m_fV[i][nlevels-1],1,
                                 m_fV[i][nlevels-1],1);
                    Vmath::Svtvp(npts, 
                                 Betaq_Coeffs[nint-1][n],
                                 m_fA[i][n],1,m_fA[i][nlevels-1],1,
                                 m_fA[i][nlevels-1],1);
                }

                // Add the fictitious forces on the RHS of the equation
                Vmath::Svtvp(npts, fictdamp,m_fV[i][nlevels-1],1,
                             fces[i],1,fces[i],1);
                Vmath::Svtvp(npts, fictrho, m_fA[i][nlevels-1],1,
                             fces[i],1,fces[i],1);
            }
        }
    }
    //structural solver is implemented on the root process
    if(colrank == 0)
    {
        // Tensioned cable model is evaluated in wave space
        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {
            Newmark_betaSolver(pFields, fces[cn], m_MotionVars[cn]);
        }
    }

    Array<OneD, NekDouble> Motvars(2*2*m_np);
    // send physical coeffients to all planes of each processor
    if(!homostrip)//full resolutions
    {
        Array<OneD, NekDouble> tmp(2*2*m_np);

        if(colrank != 0)
        {
            vcomm->GetColumnComm()->Recv(0, tmp);
            Vmath::Vcopy(2*2*m_np, tmp, 1, Motvars, 1);
        }
        else
        {
            for (int i = 1; i < nproc; ++i)
            {
                for(int j = 0; j < 2; j++) //moving dimensions
                {
                    for(int k = 0; k < 2; k++) //disp. and vel.
                    {
                        for (int n = 0; n < m_np; n++)
                        {
                            tmp[j*2*m_np+k*m_np+n] = m_MotionVars[j][k*npts+i*m_np+n];
                        }
                    }
                }
                vcomm->GetColumnComm()->Send(i, tmp);
            }

            for(int j = 0; j < 2; j++)
            {
                for(int k = 0; k < 2; k++)
                {
                    for(int n = 0; n < m_np; n++)
                    {
                        tmp[j*2*m_np+k*m_np+n] = m_MotionVars[j][k*npts+n];
                    }
                    Vmath::Vcopy(2*2*m_np, tmp, 1, Motvars, 1);
                }
            }
        }
    }
    else //strip modelling
    {
        Array<OneD, NekDouble> tmp(2*2);

        if(colrank != 0)
        {
            for (int j = 1; j < nproc/nstrips; j++)
            {
                if(colrank == j*nstrips)
                {
                    vcomm->GetColumnComm()->Recv(0, tmp);

                    for(int plane = 0; plane < m_np; plane++)
                    {
                        for(int var = 0; var < 2; var++)
                        {
                            for(int k = 0; k < 2; k++)
                            {
                                Motvars[var*2*m_np+k*m_np+plane]= tmp[var*2+k];
                            }
                        }
                    }
                }
            }

            for(int i = 1; i < nstrips; i++)
            {
                for (int j = 0; j < nproc/nstrips; j++)
                {
                    if(colrank == i+j*nstrips)
                    {
                        vcomm->GetColumnComm()->Recv(0, tmp);

                        for(int plane = 0; plane < m_np; plane++)
                        {
                            for(int var = 0; var < 2; var++)
                            {
                                for(int k = 0; k < 2; k++)
                                {
                                    Motvars[var*2*m_np+k*m_np+plane] = tmp[var*2+k];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for(int j = 0; j < 2; ++j)
            {
                for(int k = 0; k < 2; ++k)
                {
                    tmp[j*2+k] = m_MotionVars[j][k*npts];
                }
            }

            for (int j = 1; j < nproc/nstrips; j++)
            {
                vcomm->GetColumnComm()->Send(j*nstrips, tmp);
            }

            for(int plane = 0; plane < m_np; plane++)
            {
                for(int j = 0; j < 2; j++)
                {
                    for(int k = 0; k < 2; ++k)
                    {
                        Motvars[j*2*m_np+k*m_np+plane] = m_MotionVars[j][k*npts];
                    }
                }
            }

            for(int i = 1; i < nstrips; ++i)
            {
                for(int j = 0; j < 2; ++j)
                {
                    for(int k = 0; k < 2; ++k)
                    {
                        tmp[j*2+k] = m_MotionVars[j][i+k*npts];
                    }
                }

                for (int j = 0; j < nproc/nstrips; j++)
                {
                    vcomm->GetColumnComm()->Send(i+j*nstrips, tmp);
                }
            }
        }
    }

    // Set the m_forcing term based on the motion of the cable
    for(int var = 0; var < 2; var++)
    {
        for(int plane = 0; plane < m_np; plane++)
        {
            int n = pFields[0]->GetPlane(plane)->GetTotPoints();

            Array<OneD, NekDouble> tmp;

            int offset  = plane * n;
            int xoffset = var * m_np+plane;
            int yoffset = 2*m_np + xoffset;

            Vmath::Fill(n, Motvars[xoffset], tmp = m_zta[var] + offset, 1);
            Vmath::Fill(n, Motvars[yoffset], tmp = m_eta[var] + offset, 1);
        }
    }
}

 
/**
 *
 */
void ForcingMovingBody::Newmark_betaSolver(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &HydroForces,
              Array<OneD, NekDouble> &BodyMotions)
{  
    std::string supptype = m_session->GetSolverInfo("SupportType");

    int npts = HydroForces.num_elements();
        
    Array<OneD, Array<OneD, NekDouble> > fft_i(4);
    Array<OneD, Array<OneD, NekDouble> > fft_o(4);

    for(int i = 0; i < 4; i++)
    {
        fft_i[i] = Array<OneD, NekDouble>(npts, 0.0);
        fft_o[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    Vmath::Vcopy(npts, HydroForces, 1, fft_i[0], 1);
    for(int i = 0; i < 3; i++)
    {
        Vmath::Vcopy(npts, BodyMotions+i*npts, 1, fft_i[i+1], 1);   
    }
    
    // Implement Fourier transformation of the motion variables
    if(boost::iequals(supptype, "Free-Free"))
    {
        for(int j = 0 ; j < 4; ++j)
        {
            m_FFT->FFTFwdTrans(fft_i[j], fft_o[j]);
        }
    }
    else if(boost::iequals(supptype, "Pinned-Pinned"))
    {
        //TODO:
        int N = fft_i[0].num_elements();

        for(int var = 0; var < 4; var++)
        {
            for(int i = 0; i < N; i++)
            {
                fft_o[var][i] = 0;

                for(int k = 0; k < N; k++)
                {
                    fft_o[var][i] +=
                        fft_i[var][k]*
                            sin(M_PI/(N)*(k+1/2)*(i+1));
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,
                    "Unrecognized support type for cable's motion");
    }

    // solve the ODE in the wave space
    for(int i = 0; i < npts; ++i)
    {
        int nrows = 3;

        Array<OneD, NekDouble> tmp0,tmp1,tmp2;
        tmp0 = Array<OneD, NekDouble> (3,0.0);
        tmp1 = Array<OneD, NekDouble> (3,0.0);
        tmp2 = Array<OneD, NekDouble> (3,0.0);

        for(int var = 0; var < 3; var++)
        {
            tmp0[var] = fft_o[var+1][i];
        }

        tmp2[0] = fft_o[0][i];

        Blas::Dgemv('N', nrows, nrows, 1.0,
                    &(m_CoeffMat_B[i]->GetPtr())[0],
                    nrows, &tmp0[0], 1,
                    0.0,   &tmp1[0], 1);
        Blas::Dgemv('N', nrows, nrows, 1.0/m_structrho,
                    &(m_CoeffMat_A[i]->GetPtr())[0],
                    nrows, &tmp2[0], 1,
                    1.0,   &tmp1[0], 1);

        for(int var = 0; var < 3; var++)
        {
            fft_i[var][i] = tmp1[var];
        }
    }

    // get physical coeffients via Backward fourier transformation of wave
    // coefficients
    if(boost::iequals(supptype, "Free-Free"))
    {
        for(int var = 0; var < 3; var++)
        {
            m_FFT->FFTBwdTrans(fft_i[var], fft_o[var]);
        }
    }
    else if(boost::iequals(supptype, "Pinned-Pinned"))
    {
        //TODO:
        int N = fft_i[0].num_elements();

        for(int var = 0; var < 3; var++)
        {
            for(int i = 0; i < N; i++)
            {
                fft_o[var][i] = 0;

                for(int k = 0; k < N; k++)
                {
                    fft_o[var][i] +=
                        fft_i[var][k]*
                            sin(M_PI/(N)*(k+1)*(i+1/2))*2/N;
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,
                    "Unrecognized support type for cable's motion");
    }
    

    for(int i = 0; i < 3; i++)
    {
        Array<OneD, NekDouble> tmp(npts,0.0);
        Vmath::Vcopy(npts, fft_o[i], 1, tmp = BodyMotions+i*npts, 1);
    }
}
  
/**
 *
 */
void ForcingMovingBody::InitialiseCableModel(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    m_movingBodyCalls = 0;
    m_session->LoadParameter("Kinvis",m_kinvis);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);

    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();

    //number of structral modes
    int npts;
    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);

   if(!homostrip) //full resolutions
    {
        npts = m_session->GetParameter("HomModesZ");
    }
    else
    {
        m_session->LoadParameter("HomStructModesZ", npts);
    }
  
    m_MotionVars = Array<OneD, Array<OneD, NekDouble> > (2);
    m_MotionVars[0] = Array<OneD, NekDouble>(3*npts,0.0);
    m_MotionVars[1] = Array<OneD, NekDouble>(3*npts,0.0);

    Array<OneD, unsigned int> ZIDs;
    ZIDs = pFields[0]->GetZIDs();
    m_np = ZIDs.num_elements();

    std::string vibtype = m_session->GetSolverInfo("VibrationType");

    if(boost::iequals(vibtype, "Constrained"))
    {
        m_vdim = 1;
    }
    else if (boost::iequals(vibtype, "Free"))
    {
        m_vdim = 2;
    }
    else if (boost::iequals(vibtype, "Forced"))
    {
        m_vdim = 0;
        return;
    }

    if(!homostrip)
    {
        m_session->LoadParameter("LZ", m_lhom);
        int nplanes = m_session->GetParameter("HomModesZ");
        m_FFT = 
            LibUtilities::GetNektarFFTFactory().CreateInstance(
                                            "NekFFTW", nplanes);
    }
    else
    {
        int nstrips;
        NekDouble DistStrip;

        m_session->LoadParameter("DistStrip", DistStrip);
        m_session->LoadParameter("Strip_Z", nstrips);
        m_lhom = nstrips * DistStrip;
        m_FFT = 
            LibUtilities::GetNektarFFTFactory().CreateInstance(
                                            "NekFFTW", nstrips);
    }

    // load the structural dynamic parameters from xml file
    m_session->LoadParameter("StructRho",  m_structrho);
    m_session->LoadParameter("StructDamp", m_structdamp, 0.0);

    // Identify whether the fictitious mass method is active for explicit
    // coupling of fluid solver and structural dynamics solver
    bool fictmass;
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                fictmass, false);
    if(fictmass)
    {
        NekDouble fictrho, fictdamp;
        m_session->LoadParameter("FictMass", fictrho);
        m_session->LoadParameter("FictDamp", fictdamp);
        m_structrho  += fictrho;
        m_structdamp += fictdamp;

        // Storage array of Struct Velocity and Acceleration used for
        // extrapolation of fictitious force
        m_fV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        m_fA = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        for (int i = 0; i < m_motion.num_elements(); ++i)
        {
            m_fV[i] = Array<OneD, Array<OneD, NekDouble> > (2);
            m_fA[i] = Array<OneD, Array<OneD, NekDouble> > (2);

            for(int n = 0; n < 2; ++n)
            {
                m_fV[i][n] = Array<OneD, NekDouble>(npts, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(npts, 0.0);
            }
        }
    }

    // Setting the coefficient matrices for solving structural dynamic ODEs
    SetDynEqCoeffMatrix(pFields);

    // Set initial condition for cable's motion
    int cnt = 0;
    
    for(int j = 0; j < m_funcName.num_elements(); j++)
    {
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
            std::ifstream inputStream;
            int nzpoints;

            if(homostrip)
            {
                m_session->LoadParameter("HomStructModesZ", nzpoints);
            }
            else
            { 
                nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
            }

            if (vcomm->GetRank() == 0)
            {
                std::string filename = m_session->GetFunctionFilename(
                    m_funcName[j], m_motion[0]);
                
                // Open intputstream for cable motions
                inputStream.open(filename.c_str());

                // Import the head string from the file
                Array<OneD, std::string> tmp(9);
                for(int n = 0; n < tmp.num_elements(); n++)
                {
                    inputStream >> tmp[n];
                }

                NekDouble time, z_cds;
                // Import the motion variables from the file
                for (int n = 0; n < nzpoints; n++)
                {
                    inputStream >> setprecision(6) >> time;
                    inputStream >> setprecision(6) >> z_cds;
                    inputStream >> setprecision(8) >> m_MotionVars[0][n];
                    inputStream >> setprecision(8) >> m_MotionVars[0][n+nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[0][n+2*nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n+nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n+2*nzpoints];
                }
                // Close inputstream for cable motions
                inputStream.close();
            }
            cnt = cnt + 2;
        }
        else //Evaluate from the functions specified in xml file
        {
            if(!homostrip)
            {
                int nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
                Array<OneD, NekDouble> z_coords(nzpoints,0.0);
                Array<OneD, const NekDouble> pts
                                    = pFields[0]->GetHomogeneousBasis()->GetZ();

                Vmath::Smul(nzpoints,m_lhom/2.0,pts,1,z_coords,1);
                Vmath::Sadd(nzpoints,m_lhom/2.0,z_coords,1,z_coords,1);

                Array<OneD, NekDouble> x0(m_np, 0.0);
                Array<OneD, NekDouble> x1(m_np, 0.0);
                Array<OneD, NekDouble> x2(m_np, 0.0);
                for (int plane = 0; plane < m_np; plane++)
                {
                    x2[plane] = z_coords[ZIDs[plane]];
                }

                Array<OneD, NekDouble> tmp (m_np,0.0);
                Array<OneD, NekDouble> tmp0(m_np,0.0);
                Array<OneD, NekDouble> tmp1(m_np,0.0);
                LibUtilities::EquationSharedPtr ffunc0,ffunc1;

                ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
                ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);
                ffunc0->Evaluate(x0, x1, x2, 0.0, tmp0);
                ffunc1->Evaluate(x0, x1, x2, 0.0, tmp1);

                int offset = j*npts;
                Vmath::Vcopy(m_np, tmp0, 1, tmp = m_MotionVars[0]+offset,1);
                Vmath::Vcopy(m_np, tmp1, 1, tmp = m_MotionVars[1]+offset,1);

                if(colrank == 0)
                {
                    for (int i = 1; i < nproc; ++i)
                    {
                        vcomm->GetColumnComm()->Recv(i, tmp0);
                        vcomm->GetColumnComm()->Recv(i, tmp1);
                        Vmath::Vcopy(m_np, tmp0, 1, tmp = m_MotionVars[0]+offset+i*m_np,1);
                        Vmath::Vcopy(m_np, tmp1, 1, tmp = m_MotionVars[1]+offset+i*m_np,1);
                    }
                }
                else
                {
                    vcomm->GetColumnComm()->Send(0, tmp0);
                    vcomm->GetColumnComm()->Send(0, tmp1);
                }
            }
            else
            {  
                if(colrank == 0)
                { 
                    int nstrips;
                    m_session->LoadParameter("Strip_Z", nstrips);

                    ASSERTL0(m_session->DefinesSolverInfo("USEFFT"),
                            "Fourier transformation of cable motion is currently "
                            "implemented only for FFTW module.");
                
                    NekDouble DistStrip;
                    m_session->LoadParameter("DistStrip", DistStrip);

                    Array<OneD, NekDouble> x0(npts, 0.0);
                    Array<OneD, NekDouble> x1(npts, 0.0);
                    Array<OneD, NekDouble> x2(npts, 0.0);
                    Array<OneD, NekDouble> tmp (npts,0.0);
                    Array<OneD, NekDouble> tmp0(npts,0.0); 
                    Array<OneD, NekDouble> tmp1(npts,0.0);
                    for (int plane = 0; plane < npts; plane++)
                    {
                        x2[plane] = plane*DistStrip;
                    }
                    LibUtilities::EquationSharedPtr ffunc0,ffunc1;
                    ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
                    ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);
                    ffunc0->Evaluate(x0, x1, x2, 0.0, tmp0);
                    ffunc1->Evaluate(x0, x1, x2, 0.0, tmp1);

                    int offset = j*npts;
                    Vmath::Vcopy(npts, tmp0, 1, tmp = m_MotionVars[0]+offset,1);
                    Vmath::Vcopy(npts, tmp1, 1, tmp = m_MotionVars[1]+offset,1);     
                }      
            }

           cnt = cnt + 2;
        }
    }
}


/**
 *
 */
void ForcingMovingBody::SetDynEqCoeffMatrix(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    int nplanes;

    bool homostrip;
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);

    if(!homostrip)
    {
        nplanes = m_session->GetParameter("HomModesZ"); 
    }
    else
    {
        m_session->LoadParameter("Strip_Z", nplanes);
    }

    m_CoeffMat_A = Array<OneD, DNekMatSharedPtr>(nplanes);
    m_CoeffMat_B = Array<OneD, DNekMatSharedPtr>(nplanes);

    NekDouble tmp1, tmp2, tmp3;
    NekDouble tmp4, tmp5, tmp6, tmp7;

    // load the structural dynamic parameters from xml file
    NekDouble cabletension;
    NekDouble bendingstiff;
    NekDouble structstiff;
    m_session->LoadParameter("StructStiff",  structstiff,  0.0);
    m_session->LoadParameter("CableTension", cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", bendingstiff, 0.0);

    tmp1 =   m_timestep * m_timestep;
    tmp2 =  structstiff / m_structrho;
    tmp3 = m_structdamp / m_structrho;
    tmp4 = cabletension / m_structrho;
    tmp5 = bendingstiff / m_structrho;

    // solve the ODE in the wave space for cable motion to obtain disp, vel and
    // accel

    std::string supptype = m_session->GetSolverInfo("SupportType");

    for(int plane = 0; plane < nplanes; plane++)
    {
        int nel = 3;
        m_CoeffMat_A[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
        m_CoeffMat_B[plane]
                = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

        // Initialised to avoid compiler warnings.
        unsigned int K = 0;
        NekDouble beta = 0.0;

        if (boost::iequals(supptype, "Free-Free"))
        {
            K = plane/2;
            beta = 2.0 * M_PI/m_lhom;
        }
        else if(boost::iequals(supptype, "Pinned-Pinned"))
        {
            K = plane+1;
            beta = M_PI/m_lhom;
        }
        else
        {
            ASSERTL0(false,
                        "Unrecognized support type for cable's motion");
        }

        tmp6 = beta * K;
        tmp6 = tmp6 * tmp6;
        tmp7 = tmp6 * tmp6;
        tmp7 = tmp2 + tmp4 * tmp6 + tmp5 * tmp7;

        (*m_CoeffMat_A[plane])(0,0) = tmp7;
        (*m_CoeffMat_A[plane])(0,1) = tmp3;
        (*m_CoeffMat_A[plane])(0,2) = 1.0;
        (*m_CoeffMat_A[plane])(1,0) = 0.0;
        (*m_CoeffMat_A[plane])(1,1) = 1.0;
        (*m_CoeffMat_A[plane])(1,2) =-m_timestep/2.0;
        (*m_CoeffMat_A[plane])(2,0) = 1.0;
        (*m_CoeffMat_A[plane])(2,1) = 0.0;
        (*m_CoeffMat_A[plane])(2,2) =-tmp1/4.0;
        (*m_CoeffMat_B[plane])(0,0) = 0.0;
        (*m_CoeffMat_B[plane])(0,1) = 0.0;
        (*m_CoeffMat_B[plane])(0,2) = 0.0;
        (*m_CoeffMat_B[plane])(1,0) = 0.0;
        (*m_CoeffMat_B[plane])(1,1) = 1.0;
        (*m_CoeffMat_B[plane])(1,2) = m_timestep/2.0;
        (*m_CoeffMat_B[plane])(2,0) = 1.0;
        (*m_CoeffMat_B[plane])(2,1) = m_timestep;
        (*m_CoeffMat_B[plane])(2,2) = tmp1/4.0;

        m_CoeffMat_A[plane]->Invert();
        (*m_CoeffMat_B[plane]) =
            (*m_CoeffMat_A[plane]) * (*m_CoeffMat_B[plane]);
    }
}


/**
 * Function to roll time-level storages to the next step layout.
 * The stored data associated with the oldest time-level
 * (not required anymore) are moved to the top, where they will
 * be overwritten as the solution process progresses.
 */
void ForcingMovingBody::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
{
    int nlevels = input.num_elements();

    Array<OneD, NekDouble> tmp;
    tmp = input[nlevels-1];

    for(int n = nlevels-1; n > 0; --n)
    {
        input[n] = input[n-1];
    }

    input[0] = tmp;
}

/**
 *
 */
void ForcingMovingBody::CheckIsFromFile(const TiXmlElement* pForce)
{

    m_funcName = Array<OneD, std::string> (3);
    m_motion = Array<OneD, std::string> (2);
    m_motion[0] = "x";
    m_motion[1] = "y";

    m_IsFromFile = Array<OneD, bool> (6);
    // Loading the x-dispalcement (m_zta) and the y-displacement (m_eta)
    // Those two variables are bith functions of z and t and the may come
    // from an equation (forced vibration) or from another solver which, given
    // the aerodynamic forces at the previous step, calculates the 
    // displacements.

    //Get the body displacement: m_zta and m_eta
    const TiXmlElement* funcNameElmt_D
                    = pForce->FirstChildElement("DISPLACEMENTS");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as d(z,t).");

    m_funcName[0] = funcNameElmt_D->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
             "Function '" + m_funcName[0] + "' not defined.");

    //Get the body velocity of movement: d(m_zta)/dt and d(m_eta)/dt
    const TiXmlElement* funcNameElmt_V
                    = pForce->FirstChildElement("VELOCITIES");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body velocity of movement as v(z,t).");

    m_funcName[1] = funcNameElmt_V->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
             "Function '" + m_funcName[1] + "' not defined.");


    //Get the body acceleration: dd(m_zta)/ddt and dd(m_eta)/ddt
    const TiXmlElement* funcNameElmt_A
                    = pForce->FirstChildElement("ACCELERATIONS");
    ASSERTL0(funcNameElmt_A,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body acceleration as a(z,t).");

    m_funcName[2] = funcNameElmt_A->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
             "Function '" + m_funcName[2] + "' not defined.");

    LibUtilities::FunctionType vType;

    // Check Displacement x
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[0]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[0] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[0] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Displacement y
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[1]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[1] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[1] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in y must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Velocity x
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[2] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[2] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in x must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Velocity y
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[3] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[3] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in y must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Acceleration x
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[4] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[4] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Acceleration y
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[5] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[5] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in y must be specified via an "
                        "equation <E> or a file <F>");
    }
}

/**
 *
 */
void ForcingMovingBody::InitialiseFilter(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement* pForce)
{
    // Get the outputfile name, output frequency and 
    // the boundary's ID for the cable's wall
    std::string typeStr = pForce->Attribute("TYPE");
    std::map<std::string, std::string> vParams;

    const TiXmlElement *param = pForce->FirstChildElement("PARAM");
    while (param)
    {
        ASSERTL0(param->Attribute("NAME"),
                 "Missing attribute 'NAME' for parameter in filter "
                 + typeStr + "'.");
        std::string nameStr = param->Attribute("NAME");

        ASSERTL0(param->GetText(), "Empty value string for param.");
        std::string valueStr = param->GetText();

        vParams[nameStr] = valueStr;

        param = param->NextSiblingElement("PARAM");
    }

    // Creat the filter for MovingBody, where we performed the calculation of
    // fluid forces and write both the aerodynamic forces and motion variables
    // into the output files
    m_MovBodyfilter = MemoryManager<FilterMovingBody>::
                                    AllocateSharedPtr(pSession, vParams);

    // Initialise the object of MovingBody filter
    m_MovBodyfilter->Initialise(pFields, 0.0);

}

}
