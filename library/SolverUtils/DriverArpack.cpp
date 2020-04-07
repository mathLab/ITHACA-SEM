///////////////////////////////////////////////////////////////////////////////
//
// File DriverArpack.cpp
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
// Description:  Arnoldi solver using Arpack
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/DriverArpack.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

std::string DriverArpack::arpackProblemTypeLookupIds[6] = {
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestReal"    ,0),
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestReal"   ,1),
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestImag"    ,2),
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestImag"   ,3),
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","LargestMag"     ,4),
    LibUtilities::SessionReader::RegisterEnumValue("ArpackProblemType","SmallestMag"    ,5),
};
std::string DriverArpack::arpackProblemTypeDefault = LibUtilities::SessionReader::RegisterDefaultSolverInfo("ArpackProblemType","LargestMag");
std::string DriverArpack::driverLookupId = LibUtilities::SessionReader::RegisterEnumValue("Driver","Arpack",0);

std::string DriverArpack::className = GetDriverFactory().RegisterCreatorFunction("Arpack", DriverArpack::create);

std::string DriverArpack::ArpackProblemTypeTrans[6] =
{ "LR", "SR", "LI", "SI", "LM", "SM" };


/**
 *
 */
DriverArpack::DriverArpack(const LibUtilities::SessionReaderSharedPtr pSession,
                           const SpatialDomains::MeshGraphSharedPtr pGraph)
    : DriverArnoldi(pSession, pGraph)
{
}


/**
 *
 */
DriverArpack::~DriverArpack()
{
}

/**
 *
 */
void DriverArpack::v_InitObject(ostream &out)
{
    DriverArnoldi::v_InitObject(out);

    //Initialisation of Arnoldi parameters
    m_maxn   = 1000000; // Maximum size of the problem
    m_maxnev = 200;      // maximum number of eigenvalues requested
    m_maxncv = 500;     // Largest number of basis vector used in Implicitly Restarted Arnoldi

    // Error alerts
    ASSERTL0(m_nvec <= m_maxnev,"NEV is greater than MAXNEV");
    ASSERTL0(m_kdim <= m_maxncv,"NEV is greater than MAXNEV");
    ASSERTL0(2      <= m_kdim-m_nvec,"NCV-NEV is less than 2");

    m_equ[0]->PrintSummary(out);

    ASSERTL0(m_comm->GetSize() == 1,
             "ARPACK Arnoldi solver does not support execution in parallel.");

    // Print session parameters
    out << "\tArnoldi solver type    : Arpack" << endl;

    out << "\tArpack problem type    : ";
    out << ArpackProblemTypeTrans[m_session->GetSolverInfoAsEnum<int>("ArpackProblemType")] << endl;
    DriverArnoldi::ArnoldiSummary(out);

    for( int i = 0; i < m_nequ; ++i)
    {
        m_equ[i]->DoInitialise();
    }

    // FwdTrans Initial conditions to be in Coefficient Space
    m_equ[m_nequ-1] ->TransPhysToCoeff();

}


void DriverArpack::v_Execute(ostream &out)

{
    Array<OneD, NekDouble> tmpworkd;

    int nq     = m_equ[0]->UpdateFields()[0]->GetNcoeffs(); // Number of points in the mesh
    int n      = m_nfields*nq;    // Number of points in eigenvalue calculation
    int lworkl = 3*m_kdim*(m_kdim+2); // Size of work array
    int ido ;     //REVERSE COMMUNICATION parameter. At the first call must be initialised at 0
    int info;     // do not set initial vector (info=0 random initial vector, info=1 read initial vector from session file)

    int iparam[11];
    int ipntr[14];

    Array<OneD, int> ritzSelect;
    Array<OneD, NekDouble> dr;
    Array<OneD, NekDouble> di;
    Array<OneD, NekDouble> workev;
    Array<OneD, NekDouble> z;
    NekDouble sigmar, sigmai;

    Array<OneD, NekDouble> resid(n);
    Array<OneD, NekDouble> v(n*m_kdim);
    Array<OneD, NekDouble> workl(lworkl, 0.0);
    Array<OneD, NekDouble> workd(3*n, 0.0);

    ASSERTL0(n <= m_maxn,  "N is greater than   MAXN");

    if(m_session->DefinesFunction("InitialConditions"))
    {
        out << "\tInital vector       : input file  " << endl;
        info = 1;
        CopyFieldToArnoldiArray(resid);
    }
    else
    {
        out << "\tInital vector       : random  " << endl;
        info = 0;
    }

    char B;

    iparam[0] = 1;      // strategy for shift-invert
    iparam[1] = 0;      // (deprecated)
    iparam[2] = m_nits; // maximum number of iterations allowed/taken
    iparam[3] = 1;      // blocksize to be used for recurrence
    iparam[4] = 0;      // number of converged ritz eigenvalues
    iparam[5] = 0;      // (deprecated)

    // Use generalized B matrix for coupled solver.
    if (m_timeSteppingAlgorithm)
    {
        iparam[6] = 1;      // computation mode 1=> matrix-vector prod
        B = 'I';
    }
    else {
        iparam[6] = 3;      // computation mode 1=> matrix-vector prod
        B = 'G';
    }
#if 0
    if((fabs(m_realShift) > NekConstants::kNekZeroTol)|| // use shift if m_realShift > 1e-12
       (fabs(m_imagShift) > NekConstants::kNekZeroTol))
    {
        iparam[6] = 3; // This was 3 need to know what to set it to
        B = 'G';
    }
    else
    {
        iparam[6] = 1;      // computation mode 1=> matrix-vector prod
        B = 'I';
    }
#endif
    iparam[7] = 0;      // (for shift-invert)
    iparam[8] = 0;      // number of MV operations
    iparam[9] = 0;      // number of BV operations
    iparam[10]= 0;      // number of reorthogonalisation steps

    int cycle = 0;
    const char* problem = ArpackProblemTypeTrans[m_session->GetSolverInfoAsEnum<int>("ArpackProblemType")].c_str();

    std::string name = m_session->GetSessionName() + ".evl";
    ofstream    pFile(name.c_str());

    ido     = 0;    //At the first call must be initialisedat 0

    while(ido != 99)//ido==-1 || ido==1 || ido==0)
    {
        //Routine for eigenvalue evaluation for non-symmetric operators
        Arpack::Dnaupd( ido, &B,       // B='I' for std eval problem
                        n, problem,  m_nvec,
                        m_evtol, &resid[0], m_kdim,
                        &v[0], n, iparam, ipntr, &workd[0],
                        &workl[0], lworkl, info);

        //Plotting of real and imaginary part of the
        //eigenvalues from workl
        out << "\rIteration " << cycle << ", output: " << info << ", ido=" << ido << " " << std::flush;

        if(!((cycle-1)%m_kdim)&&(cycle> m_kdim)&&(ido!=2))
        {
            pFile << "Krylov spectrum at iteration: " <<  cycle << endl;

            if(m_timeSteppingAlgorithm)
            {
                pFile << "EV  Magnitude   Angle       Growth      Frequency   Residual"    << endl;
            }
            else
            {
                pFile << "EV  Real        Imaginary   inverse real  inverse imag  Residual"   << endl;
            }

            out << endl;
            for(int k = 0; k < m_kdim; ++k)
            {
                // write m_kdim eigs to screen
                WriteEvs(pFile,k, workl[ipntr[5]-1+k],workl[ipntr[6]-1+k]);
            }
        }

        if (ido == 99) break;

        switch(ido)
        {
            case -1:
            case 1:  // Note that ido=1 we are using input x
                     // (workd[inptr[0]-1]) rather than Mx as
                     // recommended in manual since it is not
                     // possible to impose forcing directly.
                CopyArnoldiArrayToField(tmpworkd = workd + (ipntr[0]-1));

                m_equ[0]->TransCoeffToPhys();

                m_equ[0]->DoSolve();
                if(m_EvolutionOperator == eTransientGrowth)
                {
                    //start Adjoint with latest fields of direct
                    CopyFwdToAdj();

                    m_equ[1]->TransCoeffToPhys();
                    m_equ[1]->DoSolve();
                }

                if(!(cycle%m_infosteps))
                {
                    out << endl;
                    m_equ[0]->Output();
                }

                // operated fields are copied into workd[inptr[1]-1]
                CopyFieldToArnoldiArray(tmpworkd = workd + (ipntr[1]-1));

                cycle++;
                break;
            case 2: // provide y = M x (bwd trans and iproduct);
            {
                //workd[inptr[0]-1] copied into operator fields
                CopyArnoldiArrayToField(tmpworkd = workd + (ipntr[0]-1));

                m_equ[0]->TransCoeffToPhys();

                Array<OneD, MultiRegions::ExpListSharedPtr>  fields = m_equ[0]->UpdateFields();
                for (int i = 0; i < fields.size(); ++i)
                {
                    fields[i]->IProductWRTBase(fields[i]->GetPhys(),
                                               fields[i]->UpdateCoeffs());
                }

                // operated fields are copied into workd[inptr[1]-1]
                CopyFieldToArnoldiArray(tmpworkd = workd + (ipntr[1]-1));
                break;
            }
            default:
                ASSERTL0(false, "Unexpected reverse communication request.");
        }

    }

    out<< endl << "Converged in " << iparam[8] << " iterations" << endl;

    ASSERTL0(info >= 0," Error with Dnaupd");

    ritzSelect = Array<OneD, int>       (m_kdim,0);
    dr         = Array<OneD, NekDouble> (m_nvec+1,0.0);
    di         = Array<OneD, NekDouble> (m_nvec+1,0.0);
    workev     = Array<OneD, NekDouble> (3*m_kdim);
    z          = Array<OneD, NekDouble> (n*(m_nvec+1));

    if(m_negatedOp)
    {
        sigmar     = -m_realShift;
    }
    else
    {
        sigmar     = m_realShift;
    }

    // Do not pass imaginary shift to Arpack since we have not
    // used a Fortran complex number format and so processing
    // is mucked up. Need to do some processing afterwards.
    sigmai = 0;

    //Setting 'A', Ritz vectors are computed. 'S' for Shur vectors
    Arpack::Dneupd(1, "A", ritzSelect.get(), dr.get(), di.get(),
                   z.get(), n, sigmar, sigmai, workev.get(), &B,
                   n, problem, m_nvec, m_evtol, resid.get(), m_kdim,
                   v.get(), n, iparam, ipntr, workd.get(),
                   workl.get(),lworkl,info);

    ASSERTL0(info == 0, " Error with Dneupd");

    int nconv=iparam[4];

    // Subtract off complex shift if it exists
    if(m_negatedOp)
    {
        Vmath::Sadd(nconv,m_imagShift,di,1,di,1);
    }
    else
    {
        Vmath::Sadd(nconv,-m_imagShift,di,1,di,1);
    }

    WARNINGL0(m_imagShift == 0,"Complex Shift applied. "
              "Need to implement Ritz re-evaluation of"
              "eigenvalue. Only one half of complex "
              "value will be correct");


    Array<OneD, MultiRegions::ExpListSharedPtr>  fields = m_equ[0]->UpdateFields();

    out   << "Converged Eigenvalues: " << nconv << endl;
    pFile << "Converged Eigenvalues: " << nconv << endl;

    if(m_timeSteppingAlgorithm)
    {
        out   << "         Magnitude   Angle       Growth      Frequency" << endl;
        pFile << "         Magnitude   Angle       Growth      Frequency" << endl;
        for(int i= 0; i< nconv; ++i)
        {
            WriteEvs(out,i,dr[i],di[i]);
            WriteEvs(pFile,i,dr[i],di[i]);

            std::string file = m_session->GetSessionName() + "_eig_"
                + boost::lexical_cast<std::string>(i)
                + ".fld";
            WriteFld(file,z + i*n);
        }
    }
    else
    {
        out   << "         Real        Imaginary " << endl;
        pFile << "         Real        Imaginary " << endl;
        for(int i= 0; i< nconv; ++i)
        {
            WriteEvs(out,i,dr[i],di[i],
                     NekConstants::kNekUnsetDouble, false);
            WriteEvs(pFile,i,dr[i],di[i],
                     NekConstants::kNekUnsetDouble, false);

            std::string file = m_session->GetSessionName() + "_eig_"
                + boost::lexical_cast<std::string>(i)
                + ".fld";
            WriteFld(file,z + i*n);
        }
    }

    m_real_evl = dr;
    m_imag_evl = di;

    pFile.close();

    for(int j = 0; j < m_equ[0]->GetNvariables(); ++j)
    {
        NekDouble vL2Error = m_equ[0]->L2Error(j,false);
        NekDouble vLinfError = m_equ[0]->LinfError(j);
        if (m_comm->GetRank() == 0)
        {
            out << "L 2 error (variable " << m_equ[0]->GetVariable(j) << ") : " << vL2Error << endl;
            out << "L inf error (variable " << m_equ[0]->GetVariable(j) << ") : " << vLinfError << endl;
        }
    }
}

}
}
