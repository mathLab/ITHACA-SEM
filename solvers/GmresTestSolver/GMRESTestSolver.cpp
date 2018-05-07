///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolver.cpp
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
// Description: Advection Diffusion Reaction framework solver
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/VmathArray.hpp> 
#include <LibUtilities/LinearAlgebra/NekLinSysIterative.h> 
#include <MultiRegions/AssemblyMap/AssemblyMap.h>


#define GMRES_DEBUG

using namespace std;
using namespace Nektar;
//using namespace MultiRegions;

#ifdef GMRES_DEBUG        
            
    // stores the A of the linear system Ax = f
    

#endif

    class GMRESTestLinSys
    {
        public:
        //GMRESTestLinSys();
        //~GMRESTestLinSys();
        void DoMatrixMultiply(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput);
        void preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);
        void initializeLinSys(const LibUtilities::SessionReaderSharedPtr &session);
        void LinSysSolve(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput);
        unsigned int Getnlinsys() const;
        protected:
        Array<OneD,     NekDouble>  m_mat; // stores the f of the linear system Ax = f
        //Array<OneD,     NekDouble>  m_rhs;
        NekLinSysIterativeSharedPtr m_NekLinSys;
        unsigned int                m_nlinsys;
        LinSysOperators             m_oprtors;
        
    };

    void GMRESTestLinSys::DoMatrixMultiply(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput)
    {
        for(int i=0;i<m_nlinsys;++i)
        {
            pOutput[i]  =   Vmath::Dot(m_nlinsys,&m_mat[i*m_nlinsys],&pInput[0]);
        }
        return;
    }
    void GMRESTestLinSys::preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out)
    {     
        int ntotal     = inarray.num_elements();
        Vmath::Vcopy(ntotal,inarray,1,out,1);
        return;
    }
    //void initializeLinSys()
    void GMRESTestLinSys::initializeLinSys(const LibUtilities::SessionReaderSharedPtr &session)
    {
        
        m_nlinsys                   = session->GetParameter("LinSysDimens");
        int nSubMatrix              = 2;
        
        m_mat       =  Array<OneD, NekDouble>(m_nlinsys*m_nlinsys);
        //Vmath::Fill(m_nlinsys,1.0,&m_rhs[0],1);
        Vmath::Zero(m_nlinsys*m_nlinsys,&m_mat[0],1);
        
        // initial the matrix(A) in Ax = f
        // for(int ij=0; ij<m_nlinsys*m_nlinsys; ++ij)
        // {
        //     m_mat[ij]    =   10.0+ij;
        // }
        // for(int i=0; i<m_nlinsys; ++i)
        // {
        //     m_mat[i*m_nlinsys+i]    =   m_mat[i*m_nlinsys+i]*10.0;
        // }
        
        int       nnn=0;
        NekDouble sum = 0.0;
        for(int i=0; i<m_nlinsys-nSubMatrix; ++i)
        {   
            sum = 0.0;
            for(int j=0; j<m_nlinsys-nSubMatrix; ++j)
            {
                m_mat[i*m_nlinsys+j] = 1.0+sqrt(1.0*nnn);
                nnn++;
                sum += m_mat[i*m_nlinsys+j];
            }
            m_mat[i*m_nlinsys+i]    =   -m_mat[i*m_nlinsys+i]+ sum +1.0;
        }
        for(int i=m_nlinsys-nSubMatrix; i<m_nlinsys; ++i)
        {
            m_mat[i*m_nlinsys+i]    =   1.0;
        }
        
        NekDouble tmp;
        cout <<"The Matrix A is:"<<endl;
        for(int i=0; i<m_nlinsys; ++i)
        {
            cout <<"i="<<i<<"   :";
            for(int j=0; j<m_nlinsys; ++j)
            {
                tmp =   m_mat[i*m_nlinsys+j];
                cout<< tmp<<"    ";
            }
            cout <<endl;
        }

        m_oprtors.DefineMatrixMultiply(&GMRESTestLinSys::DoMatrixMultiply, this);
        m_oprtors.DefinePrecond(&GMRESTestLinSys::preconditioner,this);

        LibUtilities::CommSharedPtr                 v_Comm;

        m_NekLinSys =   NekLinSysIterative::CreateInstance(session,v_Comm);
        m_NekLinSys->setLinSysOperators(m_oprtors);

        return;
    }

    //void initializeLinSys()
    void GMRESTestLinSys::LinSysSolve(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput)
    {

        m_NekLinSys->SolveLinearSystem(m_nlinsys,pInput,pOutput,0);
        
        return;
    }

    inline unsigned int GMRESTestLinSys::Getnlinsys() const
    {
        return m_nlinsys;
    }
        

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr pppsession;
    SpatialDomains::MeshGraphSharedPtr   graph;


    //Create session reader.
    pppsession = LibUtilities::SessionReader::CreateInstance(argc, argv);
    // Read the geometry and the expansion information
    graph = SpatialDomains::MeshGraph::Read(pppsession);
    // Get some information about the session
    string       sessionName    = pppsession->GetSessionName();
    string       outFile        = sessionName + ".fld";
    
    cout <<setprecision(6)<<setw(10)<<scientific;
    
    cout << 1.0<<endl;

    GMRESTestLinSys linsys;
    
    linsys.initializeLinSys(pppsession);
    //initializeLinSys();
    int ndim = linsys.Getnlinsys();

    Array<OneD, NekDouble>  rhs(ndim,1.0);
    Array<OneD, NekDouble> soltn_cg(ndim);
    //DoConjugateGradient(ndim, , soltn_cg, NoUsePtr, 0);
    Array<OneD, NekDouble> soltn_gmres(ndim);
    linsys.LinSysSolve(rhs, soltn_gmres);

#ifdef GMRES_DEBUG        
    cout<<"the gmres solution is :"<<endl;
    for(int jj=0;jj<ndim; ++jj)
    {
        cout << soltn_gmres[jj]<<endl;
    }
#endif   

    NekDouble l2error=0.0;
    NekDouble tmp;
    for(int i; i< ndim; ++i)
    {   
        tmp = soltn_cg[i]-soltn_gmres[i];
        l2error += tmp*tmp;
    }
    cout << "The L2 difference between CG&GMRES is"<<l2error<<endl;
    // Finalise session
    // pppsession->Finalise();

    return 0;
}
