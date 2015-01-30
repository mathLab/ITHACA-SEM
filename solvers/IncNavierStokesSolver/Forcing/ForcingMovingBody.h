///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.h
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
// Description: Moving Body (Wavyness and acceleration)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGBODY

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/FFT/NektarFFT.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <IncNavierStokesSolver/Filters/FilterMovingBody.h>

namespace Nektar
{

class ForcingMovingBody : public SolverUtils::Forcing
{
    public:

        friend class MemoryManager<ForcingMovingBody>;

        /// Creates an instance of this class
        static SolverUtils::ForcingSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const unsigned int& pNumForcingFields,
                const TiXmlElement* pForce)
        {
            SolverUtils::ForcingSharedPtr p =
                                    MemoryManager<ForcingMovingBody>::
                                            AllocateSharedPtr(pSession);
            p->InitObject(pFields, pNumForcingFields, pForce);
            return p;
        }

        ///Name of the class
        static std::string className;

    protected:
        virtual void v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int&                         pNumForcingFields,
            const TiXmlElement*                         pForce);

        virtual void v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
            const Array<OneD, Array<OneD, NekDouble> >& inarray,
                  Array<OneD, Array<OneD, NekDouble> >& outarray,
            const NekDouble&                            time);

    private:
        ForcingMovingBody(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        void CheckIsFromFile();

        void InitialiseCableModel(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void UpdateMotion(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                  NekDouble time);

        void TensionedCableModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                  Array<OneD, NekDouble> &FcePhysinArray,
                  Array<OneD, NekDouble> &MotPhysinArray);

        void EvaluateStructDynModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time );

        void CalculateForcing(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);

        void EvaluateAccelaration(
            const Array<OneD, NekDouble> &input,
                  Array<OneD, NekDouble> &output,
                  int npoints);

        void SetDynEqCoeffMatrix(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void OutputStructMotion(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble &time);

        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

        int m_intSteps;
        int m_movingBodyCalls;
        int m_NumLocPlane;
        int m_VarArraysize;
        int m_NumD;
        bool m_FictitiousMass;
        bool m_homostrip;

        NekDouble m_structrho;
        NekDouble m_fictrho;
        NekDouble m_cabletension;
        NekDouble m_bendingstiff;
        NekDouble m_structdamp;
        NekDouble m_fictdamp;
        NekDouble m_structstiff;
        NekDouble m_lhom;
        NekDouble m_kinvis;
        NekDouble m_timestep;

        static NekDouble StifflyStable_Betaq_Coeffs[3][3];
        static NekDouble StifflyStable_Alpha_Coeffs[3][3];
        static NekDouble StifflyStable_Gamma0_Coeffs[3];

        LibUtilities::NektarFFTSharedPtr m_FFT;
        LibUtilities::CommSharedPtr m_comm;
        FilterMovingBodySharedPtr m_filter;

        /// free vibration or forced vibration types are available
        std::string m_vibrationtype;

        /// storage for the cable's force(x,y) variables
        Array<OneD, NekDouble> m_Aeroforces;
        /// storage for the cable's motion(x,y) variables
        Array<OneD, NekDouble> m_MotionVars;
        Array<OneD, Array<OneD, NekDouble> > m_zeta;
        Array<OneD, Array<OneD, NekDouble> > m_eta;
        Array<OneD, Array<OneD, NekDouble> > m_W;
        /// fictitious velocity storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fV;
        /// fictitious acceleration storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fA;
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_BndV;
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_A;
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_B;
        /// [0] is displacements, [1] is velocities, [2] is accelerations
        Array<OneD, std::string> m_funcName;
        /// motion direction: [0] is 'x' and [1] is 'y'
        Array<OneD, std::string> m_motion;
        /// do determine if the the body motion come from an extern file
        Array<OneD, bool>        m_IsFromFile;

};

}

#endif
