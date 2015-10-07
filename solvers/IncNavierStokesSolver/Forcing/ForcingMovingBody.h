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

        void CheckIsFromFile(const TiXmlElement* pForce);

        void InitialiseCableModel(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void InitialiseFilter(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement* pForce);

        void UpdateMotion(
            const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
            const Array<OneD, Array<OneD, NekDouble> >       &  fields,
                  NekDouble time );

        void TensionedCableModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                  Array<OneD, NekDouble> &FcePhysinArray,
                  Array<OneD, NekDouble> &MotPhysinArray);

        void EvaluateStructDynModel(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                  NekDouble time );

        void CalculateForcing(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields);

        void MappingBndConditions(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pfields,
            const Array<OneD, Array<OneD, NekDouble> >        & fields,
                  NekDouble time );

        void EvaluateAccelaration(
            const Array<OneD, NekDouble> &input,
                  Array<OneD, NekDouble> &output,
                  int npoints);

        void SetDynEqCoeffMatrix(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields);

        void RollOver(Array<OneD, Array<OneD, NekDouble> > &input);

        int m_movingBodyCalls;     ///< number of times the movbody have been called
        int m_np;                  ///< number of planes per processors
        int m_vdim;                ///< vibration dimension

        NekDouble m_structrho;     ///< mass of the cable per unit length
        NekDouble m_structdamp;    ///< damping ratio of the cable
        NekDouble m_lhom;          ///< length ratio of the cable
        NekDouble m_kinvis;        ///< fluid viscous
        NekDouble m_timestep;      ///< time step
        ///
        LibUtilities::NektarFFTSharedPtr m_FFT;
        ///
        FilterMovingBodySharedPtr m_MovBodyfilter;
        /// storage for the cable's force(x,y) variables
        Array<OneD, NekDouble> m_Aeroforces;
        /// storage for the cable's motion(x,y) variables
        Array<OneD, NekDouble> m_MotionVars;
        /// srorage for the velocity in z-direction
        Array<OneD, Array<OneD, NekDouble> > m_W;
        /// fictitious velocity storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fV;
        /// fictitious acceleration storage
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_fA;
        /// matrices in Newmart-beta method
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_A;
        /// matrices in Newmart-beta method
        Array<OneD, DNekMatSharedPtr> m_CoeffMat_B;
        /// [0] is displacements, [1] is velocities, [2] is accelerations
        Array<OneD, std::string> m_funcName;
        /// motion direction: [0] is 'x' and [1] is 'y'
        Array<OneD, std::string> m_motion;
        /// do determine if the the body motion come from an extern file
        Array<OneD, bool>        m_IsFromFile;
        /// Store the derivatives of motion variables in x-direction
        Array<OneD, Array< OneD, NekDouble> > m_zta;
        /// Store the derivatives of motion variables in y-direction
        Array<OneD, Array< OneD, NekDouble> > m_eta;
        ///
        Array<OneD, Array< OneD, NekDouble> > m_forcing;
};

}

#endif
