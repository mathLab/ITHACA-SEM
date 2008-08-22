///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.h
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
// Description: Basic Advection Diffusion Reaction Field definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_AUXILIARY_ADVECTIONDIFFUSIONREACTION_H
#define NEKTAR_SOLVERS_AUXILIARY_ADVECTIONDIFFUSIONREACTION_H

#include <MultiRegions/DisContField2D.h>

namespace Nektar
{     
    /**
     * \brief This class is the base class for the development of solvers.
     *
     * It is basically a class handling vector valued fields where every field is
     * a DisContField2D class 
     */
    
    class AdvectionDiffusionReaction
    {
    public:           
      
        /**
         * Default constructor. 
         * 
         */ 
        AdvectionDiffusionReaction();


        /**
         * Constructor.
         * /param 
         * 
         *
         */
        AdvectionDiffusionReaction(SpatialDomains::MeshGraph2D &graph2D,
                                   SpatialDomains::BoundaryConditions &bcs,
                                   int variables = 1);
      
        void SetInitialConditions(SpatialDomains::BoundaryConditions &bcs, int initialtime = 0);
      
        void GetExactSolutions(SpatialDomains::BoundaryConditions &bcs);
      
        void AdvectionOperation(const Array<OneD,const NekDouble>& a, const Array<OneD,const NekDouble>& b);

        void DiffusionOperation(const Array<OneD,const NekDouble>& a, const Array<OneD,const NekDouble>& b);
      
        void FwdTrans(const AdvectionDiffusionReaction &In);
      
        void BwdTrans(void);
     
        void BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, int field_no);

        void BwdTrans(Array<OneD, Array<OneD, NekDouble> > &outarray);
    
        /**
         * 
         * \return ...
         */
        int GetNcoeffs(void)
        {
            return m_fields[0]->GetNcoeffs();
        }
      
        int GetNpoints(void)
        {
            return m_fields[0]->GetTrace()->GetNpoints();
        }

        int GetPointsTot(void)
        {
            return m_fields[0]->GetPointsTot();
        }
      
        inline int GetNvariables(void) const
        {
            return m_nvariables;
        }
      
        inline NekDouble GetTime(void)
        {
            return m_time;
        }
      
        inline void SetTime(NekDouble time)
        {
            m_time = time;
        }

        inline NekDouble GetTimeStep(void)
        {
            return m_dt;
        }

        inline void SetTimeStep(NekDouble dt)
        {
            m_dt = dt;
        }
      
        inline int GetSteps(void)
        {
            return m_steps;
        }

        inline void SetSteps(int steps)
        {
            m_steps = steps;
        }

        inline int GetCheckSteps(void)
        {
            return m_chksteps;
        }

        inline void SetCheckSteps(int check)
        {
            m_chksteps = check;
        }

        void GetCoords(Array<OneD, NekDouble>& x0, Array<OneD, NekDouble>& x1, Array<OneD, NekDouble>& x2=NullNekDouble1DArray);
      
        const Array<OneD, const NekDouble> &GetPhys(int field_no = -1);

        void GetPhys(Array<OneD, Array<OneD, NekDouble> >&outarray);


        //-----------------------------
        /**
         *
         */
        void SetPhys(Array<OneD, Array<OneD, NekDouble> >&inarray, int field_no = -1);
        //-----------------------------
      
        void SetPhys(Array<OneD, NekDouble> &inarray, int field_no = -1);
      
        void WriteToFile(std::ofstream &out, OutputFormat format, int field_no);
      
        NekDouble L2(const MultiRegions::ExpList2D &In, int field_no);
      
        void ExtractTracePhys(Array<OneD, NekDouble> &out, int field_no);

        void GetFwdBwdTracePhys(Array<OneD, NekDouble> &Fwd, Array<OneD,NekDouble> &Bwd, int field_no);
      
        void GetTraceNormals(Array<OneD, Array<OneD, NekDouble> > &Normals);

        void UpwindTrace(Array<OneD, Array<OneD, const NekDouble> > &Vec, 
                         Array<OneD, const NekDouble> &Fwd, 
                         Array<OneD, const NekDouble> &Bwd, 
                         Array<OneD, NekDouble> &Upwind); 
      
        /*  inline const boost::shared_ptr<Array<OneD, Array< OneD, NekDouble> > >GetTraceInnerNormals(void) const */
        /* 	{ */
        /* 	  return m_fields[0]->GetTraceInnerNormals();  */
        /* 	} */

     
        void IProductWRTDerivBase(const int dir, const Array< OneD, const NekDouble > &in,
                                  Array< OneD, NekDouble > &out, int field_no);
      
        void IProductWRTBase(const Array< OneD, const NekDouble > &in,
                             Array< OneD, NekDouble > &out, int field_no);

        void MultiplyByElmtInvMass(const Array< OneD, const NekDouble > &in,
                                   Array< OneD, NekDouble > &out, int field_no);

        void AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                              Array<OneD, const NekDouble> &Fy, 
                              Array<OneD, NekDouble> &outarray,
                              int field_no);
      
        Array<OneD, NekDouble> &UpdateCoeffs(int field_no);
	   
        const Array<OneD, const NekDouble> &GetCoeffs(int field_no);


      
        /*       Array<OneD, SpatialDomains::BoundaryConditionType>  &GetBndTypes(int field_no) */
        /*     { */
        /* 	return m_fields[field_no]->GetBndTypes(); */
        /*       } */
      
      
        /*       inline Array<OneD, SpatialDomains::Equation>  &GetBndUserTypes(int field_no) */
        /* 	{ */
        /* 	  return m_fields[field_no]->GetBndUserTypes(); */
        /* 	} */
      
        /*       inline Array<OneD, SpatialDomains::Equation>  &GetBndEquations(int field_no) */
        /* 	{ */
        /* 	  return m_fields[field_no]->GetBndEquations(); */
        /* 	} */
      
        /*       Array<OneD, MultiRegions::ExpList1DSharedPtr> &GetBndConstraint(int field_no) */
        /* 	{ */
        /* 	  return m_fields[field_no]->GetBndConstraint(); */
        /* 	} */
      
        enum UpwindType
        {           ///< flux not defined
            eNotSet,  ///< averaged (or centred) flux
            eAverage, ///< simple upwind flux
            eUpwind,  ///< local Lax-Friedrich flux
            eLLF,     ///< Harten-Lax-Leer Contact corrected flux
            eHLLC,    ///< Roe flux
            eRoe,    
        };

    protected:
        Array<OneD, MultiRegions::DisContField2DSharedPtr> m_fields; ///< Array holding all dependent variables
      
    private: 
      
        int m_nvariables;  ///< no of dependent variables
        NekDouble m_time;  ///< continous time
        NekDouble m_dt;    ///< time step size
        int m_steps;       ///< number of steps to be taken during the simulation
        int m_chksteps;    ///< number of steps between dumping check files
      
    };
    
    typedef boost::shared_ptr<AdvectionDiffusionReaction> AdvectionDiffusionReactionSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADVECTIONDIFFUSIONREACTION_H

/**
* $Log: $
**/
