///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.h
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
// Description: CFL tester solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_CFLTESTER_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_CFLTESTER_H

#include <Auxiliary/UnsteadySystem.h>

namespace Nektar
{
	
	static NekDouble EigenvaluesRegMeshes[10][14] =
	{{18.849560, 37.699110, 62.254810, 92.114370, 127.199300, 167.209300, 212.189900, 261.852900, 316.311600, 375.280200, 438.911800, 506.922500, 579.489600, 656.332100},
	{32.277580, 67.919290, 113.860500, 169.431300, 234.151800, 307.678500, 389.746800, 480.143200, 578.689500, 685.233500, 799.642600, 921.799700, 1051.600000, 1188.949000},
	{56.418920, 115.463200, 191.215900, 282.444600, 388.469600, 508.597500, 642.523200, 789.745400, 950.116100, 1123.225000, 1309.009000, 1507.109000, 1717.513000, 1939.900000},
	{77.988010, 158.974000, 262.115100, 385.943600, 529.415800, 691.770100, 872.414200, 1070.867000, 1286.725000, 1519.644000, 1769.323000, 2035.496000, 2317.925000, 2616.397000},
	{102.518800, 207.512500, 340.960200, 500.891100, 686.042600, 895.336300, 1128.107000, 1383.623000, 1661.474000, 1961.100000, 2282.224000, 2624.386000, 2987.395000, 3370.856000},
	{125.653700, 253.549800, 415.661800, 609.711500, 834.071300, 1087.549000, 1369.216000, 1678.315000, 2014.214000, 2376.369000, 2764.308000, 3177.613000, 3615.910000, 4078.863000},
	{150.392200, 302.534600, 495.170300, 725.542200, 991.784100, 1292.405000, 1626.378000, 1992.722000, 2390.769000, 2819.785000, 3279.297000, 3768.717000, 4287.691000, 4835.729000},
	{174.106400, 349.636100, 571.545300, 836.751500, 1143.031000, 1488.755000, 1872.650000, 2293.681000, 2750.981000, 3243.808000, 3771.512000, 4333.522000, 4929.325000, 5558.461000},
	{198.947800, 398.852300, 651.400800, 953.050400, 1301.327000, 1694.312000, 2130.624000, 2609.014000, 3128.566000, 3688.361000, 4287.733000, 4925.959000, 5602.527000, 6316.841000},
	{222.978700, 446.558300, 728.737000, 1065.639000, 1454.438000, 1893.064000, 2379.896000, 2913.620000, 3493.132000, 4117.485000, 4785.857000, 5497.519000, 6251.820000, 7079.855200}};
	
	
    class CFLtester : public UnsteadySystem
    {
    public:
        friend class MemoryManager<CFLtester>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession) {
            EquationSystemSharedPtr p = MemoryManager<CFLtester>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        virtual ~CFLtester();
		
    protected:
        
		Array<OneD, Array<OneD, NekDouble> > m_velocity;

        CFLtester(const LibUtilities::SessionReaderSharedPtr& pSession);

        void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                      Array<OneD,  Array<OneD, NekDouble> > &outarray,
                      const NekDouble time);

        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                          Array<OneD,  Array<OneD, NekDouble> > &outarray,
                          const NekDouble time);

        virtual void v_InitObject();

        // DG Advection routines
        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

        // Print Summary
        virtual void v_PrintSummary(std::ostream &out);
    private:

		virtual NekDouble v_GetTimeStep(const Array<OneD,int> ExpOrder, 
										const Array<OneD,NekDouble> CFL, NekDouble timeCFL);
		
		virtual NekDouble v_GetTimeStep(int ExpOrder, NekDouble CFL, NekDouble TimeStability);

    };
}

#endif
