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

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

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
	
	static NekDouble EigenvaluesAnaMeshesAB2[6][14] =
	{{18.849560, 37.699110, 62.254810, 92.114370, 127.199300, 167.209300, 212.189900, 261.852900, 316.311600, 375.280200, 438.911800, 506.922500, 579.489600, 656.332100}, 
	{625.358147, 1196.328877, 1817.986152, 2475.195638, 3162.726621, 3867.689152, 4595.324620, 5340.967597, 6188.286000, 7399.972000, 8708.072000, 10111.160000, 11608.000000, 13197.470000}, 
	{705.738221, 1268.022370, 1883.167893, 2533.621042, 3215.460527, 3915.574276, 4638.993612, 5380.985954, 6153.131100, 6914.337599, 7721.969961, 8529.801758, 9326.138300, 10138.712574}, 
	{733.261678, 1293.053007, 1905.511959, 2553.681869, 3233.562783, 3932.056639, 4654.037626, 5394.854534, 6165.912400, 6926.174929, 7732.912435, 8509.531270, 9335.504451, 10147.409120}, 
	{747.053764, 1305.602041, 1916.689465, 2563.704820, 3242.557000, 3940.217950, 4661.410438, 5401.601883, 6172.029656, 6931.769221, 7737.958233, 8514.117977, 9339.576627, 10151.066350},
	{755.354971, 1313.162315, 1923.399136, 2569.711788, 3247.904628, 3945.042610, 4665.705043, 5405.487354, 6175.462753, 6934.840253, 7740.605266, 8518.215978, 9343.872193, 10156.126788}};
	
	static NekDouble EigenvaluesAnaMeshesRK2[6][14] =
	{{18.849560, 37.699110, 62.254810, 92.114370, 127.199300, 167.209300, 212.189900, 261.852900, 316.311600, 375.280200, 438.911800, 506.922500, 579.489600, 656.332100}, 
	{557.050244, 1050.871071, 1626.248117, 2260.761693, 2948.760275, 3677.949231, 4461.704082, 5269.668301, 6188.286000, 7399.972000, 8708.072000, 10111.160000, 11608.000000, 13197.470000}, 
	{628.920938, 1115.037602, 1684.555318, 2314.125522, 2997.926601, 3723.485221, 4504.103288, 5309.152432, 6144.425035, 6997.224333, 7899.393530, 8799.536577, 9741.728144, 10730.151339}, 
	{653.448558, 1137.048335, 1704.542817, 2332.448417, 3014.804194, 3739.158997, 4518.709861, 5322.835874, 6157.188251, 7009.203565, 7910.587422, 8810.047856, 9751.511668, 10739.355196}, 
	{665.739421, 1148.083349, 1714.541462, 2341.603048, 3023.189930, 3746.919933, 4525.868290, 5329.493149, 6163.296851, 7014.864920, 7915.749155, 8814.796544, 9755.765307, 10668.230508}, 
	{673.137069, 1154.731489, 1720.543482, 2347.089614, 3028.175778, 3751.507906, 4508.994826, 5306.299694, 6166.725091, 7017.972766, 7918.457008, 8818.796544, 9759.765307, 10672.230508}};
	
	static NekDouble EigenvaluesAnaMeshesRK4[6][14] =
	{{18.849560, 37.699110, 62.254810, 92.114370, 127.199300, 167.209300, 212.189900, 261.852900, 316.311600, 375.280200, 438.911800, 506.922500, 579.489600, 656.332100}, 
	{606.168101, 1187.787270, 1793.615712, 2391.865039, 2950.786843, 3487.596192, 4060.011000, 5074.723000, 6188.286000, 7399.972000, 8708.072000, 10111.160000, 11608.000000, 13197.470000}, 
	{684.081594, 1260.315774, 1857.923680, 2448.327633, 2999.987783, 3530.775456, 4073.368829, 4591.645387, 5121.018725, 5647.554500, 6330.380000, 7341.464000, 8418.710000, 9561.240000}, 
	{710.760453, 1285.192482, 1879.968339, 2467.708932, 3016.876152, 3545.638033, 4086.578922, 4603.479913, 5131.656118, 5657.223101, 6197.802683, 7022.172000, 8034.773000, 9107.511000}, 
	{724.129308, 1297.665239, 1890.995890, 2477.394447, 3025.267650, 3552.997300, 4093.052410, 4609.237151, 5136.747279, 5661.792457, 6397.687000, 7383.760000, 8431.678000, 9540.598000}, 
	{732.175780, 1305.179484, 1897.615616, 2483.199183, 3030.256925, 3557.347822, 4096.823380, 4612.552660, 5139.604517, 5664.300842, 6771.812000, 7388.760000, 8436.678000, 9546.598000}};
	
	class CFLtester : public SolverUtils::AdvectionSystem
    {
    public:
        friend class MemoryManager<CFLtester>;

        static EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
        {
            EquationSystemSharedPtr p = MemoryManager<CFLtester>
                ::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }

        static std::string className;

        virtual ~CFLtester();
        		
    protected:
        SolverUtils::RiemannSolverSharedPtr     m_riemannSolver;
        Array<OneD, Array<OneD, NekDouble> > m_velocity;
        Array<OneD, NekDouble>               m_traceVn;

        CFLtester(const LibUtilities::SessionReaderSharedPtr& pSession,
                  const SpatialDomains::MeshGraphSharedPtr& pGraph);

        void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                      Array<OneD,  Array<OneD, NekDouble> > &outarray,
                      const NekDouble time);

        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                          Array<OneD,  Array<OneD, NekDouble> > &outarray,
                          const NekDouble time);
		
        /// Get the normal velocity
        Array<OneD, NekDouble> &GetNormalVelocity();

        virtual void v_InitObject();

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        virtual void v_GenerateSummary(SummaryList& s);
    private:

		virtual NekDouble v_GetTimeStep(const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFL, NekDouble timeCFL);
		
		virtual NekDouble v_GetTimeStep(int ExpOrder, NekDouble CFL, NekDouble TimeStability);
        
        // Mapping of the real convective field on the standard element.
        // This function gives back the convective filed in the standard
        // element to calculate the stability region of the problem in a
        // unique way.
        Array<OneD, NekDouble> GetStdVelocity(
            const Array<OneD, Array<OneD,NekDouble> > inarray);
    };
}

#endif
