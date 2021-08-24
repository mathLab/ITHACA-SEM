
"""

The MIT License (MIT)

Copyright (c) 2018 ITHACA-SEM contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


main.py -- this is part of the ITHACA-SEM software framework

"""


import numpy as np
#from ITHACA_SEM_code import Generate_Matrices
#from ITHACA_SEM_code import Nektar_data
#from ITHACA_SEM_code import ROM_Oseen_Iteration
from ITHACA_SEM_code import *


#(Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, physglodof, LocGloSignA, LocGloMapA, LocGloSign, LocGloMap, snap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, tt_cd, c_f_bnd, c_f_p, c_f_int, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap) = Generate_Matrices.gen()

try:
	print("Looking for ITHACA data.")
	(Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloSign, LocGloMap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, LocGloMapMatA, init_bnd, forcing) = Generate_Matrices.gen()
	print("ITHACA data found.")
except Exception as e:
	print(e)
	try:
		print("No ITHACA data found, looking for Nektar data.")
		Nektar_data.convert()
		print("Nektar data converted to ITHACA - rerun to use ITHACA data.")
		raise SystemExit
	except Exception as e: 
		print(e)
		print("No Nektar data found, stopping.")
		raise SystemExit
	pass


print("going on...")

# double check the Dbnd and Dint -- va bene

f_snap1 = open('Jet_Generate_Nektar_Data_Oseen_p00255.txt', 'r')
#f_snap1 = open('cavity_poi_Oseen_ROM_1.txt', 'r')
f_snap2 = open('Jet_Generate_Nektar_Data_Oseen_p00725.txt', 'r')
#f_snap2 = open('cavity_poi_Oseen_ROM_2.txt', 'r')
snap1 = np.transpose(np.loadtxt(f_snap1))
snap2 = np.transpose(np.loadtxt(f_snap2))
snap = np.transpose(np.c_[snap1, snap2])
print("snap.shape ", snap.shape)

# also need to convert the snapshots here 
snap_x = snap[:, 0:ngc]
snap_y = snap[:, ngc:2*ngc]
#snap_x = snap[ 0:ngc]
#snap_y = snap[ ngc:2*ngc]

glodofphys_old = np.load('glodofphys.npy', 'r')
physglodof_old = np.load('physglodof.npy', 'r')

#print("snap_y.shape ", snap_y.shape)
#print("glodofphys.shape ", glodofphys.shape)
#print("glodofphys_old.shape ", glodofphys_old.shape)
#print("np.linalg.norm(glodofphys) ", np.linalg.norm(glodofphys))
#print("np.linalg.norm(glodofphys_old) ", np.linalg.norm(glodofphys_old))
#print("np.linalg.norm(glodofphys_old - glodofphys) ", np.linalg.norm(glodofphys_old - glodofphys))

nn = np.dot(snap_y, glodofphys)
print(nn[:,378])
nn = np.dot(snap_x, glodofphys)
print(nn[:,1296])
#print(nn[378])

#param_list = [  0.00255, 0.007125]
param_list = np.matrix([[ 1, 1], [20/20, 20/20]])
param_list = np.matrix([[ 0.00255, 0.007125], [1, 1]])
print(param_list)
for i in range(0, snap.shape[0]):
	curr_snap_x = snap_x[i,:]
	curr_snap_y = snap_y[i,:]
	curr_param = param_list[:,i]
	curr_snap_x_loc = np.dot(curr_snap_x, glodofphys)
	curr_snap_y_loc = np.dot(curr_snap_y, glodofphys)
	print("init qoi ", curr_snap_y_loc[378])
	print("at param ", curr_param)
	(u_x_new, u_y_new, f_bnd, f_p, f_int) = Oseen_step_trafo_p.Oseen_step_rom(curr_snap_x_loc, curr_snap_y_loc, curr_param,Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloSign, LocGloMap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, LocGloMapMatA, init_bnd, forcing)
	if i == 0:
		c_f_bnd = f_bnd
		c_f_p = f_p
		c_f_int = f_int
	else:
		c_f_bnd = np.c_[c_f_bnd, f_bnd]
		c_f_p = np.c_[c_f_p, f_p]
		c_f_int = np.c_[c_f_int, f_int]


# set the paramter vector of interest

#mu_vec = [Re, Gr]

mu_vec = [0.005, 1]

# what should this actually do?
# I can have a fixed-point algorithm class
# and a 'model' class specified as steady-state Navier-Stokes
# then 'attach' the Oseen to the cavity/channel model and do a solve, giving it the parameter of interest (Re, Gr)

# for this primo non-class approach:
# give the projection space in as well
# and do only evaluate for one parameter here, the sweep over the parameter domain can happen somewhere else
# should pre-compute the trilinear form as well and then
# append everything to the model (class type)

# where to collect Nektar data from?
# nektar_dev
# nektar_dev_ithaca  --> 
#	library/MultiRegions/AssemblyMap/AssemblyMapCG.cpp 
#	library/MultiRegions/GlobalLinSysDirectFull.cpp 
#	library/MultiRegions/GlobalLinSysDirect.cpp 
#	library/MultiRegions/ContField2D.h 
#	library/MultiRegions/ContField2D.cpp
#	library/MultiRegions/ExpList.cpp
#	library/LibUtilities/LinearAlgebra/NekLinSys.hpp
#	library/SolverUtils/UnsteadySystem.cpp
#	solvers/IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp
#	still some missing... -- trying to get it to compile at all by moving whole dirs
#	solvers/IncNavierStokesSolver/EquationSystems/
#	library/LibUtilities
#	library/MultiRegions
#	library/SolverUtils
#	library/SpatialDomains
#	library/StdRegions
#	library/LocalRegions
#	library/GlobalMapping
#	library/Collections
#	library/
#	utilities/
# nektar_thaca
# nektar_ithaca_visualize

# compile order takes LibUtilities first

(solution_x, solution_y) = ROM_Oseen_Iteration.Oseen_Iteration(mu_vec, Dbnd, Dint, sing_A, sing_B, sing_Btilde, sing_C, M_trafo_no_pp_incl_dbc, glodofphys, LocGloMapMatA, LocGloSign, LocGloMap, snap, ngc, nlc, bmap, imap, num_elem, nsize_bndry, nsize_int, nsize_p, nvel, c_f_bnd, c_f_p, c_f_int, NumDirBCs, nGlobHomBndDofs, no_loc_dbc, elem_loc_dbc, no_not_loc_dbc, elem_not_loc_dbc, IP, IPp, bwdtrans, cartmap0, cartmap1, LocGloBndSign, LocGloBndMap, BndCondCoeffsToGlobalCoeffsMap, init_bnd, forcing )



"""
[  2%] Built target zlib-1.2.7
[  3%] Built target boost
[  4%] Built target modmetis-5.1.0
[  5%] Built target tinyxml-2.6.2
[ 15%] Built target LibUtilities
[ 18%] Built target StdRegions
[ 22%] Built target SpatialDomains
[ 24%] Built target LocalRegions
[ 26%] Built target Collections
[ 32%] Built target MultiRegions
[ 33%] Built target GlobalMapping
[ 40%] Built target FieldUtils
[ 43%] Built target UnitTests
[ 49%] Built target SolverUtils
[ 52%] Built target NekMeshUtils
[ 54%] Built target LibUtilitiesUnitTests
[ 55%] Built target LocalRegionsUnitTests
[ 56%] Built target CollectionsUnitTests
[ 58%] Built target LinearAlgebraUnitTests
[ 58%] Built target FoundationDemo
[ 58%] Built target NodalDemo
[ 58%] Built target PartitionAnalyse
[ 58%] Built target TimeIntegrationDemo
[ 58%] Built target StdEquiToCoeff2D
[ 58%] Built target StdProject0D
[ 58%] Built target StdProject1D
[ 58%] Built target StdProject2D
[ 58%] Built target StdProject3D
[ 59%] Built target StdProject_Diff1D
[ 59%] Built target StdProject_Diff2D
[ 59%] Built target StdProject_Diff3D
[ 59%] Built target LocProject1D
[ 59%] Built target LocProject2D
[ 59%] Built target LocProject3D
[ 60%] Built target LocProject_Diff1D
[ 60%] Built target LocProject_Diff2D
[ 60%] Built target LocProject_Diff3D
[ 60%] Built target CollectionTiming
[ 60%] Built target Deriv3DHomo1D
[ 60%] Built target Deriv3DHomo1D_SingleMode
[ 60%] Built target Deriv3DHomo2D
[ 60%] Built target HDGHelmholtz1D
[ 60%] Built target HDGHelmholtz2D
[ 61%] Built target HDGHelmholtz3D
[ 61%] Built target HDGHelmholtz3DHomo1D
[ 61%] Built target Helmholtz1D
[ 61%] Built target Helmholtz2D
[ 61%] Built target Helmholtz3D
[ 61%] Built target Helmholtz3DHomo1D
[ 62%] Built target Helmholtz3DHomo2D
[ 62%] Built target PostProcHDG2D
[ 62%] Built target PostProcHDG3D
[ 62%] Built target SteadyAdvectionDiffusionReaction2D
[ 63%] Built target APESolver
[ 67%] Built target CardiacEPSolver
[ 69%] Built target PrePacing
[ 71%] Built target ShallowWaterSolver
[ 72%] Built target ADRSolver
[ 74%] Built target PulseWaveSolver
[ 76%] Built target Fld2Tecplot
[ 77%] Built target DiffusionSolver
[ 77%] Built target DiffusionSolverTimeInt
[ 78%] Built target LinearElasticSolver
[ 83%] Built target CompressibleFlowSolver
[ 83%] Built target CompressibleBL
[ 83%] Built target ExtractSurface2DCFS
[ 83%] Built target ExtractSurface3DCFS
[ 86%] Built target IncNavierStokesSolver
[ 86%] Built target AddModeTo2DFld
[ 88%] Built target Aliasing
[ 90%] Built target CFLStep
[ 90%] Built target ExtractMeanModeFromHomo1DFld
[ 90%] Built target Fld2DTo2D5
[ 90%] Built target FldAddFalknerSkanBL
[ 92%] Built target NonLinearEnergy
[ 92%] Built target FldAddFld
[ 92%] Built target FldAddMultiShear
[ 93%] Built target FldAddScalGrad
[ 93%] Built target FldAddScalGrad_elmt
[ 93%] Built target FldAddWSS
[ 93%] Built target MeshConvert
[ 93%] Built target ProbeFld
[ 93%] Built target XmlToTecplot
[ 94%] Built target XmlToTecplotWireFrame
[ 94%] Built target XmlToVtk
[ 94%] Built target RectangularMesh
[ 94%] Built target Regular2DMeshGenerator
[ 94%] Built target VariableValence2DMeshGenerator
[ 94%] Built target FieldConvert
[ 99%] Built target NekMesh
[100%] Built target Tester
"""




