from NekPy.LibUtilities import SessionReader, ReduceOperator
from NekPy.StdRegions import ConstFactorMap, ConstFactorType, VarCoeffMap, VarCoeffType
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ContField2D

import numpy as np
import sys

if len(sys.argv) < 2:
    print("Usage: python HelmSolve.py in.xml [in2.xml ...]")
    exit(1)

# Load up session and create ContField2D
session = SessionReader.CreateInstance(sys.argv)
graph = MeshGraph.Read(session)

# Grab the communicator from the SessionReader.
comm = session.GetComm()

# Override polynomial order and create ContField2D.
graph.SetExpansionsToPolyOrder(10)
exp = ContField2D(session, graph, session.GetVariable(0))

# Construct factor map, using lambda from session file.
lamb = session.GetParameter("Lambda")
factors = ConstFactorMap()
factors[ConstFactorType.FactorLambda] = lamb

# Test use of variable coefficients.
coeffs = VarCoeffMap()
coeffs[VarCoeffType.VarCoeffD00] = np.ones(exp.GetNpoints())
coeffs[VarCoeffType.VarCoeffD11] = np.ones(exp.GetNpoints())

# Construct right hand side forcing term.
x, y = exp.GetCoords()
sol = np.sin(np.pi * x) * np.sin(np.pi * y)
fx = -(lamb + 2*np.pi*np.pi) * sol

# Solve Helmholtz equation.
helm_sol = exp.BwdTrans(exp.HelmSolve(fx, factors, coeffs))
L2_error = exp.L2(helm_sol, sol)
Linf_error = exp.Linf(helm_sol, sol)
Linf_error_comm = comm.AllReduce(
    np.max(np.abs(helm_sol - sol)), ReduceOperator.ReduceMax)

# Test reduction of Array types.
reduce_test = np.zeros(comm.GetSize())
reduce_test[comm.GetRank()] = 1.0
comm.AllReduce(reduce_test, ReduceOperator.ReduceSum)

# Print out some stats for debugging.
if comm.GetRank() == 0:
    print("L 2 error (variable nek)     : %.6e" % L2_error)
    print("L inf error (variable nek)   : %.6e" % Linf_error)
    print("L inf error (variable nekpy) : %.6e" % Linf_error_comm)
    print("Reduction test               : %d" % round(reduce_test.sum()))

# Clean up!
session.Finalise()
