from NekPy.LibUtilities import SessionReader
from NekPy.StdRegions import ConstFactorMap, ConstFactorType
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

# Override polynomial order and create ContField2D.
graph.SetExpansionsToPolyOrder(10)
exp = ContField2D(session, graph, session.GetVariable(0))

# Construct factor map, using lambda from session file.
lamb = session.GetParameter("Lambda")
factors = ConstFactorMap()
factors[ConstFactorType.FactorLambda] = lamb

# Construct right hand side forcing term.
x, y = exp.GetCoords()
sol = np.sin(np.pi * x) * np.sin(np.pi * y)
fx = -(lamb + 2*np.pi*np.pi) * sol

# Solve Helmholtz equation.
helm_sol = exp.BwdTrans(exp.HelmSolve(fx, factors))
print("L infinity error: %.6e" % np.max(np.abs(helm_sol - sol)))

# Clean up!
session.Finalise()
