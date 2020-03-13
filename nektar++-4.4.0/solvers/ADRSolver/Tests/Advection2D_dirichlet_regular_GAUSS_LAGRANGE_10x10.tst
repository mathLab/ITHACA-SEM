<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady WeakDG advection GAUSS_LAGRANGE, P=3, Q=5 Dirichlet bcs, regular elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_dirichlet_regular_GAUSS_LAGRANGE_10x10.xml</parameters>
    <files>
        <file description="Session File">Advection2D_dirichlet_regular_GAUSS_LAGRANGE_10x10.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 5.34741e-05 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 8.66536e-05 </value>
        </metric>
    </metrics>
</test>
