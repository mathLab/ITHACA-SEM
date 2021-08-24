<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady FRDG advection GAUSS_LAGRANGE, P=3, Dirichlet bcs, deformed elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_dirichlet_deformed_GAUSS_LAGRANGE_10x10.xml</parameters>
    <files>
        <file description="Session File">Advection2D_dirichlet_deformed_GAUSS_LAGRANGE_10x10.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.000125624 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.000262082 </value>
        </metric>
    </metrics>
</test>
