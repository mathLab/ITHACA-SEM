<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady FRSD advection GAUSS_LAGRANGE, P=10, homogeneous Dirichlet bcs, regular elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_ISO_regular_GAUSS_LAGRANGE_3x3.xml</parameters>
    <files>
        <file description="Session File">Advection2D_ISO_regular_GAUSS_LAGRANGE_3x3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.00154764 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.00960366 </value>
        </metric>
    </metrics>
</test>
