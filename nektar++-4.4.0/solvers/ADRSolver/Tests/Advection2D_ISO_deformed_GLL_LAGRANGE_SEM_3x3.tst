<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady FRSD advection GLL_LAGRANGE_SEM, P=10, homogeneous Dirichlet bcs, deformed elements</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_ISO_deformed_GLL_LAGRANGE_SEM_3x3.xml</parameters>
    <files>
        <file description="Session File">Advection2D_ISO_deformed_GLL_LAGRANGE_SEM_3x3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 0.0313807 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 0.129515 </value>
        </metric>
    </metrics>
</test>
