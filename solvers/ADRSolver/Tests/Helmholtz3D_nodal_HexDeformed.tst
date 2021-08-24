<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction P=7 Hex Deformed</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_nodal_HexDeformed.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_nodal_HexDeformed.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 2.21825e-08 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 1.00387e-07 </value>
        </metric>
    </metrics>
</test>
