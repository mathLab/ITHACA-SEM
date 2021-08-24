<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction P=7</description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_nodal.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_nodal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12"> 1.62207e-08 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12"> 5.66566e-08 </value>
        </metric>
    </metrics>
</test>
