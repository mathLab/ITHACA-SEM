<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction P=5 All Element Types </description>
    <executable>ADRSolver</executable>
    <parameters>CubeAllElements.xml</parameters>
    <files>
        <file description="Session File">CubeAllElements.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6"> 0.000902596 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6"> 0.024387 </value>
        </metric>
    </metrics>
</test>
