<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction with Periodic BCs P=5 (rotated periodic face) </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_CubePeriodic_RotateFace.xml</parameters>
    <files>
      <file description="Session File">Helmholtz3D_CubePeriodic_RotateFace.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8"> 0.000116486</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8"> 0.0453128 </value>
        </metric>
    </metrics>
</test>
