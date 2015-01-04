<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction with Periodic BCs P=5 </description>
    <executable>ADRSolver</executable>
    <parameters>Helmholtz3D_CubePeriodic_RotateFace.xml</parameters>
    <processes>5</processes>
    <files>
      <file description="Session File">Helmholtz3D_CubePeriodic_RotateFace.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="2e-5"> 0.000107671</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-4"> 0.0419256 </value>
        </metric>
    </metrics>
</test>
