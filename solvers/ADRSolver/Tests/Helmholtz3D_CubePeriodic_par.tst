<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction with Periodic BCs P=3 </description>
    <executable>ADRSolver</executable>
    <parameters>--use-scotch Helmholtz3D_CubePeriodic.xml</parameters>
    <processes>5</processes>
    <files>
      <file description="Session File">Helmholtz3D_CubePeriodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-5"> 0.000129406</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="4e-3"> 0.0461794 </value>
        </metric>
    </metrics>
</test>
