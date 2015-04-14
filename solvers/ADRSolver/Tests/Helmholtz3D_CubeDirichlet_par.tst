<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 3D Helmholtz/Steady Diffusion Reaction with Dirichlet  BCs P=5 </description>
    <executable>ADRSolver</executable>
    <parameters>--use-metis Helmholtz3D_CubeDir.xml</parameters>
    <processes>12</processes>
    <files>
      <file description="Session File">Helmholtz3D_CubeDir.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8"> 4.82538e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8"> 0.010592 </value>
        </metric>
    </metrics>
</test>
