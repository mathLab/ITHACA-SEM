<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Advection-Diffusion with Imex-DG</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_2D_ImexDG.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_2D_ImexDG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">7.67001e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">6.63362e-05</value>
        </metric>
    </metrics>
</test>
