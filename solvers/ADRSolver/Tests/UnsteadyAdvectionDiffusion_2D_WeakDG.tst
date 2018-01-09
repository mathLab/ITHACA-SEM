<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Advection-Diffusion with WeakDG (epsilon = 0.5)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_2D_WeakDG.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_2D_WeakDG.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">3.087e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.46935e-07</value>
        </metric>
    </metrics>
</test>
