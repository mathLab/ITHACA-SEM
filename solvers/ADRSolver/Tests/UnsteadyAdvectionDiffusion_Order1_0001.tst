<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Unsteady Advection-Diffusion Order1 0.001</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_Order1_0001.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_Order1_0001.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0.000135563</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0.000276572</value>
        </metric>
    </metrics>
</test>
