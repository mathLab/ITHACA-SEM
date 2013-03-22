<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Unsteady Advection-Diffusion Order2 0.001</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion_Order2_0001.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion_Order2_0001.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.5458e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">3.64252e-08</value>
        </metric>
    </metrics>
</test>
