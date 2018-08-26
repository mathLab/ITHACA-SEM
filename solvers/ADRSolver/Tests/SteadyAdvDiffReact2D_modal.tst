<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Advection-Diffusion-Reaction</description>
    <executable>ADRSolver</executable>
    <parameters>SteadyAdvDiffReact2D_modal.xml</parameters>
    <files>
        <file description="Session File">SteadyAdvDiffReact2D_modal.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0.00873717</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0.0046866</value>
        </metric>
    </metrics>
</test>
