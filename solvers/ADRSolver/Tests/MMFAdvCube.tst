<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D MMF-Advection-Diffusion-Reaction</description>
    <executable>ADRSolver</executable>
    <parameters>MMFAdvCube.xml</parameters>
    <files>
        <file description="Session File">MMFAdvCube.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">2.97266e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0.000181408</value>
        </metric>
    </metrics>
</test>
