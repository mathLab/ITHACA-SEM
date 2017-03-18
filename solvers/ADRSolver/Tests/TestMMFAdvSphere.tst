<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D MMF-Advection-Diffusion-Reaction</description>
    <executable>ADRSolver</executable>
    <parameters>TestMMFAdvSphere.xml</parameters>
    <files>
        <file description="Session File">TestMMFAdvSphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">0.0125831</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-06">0.0521416</value>
        </metric>
    </metrics>
</test>
