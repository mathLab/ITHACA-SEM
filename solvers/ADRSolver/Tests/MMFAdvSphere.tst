<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D MMF-Advection-Diffusion-Reaction</description>
    <executable>ADRSolver</executable>
    <parameters>MMFAdvSphere.xml</parameters>
    <files>
        <file description="Session File">MMFAdvSphere.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">0.0209316</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-06">0.204615</value>
        </metric>
    </metrics>
</test>
